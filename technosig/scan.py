"""Sweep a folder of TESS light curves for dimming events and catalogue them.

The approach borrows from field archaeology: record every feature first
(every dip in every star), describe it consistently, then work out which
processes explain it. Instrument systematics are the post-depositional
disturbance; eclipsing binaries and planets are the "natural" explanations.
Whatever is left after that is the interesting residue, which is usually
very little.

Usage:
    python -m technosig.scan [DATA_DIR] [--limit N] [--out results]
"""

import argparse
import glob
import os
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares

DATA_DIR = os.environ.get("TECHNOSIG_DATA", "/run/media/dylan/ELEMENTS/Lightcurves")

BIN_DAYS = 0.5 / 24       # bin 2-min cadence to 30 min
TREND_DAYS = 1.5          # running-median window; dips much longer than ~0.5 d get eaten
DIP_SIGMA = 3.0           # event core threshold
EDGE_SIGMA = 1.0          # event boundaries extend until residual recovers to this
MIN_BINS = 2              # i.e. at least 1 hour below DIP_SIGMA
GAP_DAYS = 0.5            # events this close to a data gap/edge are flagged


def load(path):
    """Read one SPOC light curve; keep only good-quality, finite cadences."""
    with fits.open(path, memmap=False) as h:
        h0, h1, d = h[0].header, h[1].header, h[1].data
        ok = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"]) & np.isfinite(d["TIME"])
        med = np.nanmedian(d["PDCSAP_FLUX"][ok])
        lc = {
            "time": np.array(d["TIME"][ok], dtype=float),
            "flux": np.array(d["PDCSAP_FLUX"][ok], dtype=float) / med,
            "bkg": np.array(d["SAP_BKG"][ok], dtype=float),
            "cen1": np.array(d["MOM_CENTR1"][ok], dtype=float),
            "cen2": np.array(d["MOM_CENTR2"][ok], dtype=float),
        }
        meta = {
            "tic": h0.get("TICID"),
            "sector": h0.get("SECTOR"),
            "camera": h0.get("CAMERA"),
            "ccd": h0.get("CCD"),
            "ra": h0.get("RA_OBJ"),
            "dec": h0.get("DEC_OBJ"),
            "tmag": h0.get("TESSMAG"),
            "teff": h0.get("TEFF"),
            "radius": h0.get("RADIUS"),
            "crowdsap": h1.get("CROWDSAP"),
            "cdpp1h_ppm": h1.get("CDPP1_0"),
            "n_good": int(ok.sum()),
            "file": os.path.basename(path),
        }
    return lc, meta


def bin_uniform(t, cols, width=BIN_DAYS):
    """Average onto a uniform time grid; empty bins become NaN so gaps stay visible."""
    edges = np.arange(t.min(), t.max() + width, width)
    idx = np.digitize(t, edges) - 1
    n = len(edges) - 1
    counts = np.bincount(idx, minlength=n)[:n]
    out = {"time": edges[:-1] + width / 2}
    with np.errstate(invalid="ignore", divide="ignore"):
        for k, v in cols.items():
            good = np.isfinite(v)
            s = np.bincount(idx[good], weights=v[good], minlength=n)[:n]
            c = np.bincount(idx[good], minlength=n)[:n]
            out[k] = np.where(c > 0, s / c, np.nan)
    out["counts"] = counts
    return out


def robust_sigma(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


def find_events(t, resid, sigma):
    """Contiguous runs below -DIP_SIGMA, widened to where flux recovers past -EDGE_SIGMA."""
    core = resid < -DIP_SIGMA * sigma
    shallow = resid < -EDGE_SIGMA * sigma
    events, i, n = [], 0, len(resid)
    while i < n:
        if not core[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and core[j + 1]:
            j += 1
        if j - i + 1 >= MIN_BINS:
            a, b = i, j
            while a > 0 and shallow[a - 1]:
                a -= 1
            while b + 1 < n and shallow[b + 1]:
                b += 1
            if events and a <= events[-1][1] + 1:
                events[-1] = (events[-1][0], max(b, events[-1][1]))
            else:
                events.append((a, b))
        i = j + 1
    return events


def describe_event(b, a, z, sigma, gap_mask):
    """Measure one event the same way every time: depth, duration, shape, context."""
    seg = slice(a, z + 1)
    r = b["resid"][seg]
    t = b["time"][seg]
    k = int(np.nanargmin(r))
    t_min = t[k]
    t0, t1 = t[0] - BIN_DAYS / 2, t[-1] + BIN_DAYS / 2
    dur = t1 - t0
    ingress, egress = t_min - t0, t1 - t_min
    # how far the background and centroid move during the dip, in units of their own scatter
    def excursion(key):
        v = b[key]
        s = robust_sigma(v)
        if not np.isfinite(s) or s == 0:
            return np.nan
        return float((np.nanmean(v[seg]) - np.nanmedian(v)) / s)

    return {
        "t_min": t_min,
        "t_start": t0,
        "t_end": t1,
        "duration_h": dur * 24,
        "depth_ppm": float(-r[k] * 1e6),
        "snr": float(-r[k] / sigma),
        "integrated_snr": float(-np.nansum(r) / (sigma * np.sqrt(np.isfinite(r).sum()))),
        "asymmetry": float((egress - ingress) / dur),  # >0: fast drop, slow recovery
        "n_bins": int(np.isfinite(r).sum()),
        "near_gap": bool(gap_mask[max(a - 1, 0)] or gap_mask[min(z + 1, len(gap_mask) - 1)]),
        "bkg_excursion": excursion("bkg"),
        "cen1_excursion": excursion("cen1"),
        "cen2_excursion": excursion("cen2"),
    }


def _bls_best(t, y, periods):
    res = BoxLeastSquares(t, y).power(periods, [0.04, 0.08, 0.16, 0.25])
    i = int(np.argmax(res.power))
    P, t0, dur = res.period[i], res.transit_time[i], res.duration[i]
    n_in = int((np.abs(((t - t0) / P + 0.5) % 1 - 0.5) * P < dur / 2).sum())
    snr = res.depth[i] / (robust_sigma(y) / np.sqrt(max(n_in, 1)))
    return P, t0, dur, res.depth[i], snr


def periodic_search(t, resid):
    """BLS on the residuals, plus the same search on the inverted residuals.

    Correlated noise produces BLS "detections" in both directions; a real periodic
    dimming is much stronger as a dip than anything the inverted light curve offers.
    """
    m = np.isfinite(resid)
    if m.sum() < 100:
        return {}
    t, y = t[m], resid[m]
    periods = np.exp(np.linspace(np.log(0.3), np.log(13.0), 4000))
    P, t0, dur, depth, snr = _bls_best(t, y, periods)
    *_, inv_snr = _bls_best(t, -y, periods)
    return {
        "bls_period": float(P),
        "bls_t0": float(t0),
        "bls_duration_h": float(dur * 24),
        "bls_depth_ppm": float(depth * 1e6),
        "bls_snr": float(snr),
        "bls_inverted_snr": float(inv_snr),
    }


def scan_one(path):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lc, meta = load(path)
            if meta["n_good"] < 500:
                return meta, []
            b = bin_uniform(lc["time"], {k: lc[k] for k in ("flux", "bkg", "cen1", "cen2")})
            f = pd.Series(b["flux"])
            win = int(round(TREND_DAYS / BIN_DAYS)) | 1
            trend = f.rolling(win, center=True, min_periods=win // 3).median().to_numpy()
            b["resid"] = b["flux"] - trend
            sigma = robust_sigma(b["resid"])
            # a bin is "gap" if it, or anything within GAP_DAYS of it, has no data
            empty = ~np.isfinite(b["flux"])
            g = int(round(GAP_DAYS / BIN_DAYS))
            gap_mask = np.convolve(np.r_[np.ones(g), empty, np.ones(g)], np.ones(2 * g + 1), "same")[g:-g] > 0

            meta.update({
                "noise_ppm_30min": sigma * 1e6,
                "trend_rms_ppm": float(np.nanstd(trend) * 1e6),  # rotation/pulsation strength
                # how far the star rises above its trend; dips no deeper than this are just flicker
                "upward_p995_ppm": float(np.nanpercentile(b["resid"], 99.5) * 1e6),
            })
            meta.update(periodic_search(b["time"], b["resid"]))

            events = []
            for a, z in find_events(b["time"], b["resid"], sigma):
                e = describe_event(b, a, z, sigma, gap_mask)
                if "bls_period" in meta:
                    ph = ((e["t_min"] - meta["bls_t0"]) / meta["bls_period"] + 0.5) % 1 - 0.5
                    e["bls_phase_h"] = ph * meta["bls_period"] * 24
                e["tic"] = meta["tic"]
                events.append(e)
            covered = b["time"][np.isfinite(b["flux"])]
            P, m = event_comb([e["t_min"] for e in events], [e["duration_h"] for e in events], covered)
            meta["comb_period"] = P
            for e, hit in zip(events, m):
                e["comb_match"] = bool(hit)
            meta["n_events"] = len(events)
            return meta, events
    except Exception as exc:  # one bad file shouldn't stop the sweep
        return {"file": os.path.basename(path), "error": repr(exc)}, []


def event_comb(t_min, dur_h, covered):
    """Do a star's events line up on a regular ephemeris (allowing secondaries at phase 0.5)?

    `covered` is the array of time bins that actually have data. A trial ephemeris
    must account for most of the primaries it predicts inside that coverage, which
    stops fractions of the true period (P/2, P/3...) matching by coincidence.
    Returns (period, matched mask). Catches eclipsing binaries whose eclipses are
    too long or too unequal for BLS to lock onto.
    """
    t = np.asarray(t_min)
    n = len(t)
    best, best_key = (np.nan, np.zeros(n, bool)), (0, 0.0)
    if n < 3:
        return best
    tol = max(0.1, 0.5 * np.median(dur_h) / 24)
    lo, hi = covered.min(), covered.max()
    diffs = np.abs(t[:, None] - t[None, :])[np.triu_indices(n, 1)]
    cands = {round(d / k, 4) for d in diffs for k in (1, 2, 3, 4) if d / k >= 0.3}
    for P in cands:
        for t0 in t:
            ph = ((t - t0) / P) % 1
            prim = np.minimum(ph, 1 - ph) * P < tol
            if prim.sum() < 2:
                continue
            epochs = t0 + P * np.arange(np.ceil((lo - t0) / P), np.floor((hi - t0) / P) + 1)
            seen = np.abs(covered[None, :] - epochs[:, None]).min(axis=1) < BIN_DAYS
            if prim.sum() < 0.7 * seen.sum():
                continue
            m = prim | (np.abs(ph - 0.5) * P < tol)
            key = (m.sum(), P)  # most events explained, then the longest such period
            if m.sum() >= 3 and key > best_key:
                best, best_key = (P, m), key
    return best


def classify(stars, events):
    """Assign each event to the most likely explanation.

    Precedence runs from "the star itself explains it" (periodic, variable) through
    "the telescope explains it" (common-mode, gaps, background, centroid). Whatever no
    rule claims is left as unexplained.
    """
    if events.empty:
        return events
    events = events.copy()

    # Periodic: either the event comb or a strong BLS signal the event sits on.
    comb = events["comb_match"].fillna(False).to_numpy(bool)
    s = stars.set_index("tic")
    ev_star = s.reindex(events["tic"])
    bls_hit = (
        (ev_star["bls_snr"].to_numpy() > 10)
        & (ev_star["bls_snr"].to_numpy() > 2 * ev_star["bls_inverted_snr"].to_numpy())
        & (np.abs(events["bls_phase_h"].to_numpy()) < np.maximum(ev_star["bls_duration_h"].to_numpy(), 1.0))
    )
    periodic = comb | bls_hit

    # Variable star: the dip is no bigger than the star's own slow wandering, or the
    # star rises above its trend about as far as it dips below it (symmetric flicker).
    depth = events["depth_ppm"].to_numpy()
    variability = (depth < 2 * ev_star["trend_rms_ppm"].to_numpy()) | (
        depth < 1.5 * ev_star["upward_p995_ppm"].to_numpy()
    )

    # Common-mode: far more stars dipping in one quarter-day slot than the sector-wide rate.
    slot = np.floor(events["t_min"] / 0.25)
    per_slot = events.groupby(slot)["tic"].nunique()
    base = per_slot.median()
    events["simultaneous_stars"] = slot.map(per_slot).to_numpy()
    common = (events["simultaneous_stars"] > base + 4 * np.sqrt(base) + 1).to_numpy()

    centroid = (events[["cen1_excursion", "cen2_excursion"]].abs().max(axis=1) > 3).to_numpy()
    background = (events["bkg_excursion"] > 3).to_numpy()
    gap = events["near_gap"].to_numpy(bool)

    cls = np.full(len(events), "unexplained", dtype=object)
    claimed = np.zeros(len(events), bool)
    for mask, label in [
        (periodic, "periodic"),
        (variability, "variability"),
        (common, "common_mode"),
        (gap, "gap_edge"),
        (background, "background"),
        (centroid, "centroid_shift"),
    ]:
        take = mask & ~claimed
        cls[take] = label
        claimed |= take
    events["class"] = cls
    return events


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("data_dir", nargs="?", default=DATA_DIR)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--out", default="results")
    p.add_argument("--workers", type=int, default=os.cpu_count())
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.data_dir, "*_lc.fits")))[: args.limit]
    print(f"Scanning {len(files)} light curves from {args.data_dir}")
    with Pool(args.workers) as pool:
        out = pool.map(scan_one, files, chunksize=8)

    stars = pd.DataFrame([m for m, _ in out])
    events = pd.DataFrame([e for _, ev in out for e in ev])
    events = classify(stars, events)

    os.makedirs(args.out, exist_ok=True)
    stars.to_csv(os.path.join(args.out, "stars.csv"), index=False)
    events.to_csv(os.path.join(args.out, "events.csv"), index=False)

    print(f"{len(stars)} stars, {int(stars.get('error', pd.Series()).notna().sum())} errors, {len(events)} events")
    if not events.empty:
        print(events["class"].value_counts().to_string())


if __name__ == "__main__":
    main()
