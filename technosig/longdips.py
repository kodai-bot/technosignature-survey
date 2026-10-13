"""Second pass: look for long, irregular dimming of the Boyajian's-star kind.

The first pass (scan.py) smooths each light curve over 1.5 days, which erases any
dip longer than about half a day. Tabby-style dips last from half a day to weeks,
so this pass works differently:

* It uses the raw SAP flux as well as PDCSAP, because the PDC pipeline can
  remove slow astrophysical signals along with the systematics.
* Each CCD is its own "reference trench". The median behaviour of every star on
  the same CCD is the shared instrumental signal (scattered light, pointing drift),
  and dividing it out leaves what is particular to each star.
* The baseline is an upper envelope (a rolling high percentile), because dimming
  events only ever go down. A symmetric variable star rises as far as it falls; a
  dip-dominated one doesn't.

Usage:
    python -m technosig.longdips [DATA_DIR] [--out results]
"""

import argparse
import glob
import os
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from astropy.io import fits

from .scan import DATA_DIR, robust_sigma

BIN_DAYS = 2 / 24          # 2-hour bins
ENVELOPE_DAYS = 5.0        # rolling window for the upper-envelope baseline
ENVELOPE_PCTL = 80         # which percentile counts as "the star's normal brightness"
DIP_SIGMA = 3.0            # against point-to-point noise
VAR_SIGMA = 4.0            # against the star's own variability
MIN_DUR_DAYS = 0.5         # only interested in dips at least 12 h long
MIN_STARS_PER_CCD = 15
EDGE_DAYS = 0.75           # dips this close to a data gap or edge are flagged


def load_binned(path, t_grid):
    """Read one light curve and average SAP and PDCSAP onto a shared 2-h grid."""
    with fits.open(path, memmap=False) as h:
        h0, d = h[0].header, h[1].data
        ok = (d["QUALITY"] == 0) & np.isfinite(d["TIME"])
        t = np.array(d["TIME"][ok], float)
        out = {"tic": h0.get("TICID"), "camera": h0.get("CAMERA"), "ccd": h0.get("CCD"),
               "tmag": h0.get("TESSMAG"), "file": os.path.basename(path)}
        idx = np.digitize(t, t_grid) - 1
        n = len(t_grid) - 1
        for col in ("SAP_FLUX", "PDCSAP_FLUX"):
            v = np.array(d[col][ok], float)
            good = np.isfinite(v) & (idx >= 0) & (idx < n)
            s = np.bincount(idx[good], weights=v[good], minlength=n)[:n]
            c = np.bincount(idx[good], minlength=n)[:n]
            with np.errstate(invalid="ignore", divide="ignore"):
                f = np.where(c >= 20, s / c, np.nan)  # need most of the 60 cadences per bin
            out[col.lower()] = f / np.nanmedian(f)
    return out


def _safe_load(args):
    path, t_grid = args
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return load_binned(path, t_grid)
    except Exception as exc:
        return {"file": os.path.basename(path), "error": repr(exc)}


def upper_envelope(f, window_bins):
    s = pd.Series(f)
    return s.rolling(window_bins, center=True, min_periods=window_bins // 3).quantile(ENVELOPE_PCTL / 100).to_numpy()


def find_long_dips(t, f):
    """Dips below the upper envelope lasting at least MIN_DUR_DAYS."""
    win = int(round(ENVELOPE_DAYS / BIN_DAYS)) | 1
    base = upper_envelope(f, win)
    r = f / base - 1
    # point-to-point scatter: not inflated by slow trends or by the dips themselves
    sigma = robust_sigma(np.diff(r[np.isfinite(r)])) / np.sqrt(2)
    # ...and the star's own wandering about its envelope (rotation, pulsation). A dip
    # has to stand out from both, otherwise every slow variable is "dipping".
    sigma_var = robust_sigma(r)
    thresh = max(DIP_SIGMA * sigma, VAR_SIGMA * sigma_var)
    below = r < -thresh
    dips, i, n = [], 0, len(r)
    while i < n:
        if not below[i]:
            i += 1
            continue
        j = i
        # allow single-bin interruptions (noise) inside a long dip
        while j + 1 < n and (below[j + 1] or (j + 2 < n and below[j + 2])):
            j += 1
        dur = (j - i + 1) * BIN_DAYS
        if dur >= MIN_DUR_DAYS:
            seg = r[i : j + 1]
            k = int(np.nanargmin(seg))
            # local minima count: >1 means a structured, multi-component dip
            sm = pd.Series(seg).rolling(3, center=True, min_periods=1).mean().to_numpy()
            mins = int(np.sum((sm[1:-1] < sm[:-2]) & (sm[1:-1] < sm[2:]) & (sm[1:-1] < -thresh)))
            e = int(round(EDGE_DAYS / BIN_DAYS))
            near = np.r_[f[max(i - e, 0) : i], f[j + 1 : j + 1 + e]]
            near_edge = bool(i - e < 0 or j + e >= n or (~np.isfinite(near)).any())
            dips.append({
                "t_start": t[i] - BIN_DAYS / 2,
                "t_end": t[j] + BIN_DAYS / 2,
                "t_min": t[i + k],
                "duration_d": dur,
                "depth_pct": float(-seg[k] * 100),
                "snr": float(-seg[k] / sigma),
                "asymmetry": float(((j - (i + k)) - k) / max(j - i, 1)),
                "n_minima": mins,
                "near_edge": near_edge,
            })
        i = j + 1
    # dip dominance: how far it falls below normal vs how far it rises above
    lo, med, hi = np.nanpercentile(r, [1, 50, 99])
    dominance = (med - lo) / max(hi - med, 1e-9)
    return dips, sigma, sigma_var, dominance


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("data_dir", nargs="?", default=DATA_DIR)
    p.add_argument("--out", default="results")
    p.add_argument("--workers", type=int, default=os.cpu_count())
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.data_dir, "*_lc.fits")))
    with fits.open(files[0]) as h:
        t = h[1].data["TIME"]
        t_grid = np.arange(np.nanmin(t), np.nanmax(t) + BIN_DAYS, BIN_DAYS)
    t_mid = t_grid[:-1] + BIN_DAYS / 2

    print(f"Loading {len(files)} light curves onto a {BIN_DAYS * 24:.0f}-h grid")
    with Pool(args.workers) as pool:
        rows = pool.map(_safe_load, [(f, t_grid) for f in files], chunksize=8)
    rows = [r for r in rows if "error" not in r]

    # CCD ensembles: the shared instrumental behaviour of every star on that CCD
    ens = {}
    for key in {(r["camera"], r["ccd"]) for r in rows}:
        members = [r for r in rows if (r["camera"], r["ccd"]) == key]
        if len(members) >= MIN_STARS_PER_CCD:
            ens[key] = np.nanmedian(np.vstack([m["sap_flux"] for m in members]), axis=0)

    stars, dips = [], []
    for r in rows:
        rec = {k: r[k] for k in ("tic", "camera", "ccd", "tmag", "file")}
        variants = {"pdcsap": r["pdcsap_flux"]}
        key = (r["camera"], r["ccd"])
        if key in ens:
            variants["sap_ens"] = r["sap_flux"] / ens[key]
        for name, f in variants.items():
            if np.isfinite(f).sum() < 100:
                continue
            ds, sigma, sigma_var, dom = find_long_dips(t_mid, f)
            rec[f"{name}_noise_pct"] = sigma * 100
            rec[f"{name}_var_pct"] = sigma_var * 100
            rec[f"{name}_dominance"] = dom
            rec[f"{name}_n_dips"] = len(ds)
            for d in ds:
                d.update(tic=r["tic"], flux=name, camera=r["camera"], ccd=r["ccd"])
                dips.append(d)
        stars.append(rec)

    stars = pd.DataFrame(stars)
    dips = pd.DataFrame(dips)
    if not dips.empty:
        # common-mode on the same CCD even after the ensemble: overlapping dips in many stars
        def overlaps(row):
            same = dips[(dips.camera == row.camera) & (dips.ccd == row.ccd) & (dips.flux == row.flux)]
            hit = (same.t_start < row.t_end) & (same.t_end > row.t_start)
            return same.loc[hit, "tic"].nunique()
        dips["ccd_overlap_stars"] = dips.apply(overlaps, axis=1)
        dips["in_both"] = False
        for tic, g in dips.groupby("tic"):
            a, b = g[g.flux == "pdcsap"], g[g.flux == "sap_ens"]
            for i, row in g.iterrows():
                other = b if row.flux == "pdcsap" else a
                dips.loc[i, "in_both"] = bool(((other.t_start < row.t_end) & (other.t_end > row.t_start)).any())

    os.makedirs(args.out, exist_ok=True)
    stars.to_csv(os.path.join(args.out, "longdip_stars.csv"), index=False)
    dips.to_csv(os.path.join(args.out, "longdips.csv"), index=False)
    np.save(os.path.join(args.out, "ccd_ensembles.npy"), {"t": t_mid, **{f"{k[0]}-{k[1]}": v for k, v in ens.items()}})
    print(f"{len(stars)} stars, {len(ens)} CCD ensembles, {len(dips)} long dips")
    if not dips.empty:
        print(dips.groupby("flux").size().to_string())


if __name__ == "__main__":
    main()
