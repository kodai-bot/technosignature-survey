"""Diagnostic plots: one star's light curve with its catalogued events marked."""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .scan import BIN_DAYS, DATA_DIR, TREND_DAYS, bin_uniform, load

CLASS_COLOURS = {
    "unexplained": "crimson",
    "periodic": "tab:blue",
    "common_mode": "0.6",
    "gap_edge": "tab:olive",
    "background": "tab:orange",
    "centroid_shift": "tab:purple",
    "variability": "tab:green",
}


def plot_star(file, events=None, ax=None, data_dir=DATA_DIR, title=None):
    lc, meta = load(os.path.join(data_dir, file))
    b = bin_uniform(lc["time"], {k: lc[k] for k in ("flux", "bkg")})
    win = int(round(TREND_DAYS / BIN_DAYS)) | 1
    trend = pd.Series(b["flux"]).rolling(win, center=True, min_periods=win // 3).median()

    if ax is None:
        fig, ax = plt.subplots(2, 1, figsize=(12, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    a0, a1 = ax
    a0.plot(lc["time"], lc["flux"], ",", color="0.75")
    a0.plot(b["time"], b["flux"], ".", ms=2, color="k")
    a0.plot(b["time"], trend, "-", lw=0.8, color="tab:green")
    a1.plot(b["time"], b["bkg"], "-", lw=0.8, color="tab:orange")
    a1.set_ylabel("bkg")
    if events is not None:
        for _, e in events.iterrows():
            c = CLASS_COLOURS.get(e["class"], "k")
            for a in (a0, a1):
                a.axvspan(e["t_start"], e["t_end"], color=c, alpha=0.25, lw=0)
    lo, hi = np.nanpercentile(b["flux"], [0.5, 99.5])
    pad = (hi - lo) * 0.3
    a0.set_ylim(lo - pad, hi + pad)
    a0.set_ylabel("norm. PDCSAP flux")
    a1.set_xlabel("time (BTJD)")
    a0.set_title(title or f"TIC {meta['tic']}  Tmag {meta['tmag']:.1f}  cam {meta['camera']}/{meta['ccd']}", fontsize=9)
    return meta
