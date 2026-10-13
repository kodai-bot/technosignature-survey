"""Which star is the HD 21476 dip on? Fit Gaia-positioned PSFs to the out-of-dip and difference images."""
import glob, sys, numpy as np, lightkurve as lk, astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from scipy.optimize import least_squares
# usage: python followup/locate_dip_source.py WORKDIR  (WORKDIR/cut/ holds the TESScut FITS; plot is written there)
SP = sys.argv[1]
tpf = lk.read(glob.glob(f'{SP}/cut/**/*.fits', recursive=True)[0])
t = np.asarray(tpf.time.value, float); F = np.asarray(tpf.flux.value, float)
good = np.isfinite(F).all(axis=(1, 2))
t, F = t[good], F[good]
inn = (t > 2919.18) & (t < 2919.45)                    # bottom of the dip
out = ((t > 2918.4) & (t < 2918.9)) | ((t > 2919.8) & (t < 2920.3))
mo, mi = np.median(F[out], 0), np.median(F[inn], 0)
diff = mo - mi
err = np.sqrt(np.var(F[out], 0) / out.sum() + np.var(F[inn], 0) / inn.sum())
print('in', inn.sum(), 'out', out.sum())

c = SkyCoord(53.006313 * u.deg, 67.585513 * u.deg)
g = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 'Gmag'], row_limit=-1).query_region(
    c, radius=7.5 * 21 * u.arcsec, catalog='I/355/gaiadr3')[0].to_pandas()
g = g[g.Gmag < 13.5].reset_index(drop=True)
gc = SkyCoord(g.RA_ICRS.values * u.deg, g.DE_ICRS.values * u.deg)
g['px'], g['py'] = tpf.wcs.world_to_pixel(gc)
g['sep'] = c.separation(gc).arcsec
g = g[(g.px > -1) & (g.px < 15) & (g.py > -1) & (g.py < 15)].sort_values('Gmag').reset_index(drop=True)
print(g.round(2).to_string())
yy, xx = np.indices(mo.shape)

def psf(x0, y0, s):
    return np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s ** 2)) / (2 * np.pi * s ** 2)

def design(dx, dy, s):
    return np.stack([psf(x + dx, y + dy, s).ravel() for x, y in zip(g.px, g.py)] + [np.ones(mo.size)], 1)

# fit global offset + PSF width on the out-of-dip image, amplitudes linear
def resid(p):
    A = design(*p)
    amp, *_ = np.linalg.lstsq(A, mo.ravel(), rcond=None)
    return (A @ amp - mo.ravel()) / np.sqrt(np.abs(mo.ravel()) + 100)
fit = least_squares(resid, [0, 0, 0.9], bounds=([-1, -1, 0.4], [1, 1, 2.5]))
dx, dy, s = fit.x
A = design(dx, dy, s)
amp_out, *_ = np.linalg.lstsq(A, mo.ravel(), rcond=None)
print(f'offset dx {dx:.2f} dy {dy:.2f} psf sigma {s:.2f}')

# difference image: which star's PSF carries the missing flux?
w = 1 / err.ravel()
amp_d, *_ = np.linalg.lstsq(A * w[:, None], diff.ravel() * w, rcond=None)
cov = np.linalg.inv((A * w[:, None]).T @ (A * w[:, None]))
se = np.sqrt(np.diag(cov))
chi_all = np.sum(((A @ amp_d - diff.ravel()) * w) ** 2)
print('\nper-star fit of the difference image (only stars near HD 21476 shown)')
for k in range(len(g)):
    if g.sep[k] < 60:
        print(f"  G={g.Gmag[k]:5.2f} sep {g.sep[k]:5.1f}\"  missing flux {amp_d[k]:8.0f} ± {se[k]:5.0f}  -> fractional dip {amp_d[k]/amp_out[k]*100:6.2f}%")

# one-source hypotheses: all the missing flux on star k alone
print('\nsingle-source hypotheses (lower chi2 = better), dof =', diff.size - 2)
for k in range(len(g)):
    if g.sep[k] < 60:
        B = np.stack([A[:, k], A[:, -1]], 1)
        a, *_ = np.linalg.lstsq(B * w[:, None], diff.ravel() * w, rcond=None)
        chi = np.sum(((B @ a - diff.ravel()) * w) ** 2)
        print(f"  dip on G={g.Gmag[k]:5.2f} ({g.sep[k]:4.1f}\"): chi2 {chi:9.1f}  implied dip depth {a[0]/amp_out[k]*100:6.2f}%")

# Model-free check: light-weighted centre of the missing light vs. of the normal image
def centroid(img, x0, y0, r=2):
    m = (np.abs(xx - x0) <= r) & (np.abs(yy - y0) <= r)
    v = np.where(m, img, 0) - np.median(img)
    v = np.where(m, np.clip(v, 0, None), 0)
    return (v * xx).sum() / v.sum(), (v * yy).sum() / v.sum()
hx, hy = g.px[0] + dx, g.py[0] + dy
nx, ny = g.px[1] + dx, g.py[1] + dy
print(f'\nHD 21476 (G7.6) at ({hx:.2f},{hy:.2f}); neighbour (G9.2) at ({nx:.2f},{ny:.2f}); separation {np.hypot(hx-nx, hy-ny):.2f} px')
for r in (1, 2, 3):
    co = centroid(mo, hx, hy, r); cd = centroid(diff, hx, hy, r)
    print(f'box ±{r}px: normal image centre ({co[0]:.2f},{co[1]:.2f})  missing-light centre ({cd[0]:.2f},{cd[1]:.2f})'
          f'  dist to HD21476 {np.hypot(cd[0]-hx, cd[1]-hy):.2f}px, to neighbour {np.hypot(cd[0]-nx, cd[1]-ny):.2f}px')
# bootstrap the missing-light centre by resampling cadences
rng = np.random.default_rng(1)
ii, oo = np.where(inn)[0], np.where(out)[0]
bs = []
for _ in range(300):
    d = np.median(F[rng.choice(oo, oo.size)], 0) - np.median(F[rng.choice(ii, ii.size)], 0)
    bs.append(centroid(d, hx, hy, 2))
bs = np.array(bs)
print('bootstrap missing-light centre (±2px box): x %.2f ± %.2f, y %.2f ± %.2f' % (bs[:, 0].mean(), bs[:, 0].std(), bs[:, 1].mean(), bs[:, 1].std()))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
for a, img, ttl in [(ax[0], np.log10(np.clip(mo, 1, None)), 'out-of-dip image (log)'), (ax[1], diff, 'missing light (out - in)')]:
    a.imshow(img, origin='lower', cmap='viridis')
    a.plot(hx, hy, 'r+', ms=14, mew=2, label='HD 21476 (G 7.6)'); a.plot(nx, ny, 'wx', ms=12, mew=2, label='neighbour (G 9.2)')
    a.set_title(ttl); a.set_xlim(hx - 4, hx + 4); a.set_ylim(hy - 4, hy + 4)
ax[1].plot(bs[:, 0].mean(), bs[:, 1].mean(), 'o', mfc='none', mec='orange', ms=14, mew=2, label='missing-light centre')
ax[1].legend(fontsize=8, loc='upper right')
fig.tight_layout(); fig.savefig(f'{SP}/hd21476_locate.png', dpi=75)
