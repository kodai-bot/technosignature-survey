# Field notes: TESS Sector 59 technosignature sweep

**Sample:** 1,095 of the 12,996 two-minute light curves in Sector 59 (BTJD 2911–2937, late 2022), downloaded Feb 2023. One file is truncated and was skipped.

**Approach:** first record every dimming event in every star, and describe each one the same way. Then try to explain each event, starting with the most mundane explanation. Only what's left after that counts as residue worth attention.

## Pass 1: short dips, under about half a day (`technosig/scan.py`)

The scan found 1,328 events. Each was assigned the first explanation that fits:

| Explanation | Events |
|---|---|
| periodic (eclipsing binary or planet) | 791 |
| stellar variability | 389 |
| common-mode (many stars at once) | 78 |
| next to a data gap | 44 |
| background spike | 3 |
| **unexplained** | **23 events, 18 stars** |

**Checked against the TESS planet-candidate (TOI) catalogue:** 23 TOIs are in the sample. The confirmed planets XO-3 b and HAT-P-58 b were recovered as periodic. Four shortlist stars turned out to be known candidates, including TIC 3224605 = TOI-5552.01, a 28.6-day orbit caught as a single transit.

**Small planets:** TOI-1685 b, a known small planet, was missed. The sweep isn't designed for small planets.

## Follow-up: TIC 85423257, then HD 21476, then BD+67 274

1. **TIC 85423257 had a single 0.26% V-shaped dip lasting 10 h, and it wasn't in any catalogue.** The pixel difference image showed the signal wasn't on the target. It came from a much brighter star about 100″ away. (`followup/TIC85423257_difference_image.png`)
2. **That brighter star is HD 21476** (TIC 85368942, F-type, G = 7.6). It is a Gaia double-lined spectroscopic binary with a 1.295-day orbit and doesn't eclipse. In Sector 59 it shows a single 3.1% dip, about 14.5 h long and symmetric. Sectors 19, 73 and 86 show nothing comparable. (`followup/HD21476_all_sectors.png`, `followup/HD21476_s59_dip.png`)
3. **The dip is really on BD+67 274** (TIC 85368955, G = 9.15), 22″ from HD 21476, which is about one TESS pixel. Both a PSF fit and a model-free test place it there: the centre of the missing light sits 0.1–0.26 px from BD+67 274 and 0.8–0.9 px from HD 21476, with ±0.01 px uncertainty from bootstrapping. That implies a ~26% eclipse. (`../followup/locate_dip_source.py`, `followup/BD+67_274_dip_location.png`)
4. **BD+67 274 is also a Gaia double-lined spectroscopic binary.** Its orbit is P = 56.39 d, e = 0.30, with two near-identical stars (mass ratio 1.01).
   - **What fits:** a central eclipse in that orbit should last about 14–16 h, and we measured 14.5 h.
   - **What doesn't fit:** the Gaia orbit puts our dip 42° away from conjunction, where an eclipse should fall. That orbit's fit quality is poor (goodness 32).
   - No catalogue we checked lists BD+67 274 as an eclipsing binary: VSX, SIMBAD, Gaia's eclipsing-binary list, the TOI list.

**Interpretation:** a natural eclipse in a long-period F-star binary that, as far as these catalogues show, nobody has reported. It is not a technosignature. The timing mismatch with Gaia's orbit is still open: either Gaia's orbit is off, or a third body is involved.

**HD 21476's Gaia RUWE of 10.9** points to something more going on in that system as well. That is a separate question.

## Pass 2: long dips, half a day to several days (`technosig/longdips.py`)

This pass looks for Boyajian's-star-style dimming. The baseline is an upper envelope over 5 days, and each dip must stand out from both the point-to-point noise and the star's own variability. It runs on both PDCSAP flux and SAP flux corrected by a CCD ensemble.

- 265 long dips were found. After removing dips next to data gaps and dips shared across a CCD, **47 dips in 33 stars** remained.
- Looked at one by one, most SAP-only dips come from curved slow drifts that the simple CCD-ensemble correction doesn't remove. Most PDCSAP dips are slow variable stars.
- **No Boyajian-like candidates were found in this sample.**
- **Instrument note:** long dips bunch up around BTJD 2919–2920 across several cameras (23% of dips fall in 7% of the time), which points to an instrumental event there. The BD+67 274 eclipse is also at 2919.3. The pixel-level location on a single star argues against a shared systematic, but this should be kept in mind.

## Limitations

- **Only 8% of one sector was searched.**
- **PDCSAP processing may already have removed slow astrophysical dimming.** The SAP route needs proper cotrending basis vectors (the standard instrument-systematics correction) before it can be trusted.
- The thresholds were tuned by eye on this sample and haven't been tested against injected signals.
