# Technosignature survey: an archaeological approach to stellar light curves

This started in February 2023 as a short MSc exercise. It was my one attempt at "SETI as archaeology" during college: using TESS and Kepler light curves to look for dimming that natural processes don't explain, in the spirit of Boyajian's star (KIC 8462852). The repository (formerly `Planetary_Transit_Detection`) now continues that work as a small survey. It treats stellar light curves the way a site is treated, recording every feature first and then explaining it away layer by layer.

**The original 2023 work is preserved at the git tag `msc-2023`** and as notebooks in [`msc-2023/`](msc-2023/).

## Layout

| Path | What it is |
|---|---|
| [`docs/technosignature_search_criteria.md`](docs/technosignature_search_criteria.md) | The design brief: search criteria derived from an archaeological, physical definition of technology |
| [`docs/research_program.md`](docs/research_program.md) | How those criteria become work packages for this instrument |
| [`technosig/scan.py`](technosig/scan.py) | Pass 1: catalogue every short dip (under ~0.5 day) in every star, then classify each by its most mundane explanation |
| [`technosig/longdips.py`](technosig/longdips.py) | Pass 2: long, Boyajian-type dimming (0.5 day to several days) |
| [`technosig/plots.py`](technosig/plots.py) | Diagnostic plots |
| [`followup/locate_dip_source.py`](followup/locate_dip_source.py) | Works out which star in the pixels a dip really belongs to |
| [`results/FIELD_NOTES.md`](results/FIELD_NOTES.md) | Findings and limitations from the first sweep |
| [`msc-2023/`](msc-2023/) | The original notebooks, unchanged |

## Data

The first sweep used 1,095 TESS Sector 59 two-minute light curves (2.0 GB). They are kept outside the repository. Set `TECHNOSIG_DATA` to the folder that holds them (the code falls back to the path on my own machine).

FITS files are git-ignored. MAST provides bulk-download scripts for each sector at <https://archive.stsci.edu/tess/bulk_downloads.html>.

## Running

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv -r requirements.txt

.venv/bin/python -m technosig.scan        # writes results/stars.csv, results/events.csv
.venv/bin/python -m technosig.longdips    # writes results/longdip_stars.csv, results/longdips.csv
```

Each pass takes about 1–2 minutes on 4 cores for the 1,095-star sample.

## First results (Sector 59, 1,095 stars)

- **No technosignature candidates.** Not surprising in a sample this size.
- **The method checks out against known objects.** Confirmed planets came out as periodic. TOI-5552.01 was picked up from a single transit.
- **One unexplained dip, followed to its source.** The dip seen in TIC 85423257 came from the bright F star HD 21476, and then from its neighbour BD+67 274, only 22″ away. BD+67 274 is a Gaia double-lined spectroscopic binary (P = 56.4 d). The dip is a ~26% eclipse there. None of the catalogues checked list it as an eclipsing binary.
- **No Boyajian-like long dips.**

Details and limitations are in [`results/FIELD_NOTES.md`](results/FIELD_NOTES.md).

## Status

This is early-stage work. The next work packages are listed in [`docs/research_program.md`](docs/research_program.md):

- a reference library of natural slow dimming processes
- injection–recovery completeness tests
- target lists built from archaeological priors: solar siblings, stars that passed close to the Sun, and the Earth transit zone

## Acknowledgements

The 2023 MSc notebooks were learned from STScI's MAST tutorial notebooks for TESS and Kepler (then at `spacetelescope/notebooks`, now [spacetelescope/mast_notebooks](https://github.com/spacetelescope/mast_notebooks)) and from the [Lightkurve tutorials](https://lightkurve.github.io/lightkurve/tutorials/index.html).

This work uses data collected by the TESS mission, obtained from the Mikulski Archive for Space Telescopes (MAST) at the Space Telescope Science Institute. Funding for the TESS mission is provided by NASA's Science Mission Directorate.

## Licence and citation

Code: MIT, see [LICENSE](LICENSE). The written material in `docs/` and `results/` may be reused with attribution. To cite this work, see [CITATION.cff](CITATION.cff).
