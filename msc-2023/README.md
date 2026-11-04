# Original MSc work (February 2023)

These notebooks are the exercise as submitted, unchanged. The same state is kept at the git tag
`msc-2023`. They are kept for the record; the current survey code is in [`../technosig/`](../technosig/).

| Notebook | What it is |
|---|---|
| `dylanLightKurve.ipynb` | My own exercise: target-pixel and light-curve analysis of Boyajian's star (KIC 8462852), the starting point for this project |
| `WASP-96.ipynb` | Worked through from a TESS science tutorial: recreating the TESS transit of WASP-96 b |
| `ExamplePipeline2.ipynb` | Worked through from a MAST/TESS tutorial on reading light-curve FITS files |
| `deepTransit1.ipynb` | A trial of the pre-trained Deep Transit detection model; it fails because the model file path is the original author's |
| `Dubhe.ipynb` | Reading a FITS image with WCS coordinates; its input file (`92991.fits`) is no longer available |

Known problems, left as they are: some cells need data files that aren't in the repository, and
masked lightkurve arrays need `.filled(np.nan)` before plotting with current versions.
