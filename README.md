# paleo_detectors
*Is it possible to make a DM detection with a paleo detector?*

### Notebooks

*  [`Notebooks/PlotSpectra.ipynb`](Notebooks/PlotSpectra.ipynb) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tedwards2412/paleo_detectors/master?filepath=Notebooks%2FPlotSpectra.ipynb) - testing the paleopy classes (loading in minerals and displaying data), the plotting the relevant track length spectra. 
* [`Notebooks/FullAnalysis.ipynb`](Notebooks/FullAnalysis.ipynb) - generating SwordFish objects, then calculating and plotting Information Flux, Upper Limits and Discovery reach for different materials. This way, we only have to initialise the SwordFish objects once.

### SRIM data

SRIM data can be found in [`Data/dRdESRIM`](Data/dRdESRIM). There, you'll also find [`CleanSRIM.ipynb`](Data/dRdESRIM/CleanSRIM.ipynb), which lets you clean up the standard SRIM output files so that they can be read by the code (you just have to trim off the top and bottom junk from the files first, then the notebook can read them in and format them properly).

### Results

Tables of projected upper limits and discovery reach are output to [`ES/limits`](ES/limits).
