# Paleo-Detectors

*Is it possible to make a Dark Matter detection with a paleo-detector?*

[![DOI](https://zenodo.org/badge/142072044.svg)](https://zenodo.org/badge/latestdoi/142072044) [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)


Code for calculating track spectra in paleo-detectors and for reproducing the results of arXiv:1811.XXXXX.

More information about paleo-detectors can be found in [arXiv:1806.05991](http://arxiv.org/abs/1806.05991) and [arXiv:1811.06844](http://arxiv.org/abs/1811.06844).

**Authors:** Tom D P Edwards, Bradley J Kavanagh.

Please let get in touch with questions, comments or bug-reports.

### Paleopy



### Notebooks

*  [`Notebooks/PlotSpectra.ipynb`](Notebooks/PlotSpectra.ipynb) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tedwards2412/paleo_detectors/master?filepath=Notebooks%2FPlotSpectra.ipynb) - testing the paleopy classes (loading in minerals and displaying data), then plotting the relevant track length spectra. Used for generating Fig. 1 of arXiv:1811.XXXXX.
* [`Notebooks/FullAnalysis.ipynb`](Notebooks/FullAnalysis.ipynb) - generating SwordFish objects, then calculating and plotting Information Flux, Upper Limits and Discovery reach for different materials. This way, we only have to initialise the SwordFish objects once. Used for generating Fig. 2 of arXiv:1811.XXXXX.
* [`Notebooks/Systematics_Check.ipynb`](Notebooks/Systematics_Check.ipynb) - Used for generating Figs. 3 and 4 of arXiv:1811.XXXXX.
* [`Notebooks/Mass_Reconstruction.ipynb`](Notebooks/Mass_Reconstruction.ipynb) -  Used for generating Fig. 5 of arXiv:1811.XXXXX.


### SRIM data

SRIM data can be found in [`Data/dRdESRIM`](Data/dRdESRIM). There, you'll also find [`CleanSRIM.ipynb`](Data/dRdESRIM/CleanSRIM.ipynb), which lets you clean up the standard SRIM output files so that they can be read by the code (you just have to trim off the top and bottom junk from the files first, then the notebook can read them in and format them properly).

### Requirements

The code in this repo should run with Python3. Standard library requirements are in [`requirements.txt`](requirements.txt). In addition, you will need [`swordfish`](https://github.com/cweniger/swordfish) for performing the statistical analysis and [`WIMpy`](https://github.com/bradkav/WIMpy_NREFT) for calculating the DM and neutrino spectra. 

### Results

Tables of projected upper limits and discovery reach are output to [`ES/limits`](ES/limits).


