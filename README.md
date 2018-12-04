# Paleo-Detectors

*Is it possible to make a Dark Matter detection with a paleo-detector?*

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tedwards2412/paleo_detectors/master?filepath=Notebooks%2FPlotSpectra.ipynb) [![DOI](https://zenodo.org/badge/142072044.svg)](https://zenodo.org/badge/latestdoi/142072044)  [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

Code for calculating track length spectra in paleo-detectors and exploring paleo-detector sensitivity to Dark Matter (DM). This code can be used to reproduce the results of [arXiv:1811.10549](http://arxiv.org/abs/1811.10549), "*Digging for Dark Matter: Spectral Analysis and Discovery Potential of Paleo-Detectors*".

More information about paleo-detectors can also be found in [arXiv:1806.05991](http://arxiv.org/abs/1806.05991) and [arXiv:1811.06844](http://arxiv.org/abs/1811.06844).

**Authors:** Thomas D P Edwards, Bradley J Kavanagh.

Please get in touch with any questions, comments or bug-reports.

### Overview: Paleopy

The core of the code is in the [`paleopy`](paleopy) module. Data for converting recoil energies to track lengths, along with tables of background distributions are loaded from [`Data/`](Data). This then allows you to calculate all the relevant track length distributions. The currently supported minerals are Nchwaningite, Sinjarite, Halite, Olivine, Gypsum and Phlogopite (see [`Data/MineralList.txt`](Data/MineralList.txt)).

To install run:

    git clone https://github.com/tedwards2412/paleo_detectors
    cd paleo_detectors
    python3 setup.py install

Check out [`Notebooks/PlotSpectra.ipynb`](Notebooks/PlotSpectra.ipynb) for an illustration of how to use the code: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tedwards2412/paleo_detectors/master?filepath=Notebooks%2FPlotSpectra.ipynb)


### Notebooks

*  [`Notebooks/PlotSpectra.ipynb`](Notebooks/PlotSpectra.ipynb) - testing the paleopy code (loading in minerals and displaying data), then plotting the relevant track length spectra. Used for generating Fig. 1 of [arXiv:1811.10549](http://arxiv.org/abs/1811.10549).
* [`Notebooks/FullAnalysis.ipynb`](Notebooks/FullAnalysis.ipynb) - generating SwordFish objects, then calculating and plotting Information Flux, Upper Limits and Discovery reach for different materials. This way, we only have to initialise the SwordFish objects once. Used for generating Fig. 2 of [arXiv:1811.10549](http://arxiv.org/abs/1811.10549).
* [`Notebooks/Systematics_Check.ipynb`](Notebooks/Systematics_Check.ipynb) - calculating upper limits for a range of different background normalisation and shape systematic uncertainties. Used for generating Figs. 3 and 4 of [arXiv:1811.10549](http://arxiv.org/abs/1811.10549).
* [`Notebooks/Mass_Reconstruction.ipynb`](Notebooks/Mass_Reconstruction.ipynb) -  generating Euclideanised signals and projected DM mass reconstruction contours. Used for generating Fig. 5 of [arXiv:1811.10549](http://arxiv.org/abs/1811.10549).


### SRIM data

SRIM data can be found in [`Data/dRdESRIM`](Data/dRdESRIM). There, you'll also find [`CleanSRIM.ipynb`](Data/dRdESRIM/CleanSRIM.ipynb), which lets you clean up the standard SRIM output files so that they can be read by the code (you just have to trim off the top and bottom junk from the files first, then the notebook can read them in and format them properly).

### Results

Tables of projected upper limits and discovery reach are output to [`ES/limits`](ES/limits).

### Requirements

The code in this repo should run with Python3. Standard library requirements are in [`requirements.txt`](requirements.txt). In addition, you will need [`swordfish`](https://github.com/cweniger/swordfish) for performing the statistical analysis and [`WIMpy`](https://github.com/bradkav/WIMpy_NREFT) for calculating the DM and neutrino spectra. 
