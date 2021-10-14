# Skipper-CCD image analysis
A set of python (and jupyter) files with different image analysis for my thesis.
The .ipynb are all the work I've been doing so far, with all steps explained with some level of detail.

#### Analisis_Imagenes.ipynb:
General analysis of the images. First steps in understanding the data and characteristics of the CCD sensor. It contains a lot of messy information and chunks of code because it had a lot of trial and error.

#### Analisis_Imagenes_probabilidades.ipynb
Two ways for computing the probability of having one electron events in a given pixel. The main idea is, suppossing a Poissonian distribution for the random variable "number of events per pixel", computing the expected value "mu" and then calculate the probability.

#### Analisis_Imagenes_clusters.ipynb
General analysis for clusters found in the images. Comparisons between clusters found and cluster with syntetic noise added, etc

#### Analisis_Imagenes_bordes.ipynb
Analysis of the border of the clusters, where there could be one electron events due to simply diffusion of the charge (genuine charge) or because of spurious events. Comparison between the real number of events in borders vs the expected number of events in borders.

#### Analisis_Imagenes_descomposicion_eventos.ipynb
Of course not all one electron events are due to spurious charges (e.g. dark currents), some are genuine events. So far, all analyses were being performed using the total effect of one electron events The idea of this notebook is to separate the contribution of spurious events from the total effect, and then compute the contribution of the real one electron events. That way, one could compute the expected number of real one electron events.

#### analyfits
Module with the functions used (first use in Analisis_Imagenes_descomposicion_eventos.ipynb), for code ordering purposes.
##### calibration.py
Python script with the current calibration to transform ADU's to electrons
##### datafits.py
Functions to extract data from .fits images
##### manipufits.py
Functions to manipulate with python .fits images
##### misc.py
Functions for diverse purposes, such as searching .fits files in a directory tree or thresholding any image.
