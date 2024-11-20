# Representational dissimilarity component analysis (ReDisCA)

This is an official implementation of the paper Alexei Ossadtchi, Ilia Semenkov, Anna Zhuravleva, 
Vladimir Kozunov, Oleg Serikov, Ekaterina Voloshina, Representational dissimilarity component analysis (ReDisCA), 
NeuroImage, Volume 301, 2024, 120868, ISSN 1053-8119, https://doi.org/10.1016/j.neuroimage.2024.120868.

This repository contains MATLAB code needed to replicate the results of the paper.

# Prerequisites

## Data

Before you can run scripts, please download [MEG data](https://osf.io/8rk67/) and put it into the ```data``` directory. 
If you are using our data please cite our paper.

## SPoC

The code in this repository relies on the SPoC algorithm. Therefore, before running scripts, please install SPoC 
[MATLAB implementation](https://github.com/bbci/bbci_public/tree/master/external) by Berlin Brain-Computer Interface

# Code 

The repository consists of 2 main MATLAB scripts:

* ```Redisca_tools_faces_3_random_norm_correct.m``` - main script, should be run first, performs ReDisCA algorithm
* ```Redisca_source_loc_for_tools_faces_3_random_.m``` - script that performs source localication

# Citation

If you use our algorithm, code or data please cite our paper:

```
@article{OSSADTCHI2024120868,
title = {Representational dissimilarity component analysis (ReDisCA)},
journal = {NeuroImage},
volume = {301},
pages = {120868},
year = {2024},
issn = {1053-8119},
doi = {https://doi.org/10.1016/j.neuroimage.2024.120868},
url = {https://www.sciencedirect.com/science/article/pii/S1053811924003653},
author = {Alexei Ossadtchi and Ilia Semenkov and Anna Zhuravleva and Vladimir Kozunov and Oleg Serikov and Ekaterina Voloshina},
keywords = {EEG and MEG, Spatial–temporal decomposition, Representational similarity analysis, Source localization},
abstract = {The principle of Representational Similarity Analysis (RSA) posits that neural representations reflect the structure of encoded information, allowing exploration of spatial and temporal organization of brain information processing. Traditional RSA when applied to EEG or MEG data faces challenges in accessing activation time series at the brain source level due to modeling complexities and insufficient geometric/anatomical data. To overcome this, we introduce Representational Dissimilarity Component Analysis (ReDisCA), a method for estimating spatial–temporal components in EEG or MEG responses aligned with a target representational dissimilarity matrix (RDM). ReDisCA yields informative spatial filters and associated topographies, offering insights into the location of ”representationally relevant” sources. Applied to evoked response time series, ReDisCA produces temporal source activation profiles with the desired RDM. Importantly, while ReDisCA does not require inverse modeling its output is consistent with EEG and MEG observation equation and can be used as an input to rigorous source localization procedures. Demonstrating ReDisCA’s efficacy through simulations and comparison with conventional methods, we show superior source localization accuracy and apply the method to real EEG and MEG datasets, revealing physiologically plausible representational structures without inverse modeling. ReDisCA adds to the family of inverse modeling free methods such as independent component analysis (Makeig, 1995), Spatial spectral decomposition (Nikulin, 2011), and Source power comodulation (Dähne, 2014) designed for extraction sources with desired properties from EEG or MEG data. Extending its utility beyond EEG and MEG analysis, ReDisCA is likely to find application in fMRI data analysis and exploration of representational structures emerging in multilayered artificial neural networks.}
}
```
