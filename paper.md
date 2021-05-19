---
title: 'beamshapes: directionality functions for various sound source models'
tags:
  - Python
  - acoustics
  - bioacoustics
authors:
  - name: Thejasvi Beleyur
    orcid: 0000-0001-5360-4383
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)

affiliations:
 - name: Centre for the Advanced Study of Collective Behaviour, University of Konstanz, Konstanz
   index: 1
 - name: Department of Biology, University of Konstanz, Konstanz
   index: 2
 - name: Acoustic and Functional Ecology, Max Planck Institute for Ornithology, Seewiesen
   index: 3
date: May 2021
bibliography: paper.bib
---

# Summary

* Directivity functions ($D_{\theta}/D_{0}$) describe the relative sound energy levels at different angular positions, relative to the on-axis level ($0^{\circ}$). 
* The directivity functions of various sound sources are described in literature [@beranek2012acoustics;beranek2019acoustics] but there aren't any openly available computational implementations 
* ```beamshapes``` implements directivity functions of multiple sound source models in one place.
* Implementing multiple models together allows an easy comparison of model and measurements, or parameter estimation [guarato2011method].


# Statement of need

* Until recently, computational power has been a limiting factor to calculating directivity functions for various models. Models with easily calculable outputs have thus been favoured (eg. piston in an infinite baffle), especially in the field of bioacoustics [mogensen1979sound;@surlykke2009echolocating]. 
* Now, with increasing computational power in hand, calculating and implementing directivity functions is much easier. However, despite the increase in computational power, older, simpler models with limited biological relevance continue to dominate the field of bio-acoustics. For instance, the piston in an infinite baffle only predicts the beam-shape for a $\pm90^{\circ}$ range off-axis, and assumes front-back symmetry. This is unrealistic for most vocalising animals, especially echolocating bats and odontocetes. However, the model continues to be a standard reference point [@jakobsen2013convergent;@macaulay2020high]. 
* `beamshapes` allows for the 


# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements
TB thanks Gaurav Dhariwal for his continual math advice and inputs, and Tim Mellow for providing Mathematica code to assist with porting to Python.

# References
