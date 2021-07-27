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
 - name: Acoustic and Functional Ecology, Max Planck Institute for Ornithology, Seewiesen
   index: 2
 - name: Max Planck Institute for Animal Behavior, Radolfzell
   index: 3
date: May 2021
bibliography: paper.bib
---

# Summary

Sound sources such as human beings or loudspeaker exhibit a 'directionality'. A
listener or microphone placed at different angles at a fixed radius will often 
pick up sometimes drastically different sound levels. The directionality of 
a sound source is often modelled as a combination of the frequency of the emitted
sound, the geometry of the sound source itself, and the combination of the physics
of sound propagation over space. The resulting pattern of sound levels with angular location
 is called the *directivity* function [@beranek2012acoustics]. 

Directivity functions describe the relative sound levels at a given angle off-axis (at angle $\theta$)
with relation to the level on-axis (at $\theta = 0$), $directivity = \frac{D_{\theta}}{D_{0}}$. Directivity
functions exist for a wide variety of sound sources that can be modelled analytically. A well known example 
of a directivity function is that of the piston in an infinite baffle. The piston is a circular surface of radius
$a$, vibrating back and forth about a hole in an infinite wall or baffle. The directivity is described by $\frac{2J_{1}(ka \times sin \theta)}{ka \times sin \theta}$, 
where $J_{1}$ is the Bessel's function of the first kind, $k$ is the wavenumber, where $k = \frac{2\pi f}{c}$, where $f$ is the frequency 
of the sound, and $c$ is the speed of sound. The angle of the receiver is $\theta$, which varies from from 0-2$\pi$ radians in the azimuth. 

Directivity functions can be used  experimentally for both prediction (to estimate the expected directionality of
a given sound source), and inference (to estimate the parameters which match an observed set of sound level measurements).





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
TB thanks Gaurav Dhariwal for his continual math advice and inputs, and Tim Mellow for 
providing Mathematica code and clarifications to assist with porting models to Python.
This work was funded by a mix of TB's private funds, the IMPRS for Organismal Biology, 
a DAAD stipend, and then finally by a CASCB Medium grant.

# References
