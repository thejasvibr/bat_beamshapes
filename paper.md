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

Sound sources such as human beings or loudspeaker often exhibit a 'directionality' in how loud they sound at different angles.
A listener or microphone placed at various angles at a fixed radius will often pick up sometimes drastically different sound levels. Oftentimes the same sound source may actually produce omnidirectional sound fields, with rather uniform levels at various angles. The directionality of a sound source can be modelled as a combination of the frequency of the emitted
sound, the geometry of the vibrating and non-vibrating parts of the sound source itself. The resulting pattern of sound radiation with angular location is called the *directivity* of the source [@beranek2012acoustics]. 

Directivity ($D$) describes the relative sound levels at a given angle $D(\theta)$
with relation to the level on-axis ($D(0)$, and thus $directivity = \frac{D_{\theta}}{D_{0}}$). Directivity
functions exist for a wide variety of sound sources that can be modelled analytically. A well known example 
of a directivity function is that of the piston in an infinite baffle. The piston is a circular surface of radius
$a$, vibrating back and forth about a hole in an infinite wall or baffle. The directivity is described by $\frac{2J_{1}(ka \times sin \theta)}{ka \times sin \theta}$, where $J_{1}$ is the Bessel's function of the first kind, $k$ is the wavenumber, where $k = \frac{2\pi f}{c}$, where $f$ is the frequency of the sound, and $c$ is the speed of sound. The angle of the receiver is $\theta$, which varies from from 0-2$\pi$ radians in the azimuth. 

Directivity functions can be used in a two-fold manner 1) to engineer deliberately engineer devices to suit particular specifications, eg. loudspeaker sound fields [@beranek2012acoustics] and 2) to infer parameters of a sound source itself having assumed a relevant model, eg. estimating direction of call emission [@guarato2011method] and mouth aperture of bat echolocation calls [@jakobsen2013convergent].


# Statement of need

A host of published directivity functions exist in the literature, but it is my experience that their computational implementations remain as inhouse scripts that are often in proprietary language platforms. To my knowledge there are no computational implementations of sound source directivities that are open-source and developed using modern software practices such as version-control and unit-testing. In this paper I present ```beamshapes```, a Python package that implements directivity patterns for sound source models. As of this publication, ```beamshapes``` version 0.2.0 implements four sound sources 1) piston in an infinite baffle, 2) point source of a sphere 3) oscillating cap of a sphere 4) piston in a sphere. 

Computational implementations of directivity functions often require long run-times due to the intensive numerical routines and arbitrary precision math required to generate results. Long run-times hinder scientific projects in reducing the number of models and parameter space that can be explored. ```beamshapes``` boasts parallelised code to generate significant speed-ups in run-times. 

The availability of openly-available directivity implementations will hopefully stir the acoustics, and specifically the bio-acoustics community to test and compare sound radiation using alternate models that may present a better description of the data. Until recently, computational power has perhaps been a limiting factor to calculating directivity functions for various models. Models with easily calculable outputs have thus been favoured (eg. piston in an infinite baffle), especially in the field of bioacoustics [@strother1970acoustical;@mogensen1979sound]. However, despite the recent availability of computational power, older, simpler models with limited biological relevance continue to dominate the field of bio-acoustics. For instance, the piston in an infinite baffle only predicts the beam-shape for a $\pm90^{\circ}$ range off-axis, and assumes front-back symmetry. This is unrealistic for most vocalising animals, especially echolocating bats and odontocetes. However, the model continues to be a standard reference point for multiple studies ranging over the past decade [@jakobsen2013convergent;@kounitsky2015bats;@macaulay2020high]. The piston in a sphere, implemented in ```beamshapes```, for instance recreates many of the highly directional central and side lobes and that are seen in bats, while also predicting sound radiation behind the bat (Figure \autoref{fig:pistonsphereinf}). 


![Comparison on directivity patterns for piston in a sphere (left) and piston in an infinite baffle (right) for various ka values. The half-angle aperture of the piston ($\alpha$) of the sphere was set to 30 degrees. There is a broad similarity in the primary lobe, though the differences become clear with increasing off-axis angles. The lack of front-back symmetry may be very informative for parameter estimation of real bats in the field for instance. \label{fig:pistonsphereinf}](paper_related/piston_sphere_baffle.png)


Future releases of ```beamshapes``` are scheduled to include directivity patterns for additional models of interest such as rectangular cap of a sphere or piston in a closed finite baffle.


# Software packages used in this work
```beamshapes``` relies on the Python open-source ecosystem and is built on the numpy, scipy, sympy, mpmath and flint libraries. 


# Figures

Figures can be included like this:

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements
TB thanks Gaurav Dhariwal for his continual math advice and inputs, and Tim Mellow for 
providing Mathematica code and clarifications to assist with porting models to Python.
This work was funded by a mix of TB's private funds, the IMPRS for Organismal Biology, 
a DAAD stipend, and then finally by a CASCB Medium grant.

# References
