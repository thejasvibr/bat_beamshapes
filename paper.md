---
title: '```beamshapes```: a Python package to generate directivity patterns for various sound source models'
authors:
- affiliation: 1, 2, 3
  name: Thejasvi Beleyur
  orcid: 0000-0001-5360-4383
date: "May 2021"
output: pdf_document
bibliography: paper.bib
tags:
- Python
- acoustics
- bioacoustics
affiliations:
- index: 1
  name: Centre for the Advanced Study of Collective Behaviour, University of Konstanz,
    Konstanz
- index: 2
  name: Acoustic and Functional Ecology, Max Planck Institute for Ornithology, Seewiesen
- index: 3
  name: Max Planck Institute for Animal Behavior, Radolfzell
---

# Summary

Sound sources such as human beings or loudspeakers often exhibit a 'directionality' in how loud they sound at different angles.
A listener or microphone placed at various angles at a fixed radius will often pick up sometimes drastically different sound levels. The same sound source may however sometimes produce omnidirectional sound fields. The directionality of a sound source can be modelled as a combination of the frequency of the emitted
sound and the geometry of the vibrating and non-vibrating parts of the sound source itself. The resulting pattern of sound radiation with angular location is called the *directivity* of the source [@beranek2012acoustics]. 

Directivity ($D$) describes the relative sound levels at a given angle $D(\theta)$
with relation to the level on-axis ($D(0)$, and thus $directivity = \frac{D_{\theta}}{D_{0}}$). Directivity
functions exist for a wide variety of sound sources that can be modelled analytically. A well-known example 
of a directivity function is of the piston in an infinite baffle. The piston is a circular surface of radius
$a$, vibrating back and forth about a hole in an infinite wall or baffle. The directivity is described by $\frac{2J_{1}(ka \times sin \theta)}{ka \times sin \theta}$, where $J_{1}$ is the Bessel's function of the first kind, $k$ is the wavenumber, where $k = \frac{2\pi f}{c}$, where $f$ is the frequency of the sound, and $c$ is the speed of sound. The angle of the receiver is $\theta$, which varies from 0-2$\pi$ radians in the azimuth. The ```beamshapes``` package aims to provide an easy interface to generate direcitivities of various sound sources.

Directivity functions can be used in a two-fold manner: 1) to deliberately engineer devices to suit particular specifications, e.g., loudspeaker sound fields [@beranek2012acoustics] and 2) to infer parameters of a sound source itself having assumed a relevant model, e.g., estimating the direction of call emission [@guarato2011method] and mouth aperture of bat echolocation calls [@jakobsen2013convergent;@kounitsky2015bats].


![Directivity patterns ($\frac{D_{\theta}}{D_{0}}$) of the currently implemented sound sources for a common set of $ka$ values for comparison, where $k$ is the wavenumber ($\frac{2\pi\:f}{c}$) and $a$ is the piston radius - or its equivalent. The directivity pattern shows the ratio between the off-axis sound level at ($\theta^{\circ}$) to the on-axis level at ($0^{\circ}$) in decibels. A) piston in an infinite baffle B) piston in a sphere with the half-angle aperture $\alpha = 30^{\circ}$, and where the piston radius $a=R\:sin \:\alpha$ C) oscillating cap of a sphere with  $\alpha=30^{\circ}$, and the equivalent piston radius is $a=R\:sin \:\alpha$ D) vibrating point on a sphere, here the $kR=ka$ for comparison with the other models. \label{fig:pistonsphereinf}](paper_related/piston_sphere_baffle.png)


# Statement of need

A host of published directivity functions exist in the literature, but it is my experience that their computational implementations remain as in-house scripts that are often in proprietary language platforms. To my knowledge, there are no computational implementations of multi-parameter sound source directivities. For instance, the ```levitate``` package [@levitate] implements directivities of 'simple' sound sources that can be defined with one parameter (e.g., piston in a baffle, circular ring). Since analytical solutions are available, the directivities of simple source models can be rapidly implemented and computed. In contrast to simple source models, the directivities of other models (as in this package) involve more parameters (e.g., piston in a sphere, oscillating cap of a sphere) and do not have analytical solutions. Their directivity calculations require numerical routines and run-time optimisation. The advantage of more involved source models is however their ability to capture aspects of experimental sound sources. In this paper, I present ```beamshapes```, a Python package that currently implements directivity patterns for two 'simple' and two 'involved' sound source models. As of this publication, ```beamshapes``` version 0.2.0 implements four sound sources: 1) piston in an infinite baffle, 2) point source on a sphere 3) oscillating cap of a sphere, and 4) piston in a sphere. 

Computational implementations of directivity functions often require long run-times due to the intensive numerical routines and arbitrary precision math involved. Long run-times hinder scientific projects in reducing the number of models and parameter space that can be explored. ```beamshapes``` boasts parallelised code to generate significant speed-ups in run-times. 

The availability of openly-available directivity functions will hopefully stir the acoustics, and specifically the bio-acoustics community to rigorously test and compare sound radiation with models. The availability of multiple implementations allows comparison of data with multiple models. Until recently, computational power has perhaps been a limiting factor in calculating directivity functions. Models with easily calculable outputs have thus been favoured (e.g., piston in an infinite baffle), especially in the field of bioacoustics [@strother1970acoustical;@mogensen1979sound]. 


Despite the recent availability of computational power, older, simpler models with limited biological relevance continue to dominate the field of bio-acoustics. For instance, the piston in an infinite baffle only predicts the beam-shape for a $\pm90^{\circ}$ range off-axis, and assumes front-back symmetry (\autoref{fig:pistonsphereinf} A). This is unrealistic for most vocalising animals, especially echolocating bats and odontocetes. However, the piston in an infinite baffle continues to be a standard reference model for multiple studies ranging over the past decade [@jakobsen2013convergent;@kounitsky2015bats;@macaulay2020high]. The piston in a sphere, for instance, recreates many of the highly directional central and side lobes and that is seen in bats, while also predicting sound radiation behind the source (\autoref{fig:pistonsphereinf} B). In contrast to the echolocation literature, the bird [@witkin1977importance;@larsen1990directionality;@brumm2002sound;@patricelli2007differences;@patricelli2008acoustic;@yorzinski2010birds]  and frog call literature [@gerhardt1975sound;@rodriguez2020sound] for instance has been dominated by quantitative characterisations of sound radiation with no attempts at directly comparing measurements to model predictions. Using model-based directivity patterns to infer source parameters allows the discovery of common parameter spaces that birds occupy and facilitates cross-species comparisons. The oscillating cap of a sphere and vibrating point on a sphere (\autoref{fig:pistonsphereinf} C,D) are two other models of potential relevance to bioacousticians attempting to describe the sound radiation of bird and frog calls for instance. 

Future releases of ```beamshapes``` are scheduled to include directivity patterns for additional models of interest such as rectangular cap of a sphere or piston in a closed finite baffle.


# Software packages used in this work
```beamshapes``` relies on the Python open-source ecosystem and is built on the ```numpy``` [@2020NumPy], ```scipy``` [@2020SciPy], ```sympy``` [@meurer2017sympy], ```mpmath``` [@mpmath] and ```flint``` [@hart2011flint] libraries. 

# Acknowledgements
TB thanks Gaurav Dhariwal for his continual math advice and inputs, and Tim Mellow for 
strong support through sharing his Mathematica code and timely clarifications. TB thanks Lasse Jakobsen and Holger Goerlitz for productive discussions leading to this project. This work was executed through a combination of TB's private time, a DAAD stipend, the IMPRS for Organismal Biology,  and then finally by a CASCB Medium grant.

# References
