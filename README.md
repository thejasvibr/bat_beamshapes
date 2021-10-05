[![Build Status](https://travis-ci.com/thejasvibr/bat_beamshapes.svg?branch=dev)](https://travis-ci.com/thejasvibr/bat_beamshapes)
[![status](https://joss.theoj.org/papers/820fd37f8255a8c533d6cc4c9475ecb5/status.svg)](https://joss.theoj.org/papers/820fd37f8255a8c533d6cc4c9475ecb5)

# Beamshapes

Calculate directivity of various sound source models and get their beamshapes.

This package is released under an MIT license. 

## Getting started with *beamshapes*

```
>>> import matplotlib.pyplot as plt 
>>> import numpy as np 
>>> import beamshapes
>>> from beamshapes import piston_in_infinite_baffle_directivity as PIB # short alias
>>> input_parameters = {'k':50, 'a':0.1}
>>> angles = np.linspace(-np.pi/2,np.pi/2,50)
>>> _, directionality = PIB(angles, input_parameters) # output the dB(D(theta)/D(on-axis))
>>> plt.figure();a0 = plt.subplot(111, projection='polar')
>>> plt.plot(angles, directionality) # plot the beamshape !!
```

For more detailed use-cases, check out the [example gallery online](https://beamshapes.readthedocs.io/en/latest/gallery_examples/index.html)!

## Installation 

*PyPi installation (>=version 0.2.1)*

```pip install beamshapes```

*Local installation instructions*
For the steps below to work you need to have a working Python installation that you can access from the command line. It is recommended to do the installation in a  [virtual environment](https://realpython.com/effective-python-environment/#virtual-environments). 

1. Clone the Github repository ```git clone https://github.com/thejasvibr/bat_beamshapes.git```
1. Change directories to the downloaded repo, and switch to the *dev* branch: ```git checkout dev``` 
1. Install the dependencies with ```pip install -r beamshapes/tests/requirements_test.txt```
1. Install *beamshapes* with ```python setup.py install```


## Detailed documentation 
For more details on the concepts and source documentation - please check out the [online docs](beamshapes.rtfd.io).


## Citation information 
If you use this package - please do consider citing the preprint: 

*Beleyur, T. (2021, August 5). beamshapes : a Python package to generate directivity patterns for various sound sources. https://doi.org/10.31219/osf.io/zc52t, version {VERSIONNUMBER-HERE} *

You can access the version number of the package being used with 
```
>>> import beamshapes
>>> print(beamshapes.__version__)
```

Bibtext format: 

```
@misc{beleyur_2021,
 title={beamshapes : a Python package to generate directivity patterns for various sound sources},
 url={osf.io/zc52t},
 DOI={10.31219/osf.io/zc52t},
 publisher={OSF Preprints},
 author={Beleyur, Thejasvi},
 year={2021},
 month={Aug}
}
```



## Future implementations
* Piston on a cylinder
* Rectangular piston on a prolate spheroid ([paper](https://asa.scitation.org/doi/pdf/10.1121/1.1778840?casa_token=wDAHTxJBISUAAAAA:MW-OSeGIkft-mces_mJgFBuyOhzI1qpPbc_7Xuu9EhDDD8CF8vnCIYaGyVivUb2qOpFda4GkPWto))


TODO:
* Add more examples 
* Remove 'a' parameter from piston in sphere directivity?  (next releases, this is backwards incompatible!)

