[![Build Status](https://travis-ci.org/thejasvibr/bat_beamshapes.svg?branch=dev)](https://travis-ci.org/thejasvibr/bat_beamshapes)

# Beamshapes

Calculate directivity of various sound source models and get the beamshapes!
This package is released under an MIT license. 

## Getting started with *beamshapes*

```
>>> import matplotlib.pyplot as plt 
>>> import numpy as np 
>>> import bat_beamshapes
>>> from bat_beamshapes import piston_in_infinite_baffle_directionality as PIB # short alias
>>> input_parameters = {'k':50, 'a':0.1}
>>> angles = np.linspace(-np.pi/2,np.pi/2,50)
>>> _, directionality = PIB(angles, input_parameters) # output the dB(D(theta)/D(on-axis))
>>> plt.figure();a0 = plt.subplot(111, projection='polar')
>>> plt.plot(angles, directionality) # plot the beamshape !!
```

For more detailed use-cases, check out the [example gallery online](https://beamshapes.readthedocs.io/en/latest/gallery_examples/index.html)!

## Installation 

*Pre-PyPi installation instructions*

1. Clone the Github repository ```git clone https://github.com/thejasvibr/bat_beamshapes.git```
1. Change directories to the downloaded repo, and wwitch to the *dev* branch: ```git checkout dev``` 
1. Install the package with ```pip install ./```



## Detailed documentation 
For more details on the concepts and source documentation - please check out the [online docs](beamshapes.rtfd.io).


## Citation information 
If you use this package - please do consider citing it! 

*Beleyur, Thejasvi, , 2021, beamshapes: directivity patterns for multiple sound sources {version number here}, https://github.com/thejasvibr/bat_beamshapes*

You can access the version number of the package being used with 
```
>>> import bat_beamshapes
>>> print(bat_beamshapes.__version__)
```

## Future implementations
* Piston on a cylinder
* Rectangular piston on a prolate spheroid ([paper](https://asa.scitation.org/doi/pdf/10.1121/1.1778840?casa_token=wDAHTxJBISUAAAAA:MW-OSeGIkft-mces_mJgFBuyOhzI1qpPbc_7Xuu9EhDDD8CF8vnCIYaGyVivUb2qOpFda4GkPWto))


TODO:
* Verification tests to check output with Beranek & Mellow 2012 graphs:
    * piston in a sphere pending.s
* change official name from 'bat_beamshapes' to just 'beamshapes'
* Ideas for package logo:
    * butterfly - N suggests 

DONE:
* ~~Uniform output format to a np.array that has 20*log10 D_theta/D_0~~ Done.
* ~~Include option to run FLINT version of piston in sphere! This will help at least Linux users a lot.~~
* ~~Cap of sphere - implement parallelisation when ka>5~~
* ~~Include module-level docs too.~~
* ~~Check if you really need symengine and gmpy2 for the code to run -- is sruntime really affected? -- Not needed. Speed not affected much in tests.~~
