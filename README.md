# Beamshapes

Calculate directivity/beamshapes of various sound production models.


TODO:
* Verification tests to check output with Beranek & Mellow 2012 graphs:
    * piston in a sphere pending.s
* change official name from 'bat_beamshapes' to just 'beamshapes'


DONE:
* ~~Uniform output format to a np.array that has 20*log10 D_theta/D_0~~ Done.
* ~~Include option to run FLINT version of piston in sphere! This will help at least Linux users a lot.~~
* ~~Cap of sphere - implement parallelisation when ka>5~~
* ~~Include module-level docs too.~~
