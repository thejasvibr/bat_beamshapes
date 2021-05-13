## Speeding up symbolic calculations

Right now I'm stuck with SymPy calcultions that are taking for ever. 
I have 3 options right now, 

## What are other people using for the same task?


1. Mellow, Tim: *Mathematica*

1. Aarts Ronald, Janssen, Augustus J. E. M.: ronald.m.aarts@philips.com : *NO reply*
    * Comparing Sound Radiation from a Loudspeaker with that from a Flexible Spherical Cap on a Rigid Sphere 2011

1. Wojciech P.Rdzanek: wprdzank@ur.edu.pl :*replied, uses Mathematica*
    * The acoustic power of a vibrating clamped circularplate revisited in the wide low frequency rangeusing expansion into the radial polynomials 2016
1. Mark Poletti, NZ : *replied, uses MATLAB*


## Trying out FriCAS

* I've been getting a hang of FriCAS - it seems to be a very powerful, yet low level language. 
* Write and read text from/to files using the writeFile syntax http://fricas-wiki.math.uni.wroc.pl/SandBoxTextFiles

## I've been messing up the plain Hankel function of the 2nd kind nd the *spherical Hankel function of the 2nd kind!!!*
## Not having a formal math education means, when I encounter the derivative of a legendre function (Pn'(x)) -- I have no idea whether the derivative is with respect to n or x :P

* FriCAS is pretty cool as a system, though the learning curve made me stick to Python in the end. The confusion between multiple derivatives. 


## Troubleshooting 

*copy from ```ts_pistonsphere.py```*
### Troubleshooting piston in a sphere

> For ka>3 the calculations show values that deviate from the 
textbook groundtruth by 2-5 dB -- which is not ignorable!!
> I suspected the problem may come from:
    * low mpmath.mp.dps (decimal places) -- changing from 50-400 had no
    effects, which is odd
    * low N (matrix size/number of terms calculcated) -- changing from 
    baseline of 12+f(ka) --> 15+f(ka) had no effect
    * the directionality calculations were done with numpy pre 5th may, 
    and then changed to mpmath backend --> no effect. 

> 2021-05-65: I now suspect the problem lies perhaps with the quadrature 
terms. What if the quadrature is not 'accurate' enough? Here I'll test this idea
    * Some points to support this idea. The default quadrature method behind
    mpmath.quad is the 'tanh-sinh' algorithm. Instead of directly lambdifying 
    the Imn term into a standard mpmath.quad function I 'manually' made a 
    quadrature function for it to manipulate the options. 
    * For dps 200. Using the default 'tanh-sinh' leads to an estimated error 
    of e-203, while using 'gauss-legendre' leads to an estimated error of e-382. 
    Perhaps this is where the error is arising from. I noticed the integration 
    error increases with increasing m,n values. Perhaps this is why for bigger ka's, 
    (ka>3), the predictions get messier than for small ka's? There is at least 
    a connection here. 

The exact quadrature algorithm used seems to play a big difference. 

*The Imn term quadrature algorithm is NOT the problem* -- though choosing gauss-legendre
reduced execution time by 1/2!!!!

2021-05-06: Despite troubleshooting for so long on the piston in a sphere, I"m not having any luck understanding where the issue is coming from. 

At this point I have two options:
* stick with the piston in a sphere, but run all calculations with ka <= 3 
OR
* continue to troubleshoot with SAGE. 
* The challenge here is basically that 


2021-05-12:
I'm realising there is some effect of *N*, the number of terms used to get the beamshape. 
While the overall effect of increasing N for ka=3 is weak (an decrease of ~0.2 dB max|error|,
and decrease of 0.7 dB of mean|error|), the LUsolve residuals dramatically increase 2-3 orders. 

This tells me there may be issues with the numerical integration in general. Perhaps the numerical integration gets more messy with N terms.


Also -- WHAT ABOUT (PY)-FLINT. Check out  this [page](https://fredrikj.net/blog/2018/11/announcing-python-flint-0-2/), and this [one](https://fredrikj.net/python-flint/index.html)?


Issues faced while installing Flint
> first need to install MPIR. Downloaded source from Github.
    * the ```./configure``` command didn't work because the previous steps weren't complete -- needed to run ```autoreconf -i``` and then the ./configure 
> Also don't forget you have to install GMP+ MPFR too. 

One thing I realised is the installation flow for 'source' installations is 
> ./configure
> make 
> make check 
> make install --> this step typically needs permissions, and so `sudo make install` with 
the password prompt is what is actually needed. 

Installing python-flint with ```pip``` failed. Installing python-flint from Git source failed. Installing with Conda -- This worked!!! 

### What it's like to use py-flint

> Incredibly fast !!! eg. show the example of the Lm integral (eqn. 12.018) for a dps=500, it happens nearly instantaneously in py-flint, while with mpmath it takes 10's of seconds. 
2021-05-12
16:21 - I'm getting different values for P'm(cos(alpha)) between my Flint implementation using eqn. 12.98 and the native legendre_p(n,z).diff(z) version from sympy. 

16:44 -- and that's because there's a typo in eqn. 12.98 --> it should be sin(theta)^2 instead of sin(theta)

22:52 : FLINT is really, really fast now that I think of it. For dps=150, and ka=10 the Flint version of the code runs in ~9.5mins while the parallelised mpmath version takes about that long!!




