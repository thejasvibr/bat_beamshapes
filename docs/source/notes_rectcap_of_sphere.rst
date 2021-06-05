Rectangular cap of a sphere
===========================

.. |Amn| replace:: :math:`A_{mn}`
.. |A0n| replace:: :math:`A_{0n}`
.. |Imn| replace:: :math:`I_{mn}`
.. |I0n| replace:: :math:`I_{0n}`
.. |dthetaphi| replace:: :math:`D(\theta,\pi)`
.. |d00| replace:: :math:`D(0,0)`
.. |date| date::

`Last updated`: |date|

The rectangular cap of a sphere models a square/rectangular part of a sphere's surface
vibrating (in contrast to the cap of  sphere, which models a circular cap). 

Here we'll try a top-down approach to understand the relevant variables and their corresponding code implementations in further detail. 

Directivity functions |dthetaphi| , |d00|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The off-axis directivity function |dthetaphi| is given by: 

.. math::
         
    D(\theta,\pi) = -\frac{4\pi}{k^{2}S}\sum^{N}_{n=0}\sum^{n/2}_{m=0}A_{mn}j^{n}P^{2m}_{n}(cos \:\theta)cos\:2m\phi \:\:\: (12.83)

The on-axis directivity function |d00| is given by: 

.. math::

    D(0,0) = -\frac{4\pi}{k^{2}S}\sum^{N}_{n=0}A_{0n}j^{n} \:\:\: (12.85)

The variables :math:`S, A_{mn}, A_{0n}` are defined below. 

.. math::

    S = 4R^{2}\Bigg(arctan\bigg(\frac{tan\:\alpha\:tan\:\beta}{sec^2\:\alpha + \sqrt{sec^{2}\:\alpha + tan^{2}\:\beta}}\bigg) \\ + arctan\bigg(\frac{tan\:\alpha\:tan\:\beta}{sec^2\:\beta + \sqrt{sec^{2}\:\beta + tan^{2}\:\alpha}}\bigg)\Bigg) \quad (12.69) \\
    \\ 

    A_{mn} = \frac{(2n+1)^2(n-2m)!I_{mn}}{j2\pi(n+2m)!\bigg(nh^{(2)}_{n-1}(kR) - (n+1)h^{(2)}_{n+1}(kR)\bigg)} \quad (12.75)

And |A0n| is all |Amn| where :math:`m=0`. 

Definition of |Imn|:
~~~~~~~~~~~~~~~~~~~~
|Amn| has the variable |Imn|, defined as: 

.. math::
    
    I_{mn} = \\
    \int^{arctan\frac{tan\:\beta}{tan\:\alpha}}_{0} cos\:2m\phi\:\int^{arctan\frac{tan\:\alpha}{cos\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\ 
    + \int^{\frac{\pi}{2}+arctan\frac{tan\:\alpha}{tan\:\beta}}_{\frac{\pi}{2}-arctan\frac{tan\:\alpha}{tan\:\beta}} cos\:2m\phi \:\int^{arctan\frac{tan\:\beta}{sin\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi  \\
    + \int^{\pi}_{\pi-arctan \frac{tan\:\beta}{tan\:\alpha}} cos\:2m\phi\:\int^{arctan\:\frac{tan\:\alpha}{-cos\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\
    \quad (12.76)

And |I0n| is all |Imn| where :math:`m=0`, given by: 

.. math::
    
    I_{0n} = \\
    \int^{arctan \frac{tan\:\beta}{tan\:\alpha}}_{0} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg(\frac{\cos\:\phi}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}\bigg)d\phi \\
    + \int^{\pi/2+arctan\frac{tan\:\alpha}{tan\:\beta}}_{\pi/2-arctan\frac{tan\\:alpha}{tan\:\beta}} \frac{tan\:\beta}{\sqrt{sin^{2}\:\phi + tan^{2}\:\beta}}
    P^{-1}_{n}\bigg(\frac{sin\:\phi}{\sqrt{sin^{2}\:\phi + tan^{2}\:\beta}}\bigg)d\phi \\
    + \int^{\pi}_{\pi-arctan\frac{tan\:\beta}{tan\:\alpha}} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg(\frac{-cos\:\phi}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}\bigg)d\phi \\
    \quad (12.77)


Computational implementation of |Imn|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The code implementation looks different from these equations, and that is because the solutions to some of the integrals seem to 
have been included already. 

Let's examine them one-by-one. 


|Imn| is coded with the equivalent of :code:`Int[m_,n_,beta_]`:

.. math::

    I_{mn} = \\
    \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{cos\:\phi_{1}}cos\:2m\phi_{1}\:P^{2m}_{n}(cos\theta_{1})\:sin\:\theta_{1} \\
    + \frac{2}{PQ}arctan\frac{sin\:\alpha}{sin\:\beta}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\beta}{sin\:\phi_{2}}cos\:2m\phi_{2}\:P^{2m}_{n}(cos\theta_{2})\:sin\:\theta_{2} \\
    + \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{-cos\:\phi_{3}}cos\:2m\phi_{3}\:P^{2m}_{n}(cos\theta_{3})\:sin\:\theta_{3} \\

where,

.. math::
    P=100, Q=100 \\
    \phi_{1}(p,\beta) = \frac{p-1/2}{P}arctan\:\frac{sin\:\beta}{sin\:\alpha} \\
    \phi_{2}(p,\beta) = \frac{\pi}{2} - arctan\:\frac{sin\:\alpha}{sin\:\beta} + 2\frac{p-1/2}{P}arctan\:\frac{sin\:\alpha}{sin\:\beta} \\
    \phi_{3}(p.\beta) = \pi - arctan\:\frac{sin\:\beta}{sin\:\alpha} + \frac{p-1/2}{P}arctan\:\frac{sin\:\beta}{sin\:\alpha} \\
    \theta_{1}(p,q,\beta) = \frac{q-1/2}{Q}\bigg(arctan\:\frac{tan\:\alpha}{cos\:\phi_{1}} \bigg) \\
    \theta_{2}(p,q,\beta) = \frac{q-1/2}{Q}\bigg(arctan\:\frac{tan\:\beta}{sin\:\phi_{2}} \bigg) \\
    \theta_{3}(p,q,\beta) = \frac{q-1/2}{Q}\bigg(arctan\:\frac{tan\:\alpha}{-cos\:\phi_{3}} \bigg) \\

Comparing the math |Imn| definitions and code implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's begin to compare the code implementation of the solution with the original |Imn| terms:

Math |Imn| term1:

.. math::

    \int^{arctan\frac{tan\:\beta}{tan\:\alpha}}_{0} cos\:2m\phi\:\int^{arctan\frac{tan\:\alpha}{cos\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\ 
    
    implemented\:as\: :

    \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{cos\:\phi_{1}}cos\:2m\phi_{1}\:P^{2m}_{n}(cos\theta_{1})\:sin\:\theta_{1} \\

|Imn| term2:

.. math::

    \int^{\frac{\pi}{2}+arctan\frac{tan\:\alpha}{tan\:\beta}}_{\frac{\pi}{2}-arctan\frac{tan\:\alpha}{tan\:\beta}} cos\:2m\phi \:\int^{arctan\frac{tan\:\beta}{sin\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi  \\    
    
    implemented\:as\: :

    \frac{2}{PQ}arctan\frac{sin\:\alpha}{sin\:\beta}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\beta}{sin\:\phi_{2}}cos\:2m\phi_{2}\:P^{2m}_{n}(cos\theta_{2})\:sin\:\theta_{2} \\

|Imn| term 3:

.. math::

    \int^{\pi}_{\pi-arctan \frac{tan\:\beta}{tan\:\alpha}} cos\:2m\phi\:\int^{arctan\:\frac{tan\:\alpha}{-cos\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\
    
    implemented\:as\: :

     \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{-cos\:\phi_{3}}cos\:2m\phi_{3}\:P^{2m}_{n}(cos\theta_{3})\:sin\:\theta_{3} \\


Definition of |I0n|:
~~~~~~~~~~~~~~~~~~~~

|I0n| (all `m=0`) terms:

.. math:: 

    \int^{arctan\:\frac{tan\:\beta}{tan\:\alpha}}_{0}\:\frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg( \frac{cos\:\phi}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}\bigg)d\phi \\
    + \int^{\pi/2 + arctan\:\frac{tan\:\alpha}{tan\:\beta}}_{\pi/2-arctan\:\frac{tan\:\alpha}{tan\:\beta}}\:\frac{tan\:\beta}{\sqrt{sin^{2}\:\phi+tan^{2}\:\beta}}
    P^{-1}_{n}\bigg( \frac{sin\:\phi}{\sqrt{sin^{2}\:\phi+tan^{2}\:\beta}}\bigg)d\phi \\
    + \int^{\pi}_{\pi-arctan\:\frac{tan\:\beta}{tan\:\alpha}}\:\frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg( \frac{-cos\:\phi}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}\bigg)d\phi \\

Computational implementation of |I0n|:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    \frac{1}{P} arctan\:\frac{sin\:\beta}{sin\:\alpha} \sum^{P}_{p=1} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi_{1}+tan^{2}\alpha}}
    P^{-1}_{n}\bigg(\frac{cos\:\phi_{1}}{\sqrt{cos^{2}\:\phi_{1}+tan^{2}\:\alpha}}\bigg) \\
    + \frac{2}{P} arctan\:\frac{sin\:\alpha}{sin\:\beta}\sum^{P}_{p=1}\frac{tan\:\beta}{\sqrt{sin^{2}\:\phi_{2}+tan^{2}\beta}}
    P^{-1}_{n}\bigg(\frac{sin\:\phi_{2}}{\sqrt{sin^{2}\:\phi_{2}+tan^{2}\:\beta}}\bigg) \\
    \frac{1}{P} arctan\:\frac{sin\:\beta}{sin\:\alpha} \sum^{P}_{p=1} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi_{3}+tan^{2}\alpha}}
    P^{-1}_{n}\bigg(\frac{-cos\:\phi_{3}}{\sqrt{cos^{2}\:\phi_{3}+tan^{2}\:\alpha}}\bigg) \\

Comparing |I0n| in code and definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The |I0n| term doesn't have any substitutions/solutions in the code implementation it. The main difference is that the numerical integration is built into the 
implementation (an alternative approach to use a inbuilt numerical integration routine). 


Acknowledgements
----------------
Thanks to Tim Mellow for sharing the `Mathematica` code behind the textbook figures. 
