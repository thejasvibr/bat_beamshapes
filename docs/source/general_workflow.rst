Implementing a directivity function
===================================
This document describes the general workflow behind implementing the directivity 
function for a published sound radiation model. 

Pre-coding
~~~~~~~~~~
#. Find a suitable model with verifiable results (in the form of plots/numerical results)
#. Try and request underlying code for later comparison in case of deviations in results

Coding
~~~~~~
#. Implement as much of the model's components as SymPy objects. 
#. Always add the equation numbers in a comment above the variable/equation. 
#. Name the objects and intermediate variables to be as similar as possible to the names used in the publication. There will be some terms that are extremely long, split them where and when possible. Splitting long terms helps in later verification, and makes for pretty code. 
#. Convert the objects into functions using :code:`lambdify` and the appropriate backend (sympy, scipy, mpmath)
#. Put all the relevant code for the source model into one module. If necessary implement additional convenience functions in a separate model. 
#. Try switching between backends if the code throws unexpected errors or is taking too long. 
#. The final directivity function should follow the `{name-of-the-model}_directivity`, and accept 2 inputs. The first input :code:`angles` should be an array/list-like object describing the location/s at which :math:`\frac{D_{\theta}}{D_{0}}` are to be calculated, and a :code:`params` dictionary object - which parametrises the source model. 

Result verification 
~~~~~~~~~~~~~~~~~~~
#. Replicate the key plots from the publication to ensure the correctness of your implementation. 
#. In case you can't directly compare the results by running the original publication's code - use a data digitisation tool like `WebPlotDigitizer <https://apps.automeris.io/wpd/>`_ to extract data from plots directly. The digitisation itself can have its own errors - so make sure to account for it while looking at discrepancies between your output and the published data.
#. If the results don't match - check the source of discrepancy. There is a high chance of a typo over the course of entering 10's of equations and variables! 
# In case the code itself reflects the underlying equations correctly: check for 'patterns' in matching or discrepancy to try and isolate the problem. Does the match improve with the angle, does it get better at lower/higher frequencies, is it worse when the emitter is larger - or with increasing decimal precision (when using a :code:`mpmath` backend implementation). 


Optimisation
~~~~~~~~~~~~
* Certain models involve numerical methods that can take a long time to run (integration, summations, etc.). Check how to reduce run times without affecting output correctness eg. parallelising the code
* SymPy based implementations and compatible backends are strongly preferred. However, if a correct implementation takes too long (5-10 minutes each time) - consider switching to an alternative library like :code:`flint`.

To Do
~~~~~
* What when results don't match -- possible typo in the original publication/code.



