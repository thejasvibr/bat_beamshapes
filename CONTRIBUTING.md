#  Contributing 

## Be nice, patient and respectful
Like it says above, while interacting with each other, contributing or raising issues - be nice, patient and respectful. We are all in it for the joy of cool science and software - let's make it about that. 

## Issues 
If you are facing issues running the software, or suspect bugs please report the following in the issue raised:

* ```beamshapes``` version 
* the smallest possible reproducible example that consistently generates the bug in your system
* OS and version
* Your Python type (direct installation, anaconda, etc.)

## Directivity implementations

If you're in doubt about whether an implementation is of relevance to the ```beamshapes``` package, raise an issue with the proposed model and let's discuss. For now, ```beamshapes``` is specifically looking to add source models that have a 'rich' front-back directivity ie. that are not front-back symmetric. 

If you wish to contribute a directivity implementation for a source-model see the generalised workflow page [here](https://beamshapes.readthedocs.io/en/latest/general_workflow.html). 
For reasons of code homogeneity, I request that you keep the template structure suggested in the [how-to](https://beamshapes.readthedocs.io/en/latest/general_workflow.html) page. All implementations must have tests that at least check the coded implementation matches the published results of the original paper/book that described the model. (If you are planning to implement a directivity that isn't published or have calculated yourself - we still need tests, but I'd be curious to hear how to check the validity of the implementation's results!)

##  Unit-tests
All unit-tests in the ```beamshapes``` package are written using the ```unittest``` package. Be aware that running all the tests will take a decent amount of time on an average PC/laptop (anywhere upwards of 20 minutes!).

### Running unit-tests
To run the unittests, first clone and install the ```beamshapes``` package (see the [README](README.md)).To run the tests you will need to install all the packages in the ```tests/requirements_test.txt``` file. From the root of the directory, open your Python command line tool of choice and run:

```python -m unittest```  

To run specific unit tests, use the commande ```python -m unittest beamshapes/tests/{name-of-test-here.py}```, where *{name-of-test-here.py}* is the 
test you want to run.

### Writing unit-tests
The unit-tests in ```beamshapes``` should test the computational implementation for correctness based on previously published results. See [here](beamshapes/tests/testing.md) for more on the testing approach. 

### Test coverage
If you'd like to estimate test coverage:

* Install the [```coverage```](https://coverage.readthedocs.io/en/6.0.2/) package: ```pip install coverage``` 
* Change directories and enter the repo ```cd bat_beamshapes```
* Run the tests with ```coverage``` and ```unittest``` with: ```coverage run -m unittest discover```
* To generate an overall report on coverage: ```coverage report```

## Documentation

### Building docs locally
To build the documentation locally, you must first install all the required packages:

```pip install -r docs/docs_requirements.txt``` 

And then:

* Move into the ```docs``` folder with ```cd docs```
* Make an empty ```build``` folder with ```mkdir build```
* Generate the docs with ```sphinx-build source build```

All the formatted web-pages are in the ```build/``` directory. Start with the ```index.html``` to begin browsing from the landing page. 



