
from setuptools import setup, find_packages
import bat_beamshapes

 # link to test upload and fresh install on Test PyPi https://packaging.python.org/guides/using-testpypi/
 
version_number = bat_beamshapes.__version__

setup(name='bat_beamshapes',
     version=version_number,
     description='Bat beamshape modelling',
     long_description=open('README.md').read(),
     long_description_content_type="text/markdown",
     url='https://github.com/thejasvibr/bat_beamshapes.git',
     author='Thejasvi Beleyur',
     author_email='thejasvib@gmail.com',
     license='MIT',
     packages=find_packages(),
     install_requires=['numpy','sympy','mpmath','symengine','gmpy2',
        'scipy','matplotlib','tqdm'],
     zip_safe=False,
     include_package_data=True,
     classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Multimedia :: Sound/Audio :: Analysis',
        'Programming Language :: Python :: 3'
        ])
