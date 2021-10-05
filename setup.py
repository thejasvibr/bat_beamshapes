
from setuptools import setup, find_packages

from beamshapes.version import __version__

# link to test upload and fresh install on Test PyPi https://packaging.python.org/guides/using-testpypi/
 
version_number = '0.2.0'

setup(name='beamshapes',
     version=version_number,
     description='Acoustic beamshape modelling for various sources',
     long_description=open('README.md').read(),
     long_description_content_type="text/markdown",
     url='https://github.com/thejasvibr/bat_beamshapes.git',
     author='Thejasvi Beleyur',
     author_email='thejasvib@gmail.com',
     license='MIT',
     packages=find_packages(),
     install_requires=['joblib','numpy','sympy','mpmath',
        'scipy','matplotlib','tqdm'],
     zip_safe=False,
     include_package_data=True,
     classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Multimedia :: Sound/Audio :: Analysis',
        'Programming Language :: Python :: 3'
        ])