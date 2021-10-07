# Testing approach

All the source models in ```beamshapes``` are tested for their correctness by comparing
with published results. As of version 0.2, the models and their outputs are tested by
comparing the ```beamshapes``` output with plots in Beranek & Mellow 2012.

## Checking directivity patterns of source models 
Bernanek & Mellow 2012 calculate two types of directivity patterns: 1) the 'on-axis' pattern which describes the on-axis sound intensity across different *ka* values for instance, and 2) the 'full' directivity pattern, looking at the relative sound intensity (re. on-axis level) across the azimuth, and for different *ka* values. 

Tests pass if they are within the threshold of the published pattern. 