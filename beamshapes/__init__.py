'''Loads a bunch of directivities into the packages main namespace.
All others will need to be deliberately loaded. 
'''

from beamshapes.version import __version__

# from beamshapes.beamshape_predictions import * 
# from beamshapes.sim_reclevels  import * 
from beamshapes.piston_in_sphere import piston_in_sphere_directivity
from beamshapes.cap_in_sphere import cap_in_sphere_directivity
from beamshapes.piston_in_infinite_baffle import piston_in_infinite_baffle_directivity
from beamshapes.point_source_on_a_sphere import point_source_on_a_sphere_directivity
