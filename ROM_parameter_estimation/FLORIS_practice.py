# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 19:54:52 2020

@author: sanja
"""


# learn FLORIS

import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

# Initialize FLORIS model
fi = wfct.floris_utilities.FlorisInterface("example_input.json")
fi.calculate_wake()
print(fi.floris.farm.wind_direction)

print(fi.floris.farm.turbines[3].yaw_angle)


# Ok, whatever, now adjust floris
fi.floris.farm.turbines[3].yaw_angle