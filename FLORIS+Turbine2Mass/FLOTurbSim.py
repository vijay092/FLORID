# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 13:12:19 2019

@author: sanjana
"""


import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import scipy.integrate as spi
import numpy as np
from scipy.io import loadmat
import ODETurbine2Mass
from tempfile import TemporaryFile
from scipy.integrate import ode


## Initialize FLORIS model
fi = wfct.floris_utilities.FlorisInterface("example_input.json")



###############################################################


# We want to evaluate the system on 300 linearly spaced times between 
# t=0 and t=300.
t = np.linspace(0, 2, 200)
                             
                              
nTurb = 4;
nState = 3;
x0 = np.ones(nTurb*nState)



# Set the yaw angles to whatever you want 
yaw_angles = [25.0, 0, 25.0, 0]
fi.calculate_wake(yaw_angles=yaw_angles)

# New power
power_yaw = fi.get_turbine_power()

r = ode(ODETurbine2Mass.MultipleTurbs,).set_integrator('lsoda', method='bdf')
r.set_initial_value(x0, 0).set_f_params(power_yaw)
t1 = 50;
dt = 0.01; 
temp = []
while r.successful() and r.t < t1:
    v = r.integrate(r.t+dt)
    temp.append(v)
    print("t =", r.t)
xt = np.asarray(temp)
plt.plot(xt)





