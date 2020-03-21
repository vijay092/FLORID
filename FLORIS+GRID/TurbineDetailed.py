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
import ODETurbineDetailed
from tempfile import TemporaryFile
from scipy.integrate import ode


## Initialize FLORIS model
fi = wfct.floris_utilities.FlorisInterface("example_input.json")



###############################################################


# We want to evaluate the system on 300 linearly spaced times between 
t = np.linspace(0, 2, 200)                             
                              
nTurb = 4;
nState = 11;
x0 = 0.01*np.ones(nTurb*nState)

# Set the yaw angles to whatever you want 
yaw_angles = [0.0, 40.0, 50.0, 0]
fi.calculate_wake(yaw_angles=yaw_angles)

# New power
power_yaw = fi.get_turbine_power()
r = ode(ODETurbineDetailed.MultipleTurbs,).set_integrator('lsoda', method='bdf')
r.set_initial_value(x0, 0).set_f_params(power_yaw)
t1 = 20;
dt = 0.01; 
temp = [];
while r.successful() and r.t < t1:
    v = r.integrate(r.t+dt)
    line1, = plt.plot(r.t,v[0]*5,'r.',label='1st Turbine')
    line2, = plt.plot(r.t,v[11]*5,'g.',label='2nd Turbine')
    
plt.legend(handles=[line1, line2])
plt.ylabel('Power in MW')
plt.xlabel('Time in s')







