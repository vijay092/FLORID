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
import Grid
from scipy.io import loadmat
import ODEs
from tempfile import TemporaryFile
from scipy.integrate import ode


## Initialize FLORIS model
fi = wfct.floris_utilities.FlorisInterface("example_input.json")



###############################################################



# This was originally written in MATLAB
x = loadmat('2area.mat')
plt.close()
fig, ax = plt.subplots(1, 1, figsize=(8, 4))

# We want to evaluate the system on 300 linearly spaced times between 
# t=0 and t=300.
t = np.linspace(0, 2, 200)

No_mach = 4;
No_turb = 4;
No_states_mach = 4;
No_states_turb = 11;
N_full = No_mach*No_states_mach + No_states_turb*No_turb                                          
                              
# # # We simulate the system 
Ybus = x['Ybus']

x0 = np.ones(N_full)
x0[3] = 377;
x0[7] = 377;
x0[11] = 377;
x0[15] = 377;
# v = spi.odeint(Grid.NonLinSystem,x0, t, args=(fi,Ybus))
# ax.plot(t,v)

# Set the yaw angles to whatever you want 
yaw_angles = [25.0, 0, 25.0, 0]
fi.calculate_wake(yaw_angles=yaw_angles)

# New power
power_yaw = fi.get_turbine_power()

r = ode(Grid.NonLinSystem,).set_integrator('lsoda', method='bdf')
r.set_initial_value(x0, 0).set_f_params(power_yaw, Ybus)
t1 = 50
dt = 0.01; 
temp = []
while r.successful() and r.t < t1:
    v = r.integrate(r.t+dt)
    temp.append(v)
    print("t =", r.t)
xt = np.asarray(temp)
plt.plot(xt[:,3])





