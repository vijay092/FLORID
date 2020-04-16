

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:43:33 2020

@author: sanja
"""



import numpy as np
from casadi import *
import matplotlib.pyplot as plt
import ODE_fullTurb
import collocation as cl
import random
from tempfile import TemporaryFile


# An augmented state vector that includes the running cost and
# the original state
nState = 24;
y = MX.sym('y',nState +1 ); 

# The original state
x = y[1:,0]

# The input
u = MX.sym('u',1)

# cost
ell = (x[0] - 0.7)**2 

f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9,\
f_10, f_11, f_12, f_13, f_14, f_15, f_16, f_17, f_18, f_19,\
            f_20, f_21, f_22, f_23 = ODE_fullTurb.TurbineEqns(x,u) ;



# Augmented Dynamics
# This includes the running cost, as well as a dummy index for time.
t_sym = MX.sym('t_sym',1)
f_aug = vertcat(ell,f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9,\
f_10, f_11, f_12, f_13, f_14, f_15, f_16, f_17, f_18, f_19,\
            f_20, f_21, f_22, f_23)
#f_aug = vertcat(ell,f_0)
f_aug_fun = Function('f_aug',[t_sym,y,u],[f_aug]) 

# constraint
g1 = x[17] - 0.1;
g2 = x[18] - 0.1
g = vertcat(g1,g2)
g_fun = Function('g',[y,u],[g])
    

T_max = 1e-3
C_list = [np.array([0,0.5,1])]


x0 = 0.08*np.ones(25,)
NumNodes = [20]
Markers = ['o','s','*','.']
for c in C_list:

    
    # c_str = ','.join([str(ci) for ci in c])
    # fig.suptitle('c = [%s]' % c_str)
    for N,m in zip(NumNodes,Markers):
        Time = np.linspace(0,T_max,N+1)
        X_opt,U_opt = cl.collocation_optimize(f_aug_fun,g_fun,Time,x0,1,c)



ref = np.ones(len(Time),)
fig,ax = plt.subplots(1,3,figsize=(18,5))
ax[0].plot(Time,X_opt[1],'-', linewidth=3.5)
ax[0].plot(Time,0.7*ref,'r*',linewidth=3.5)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Active Power')

ax[1].plot(Time,X_opt[18],'-',linewidth=3.5)
ax[1].plot(Time,0.1*ref,'g*',linewidth=3.5)
ax[1].set_xlabel('Time')
ax[1].set_ylabel(r'$P_{avg}$')

ax[2].plot(Time,X_opt[19],'-',linewidth=3.5)
ax[2].plot(Time,0.1*ref,'g*',linewidth=3.5)
ax[2].set_xlabel('Time')
ax[2].set_ylabel(r'$Q_{avg}$')



