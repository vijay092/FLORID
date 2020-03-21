# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 10:15:10 2020

@author: sanja
"""


import numpy as np
from casadi import *
import matplotlib.pyplot as plt
import ODETurbineDetailed
import collocation as cl
import random
# An augmented state vector that includes the running cost and
# the original state

y = MX.sym('y',12); 

# The original state
x = y[1:,0]

# The input
u = MX.sym('u',1)

# cost
ell = x[0]**2

f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10 = ODETurbineDetailed.TurbineEqns(x,u) ;


# Augmented Dynamics
# This includes the running cost, as well as a dummy index for time.
t_sym = MX.sym('t_sym',1)
f_aug = vertcat(ell,f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10)
f_aug_fun = Function('f_aug',[t_sym,y,u],[f_aug]) 

# constraint
g = []
g_fun = Function('g',[y,u],[g])


T_max = 4.5
C_list =  [
          np.array(collocation_points(2,'legendre')),
          np.array([0,.5,1])]


x0 = np.ones(12)

NumNodes = [10,20,30]
Markers = ['o','s','*','.']
for c in C_list:

    fig,ax = plt.subplots(1,2,figsize=(12,4))
    
    c_str = ','.join([str(ci) for ci in c])
    fig.suptitle('c = [%s]' % c_str)
    for N,m in zip(NumNodes,Markers):
        Time = np.linspace(0,T_max,N+1)
        X_opt,U_opt = cl.collocation_optimize(f_aug_fun,g_fun,Time,x0,1,c)

        for i in range(2):
            ax[i].plot(Time,X_opt[i+1],marker=m)
            ax[i].set_xlabel('Time')
            ax[i].set_ylabel(r'$x_{%d}$' % (i+1,))




