# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:43:33 2020

@author: sanja
"""



import numpy as np
from casadi import *
import matplotlib.pyplot as plt
import ODETurbineDetailed
import ODE2Mass
import collocation as cl
import random
from tempfile import TemporaryFile

_ = outfile.seek(0)

h = np.load(outfile)
# An augmented state vector that includes the running cost and
# the original state
nState = 9;
y = MX.sym('y',nState +1 ); 

# The original state
x = y[1:,0]

# The input
u = MX.sym('u',1)

# cost
ell = (x[0] - 0.8)**2

f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8 = ODE2Mass.TurbineEqns(x,u) ;
#f_0 = ODE2Mass.TurbineEqns(x,u) ;


# Augmented Dynamics
# This includes the running cost, as well as a dummy index for time.
t_sym = MX.sym('t_sym',1)
f_aug = vertcat(ell,f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8)
#f_aug = vertcat(ell,f_0)
f_aug_fun = Function('f_aug',[t_sym,y,u],[f_aug]) 

# constraint
g = []
g_fun = Function('g',[y,u],[g])
    

T_max = 0.01
C_list = [np.array([0,0.5,1])]


x0 = h
NumNodes = [10]
Markers = ['o','s','*','.']
for c in C_list:

    fig,ax = plt.subplots(1,3,figsize=(12,4))
    
    c_str = ','.join([str(ci) for ci in c])
    fig.suptitle('c = [%s]' % c_str)
    for N,m in zip(NumNodes,Markers):
        Time = np.linspace(0,T_max,N+1)
        X_opt,U_opt = cl.collocation_optimize(f_aug_fun,g_fun,Time,x0,1,c)

        for i in [0,1,2]:
            
            ax[i].plot(Time,X_opt[i+1],marker=m)
            ax[i].set_xlabel('Time')
            ax[i].set_ylabel(r'$x_{%d}$' % (i+1,))



