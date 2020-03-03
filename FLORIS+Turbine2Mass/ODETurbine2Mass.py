# -*- coding: utf-8 -*-

import numpy as np

def Cp(lamda, beta):
    
    c1 = 0.22; c2=116; c3 = 0.4;
    c4 = c5 = 0; c6= 5; c7=12.5; c8=0.08; c9=0.035; c10=0;
    
    cp = c1*( (c2/(lamda + c8*beta)) - c2*c9/(beta**3 + 1) - c3 * beta - 
             c4 * beta**c5  - c6 )\
    * np.exp(-c7/(lamda + c8*beta)) + c10*lamda;
    
    return cp 

def MultipleTurbs(t,x,power_aero):
    nTurb = 4; nState = 3; 
    F= np.zeros(nTurb*nState)
    for i in range (nTurb):
        F[i:i+nState] = TurbODE(t,x[i:i+nState],power_aero[i]);
        
    return F 

def TurbODE(t,x,power):
    
    '''
    This function is a mapping between the control handle "Kopt" and the
    the power output of the turbine. 
    
    Inputs:
        t: time in s
        x: state in p.u.
        Kopt: Control handle
    
    Output:
        Power: Generator power output in p.u.
    
    '''
    
    F = np.zeros(len(x))
    
    # Rename variables
    wg = x[0];
    theta_tw = x[1];
    wt = x[2];
    
    # Parameters
    Ht = 4; 
    Hg = 0.1*Ht; 
    ksh = 0.3; 
    csh = 0.01;
    wel = 2*np.pi*60;
    Kopt = .5787
    
    rho = 1.225 
    R_turb = 58.6;
    V = 12;
    Prated = 5e6;
    GB = 145.5;
    #wtB = wel/(2*GB);
    lamda = wt*R_turb/V;
    beta = 0;
    
    Tmbase =  GB * Prated * 1/(wel/2);
    
    # Intermediate variables
    #Tm = 0.5 * rho * np.pi* R_turb**2 * Cp(lamda,beta) * V**3 /(wt*Tmbase);
    Tm = power/(5e6*wt);
    # Equations
    F[0] = 1/(2*Hg) * (ksh*theta_tw + \
                           csh*wel* (wt - wg) - Kopt*x[2]**2 );
    
    F[1] = wel*(wt - wg);
    
    F[2] = 1/(2*Ht) * (Tm - ksh * theta_tw \
                             - csh*wel*(wt - wg));
    
  
    return F