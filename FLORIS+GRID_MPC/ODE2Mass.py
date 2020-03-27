# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 23:07:05 2020

@author: sanja
"""


# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:58:23 2019

@author: sanja
"""

import numpy as np
import scipy.linalg as la
from casadi import *
import floris.tools as wfct
import floris.tools.visualization as vis
import floris.tools.cut_plane as cp
fi = wfct.floris_utilities.FlorisInterface("example_input.json")



def Cp(lamda, beta):
    
    c1 = 0.22; c2=116; c3 = 0.4;
    c4 = c5 = 0; c6= 5; c7=12.5; c8=0.08; c9=0.035; c10=0;
    
    cp = c1*( (c2/(lamda + c8*beta)) - c2*c9/(beta**3 + 1) - c3 * beta - 
             c4 * beta**c5  - c6 )\
    * np.exp(-c7/(lamda + c8*beta)) + c10*lamda;
    
    return cp 

def MultipleTurbs(t,x,power_aero):
    nT = [0,1,2,3];
    nState = 11; 
    F= np.zeros(len(nT)*nState)
    for i in nT:
        F[nState*i:nState*i+nState] = TurbineEqns(t,x[nState*i:nState*i+nState],power_aero[i]);
        
    return F 



def TurbineEqns(x,u):
    
    '''
    This function takes in the state and outputs 
    '''

    # Miscellaneous
    Lm = 4;
    Lss = 1.01*Lm;
    Lrr = 1.008*Lss;
    Kmrr = Lm/Lrr;
    Ls = Lss - Lm*Kmrr;
    Rs = 0.005;
    Rr = 1.1*Rs;
    R2 = Kmrr**2 * Rr;
    R1 = Rs + R2;
    Tr = 10;
    wel = 2*np.pi*60;
    Kopt = 0.5787;
    ws = 1;
    Qsref = 0;
    vqs = 1;
    vds = 0;
    N = 1;
    
    # Controller
    Kiq = -1;
    Kte = -1.5;
    Kid = -0.5;
    Kqs = 1;
    Tiq = 0.0025;
    Tte = 0.025;
    Tid = 0.005;
    Tqs = 0.05;
    
    # Mechanical
    Ht = 4; 
    Hg = 0.1*Ht; 
    ksh = 0.3; 
    csh = 0.01;

    # Rename the states of the turbine for ease of analysis


    iqs = x[0];
    ids = x[1];
    eqs = x[2];
    eds = x[3];
    phiTe = x[4];
    phiIq = x[5];
    phiQs = x[6];
    phiId = x[7];
    wg = x[8];
    
    
    R1 = R1/N;
    Ls = Ls/N;
    Kid  = Kid/N;
    Kiq = Kiq /N;
    R2 = R2/N;
    ksh = ksh*N;
    csh = csh*N;
    Ht = Ht*N;
    Hg = Hg*N;
    

    rho = 1.225 ;
    R_turb = 128/2;
    V = 12;
    Prated = 5e6;
    GB = 145.5;
    #wtB = wel/(2*GB);
    wt = wg/GB
    lamda = (wt)*R_turb/V;
    beta = 0.01;
    Tmbase =  GB * Prated * 1/(wel/2);
    Teref = u;
    # Intermediate variables
    Pm = 0.5 * rho * np.pi* R_turb**2 * \
            Cp(lamda,beta) * V**3 /(Tmbase);
 
    # Intermediate variables
    Te =    (eqs /ws) * iqs + (eds /ws) * ids;

    Qs =  -vqs * ids + vds * iqs;

    iqr  = - eds - Kmrr * iqs;
    
    idr =  eqs - Kmrr * ids;
    
    vqr =  Kiq* Kte * (Teref - Te) + Kiq * Kte /Tte * phiTe \
            - Kiq * iqr + Kiq/Tiq * phiIq; 
    
    vdr = +  Kid* Kqs * (Qsref - Qs) + Kid * Kqs /Tqs * phiQs \
            - Kid * idr + Kid/Tid * phiId; 
    
    
    # Differential equations:
    f_0 = wel / Ls * ( -R1 * iqs + ws * Ls * ids  \
                              + wg/ws * eqs - 1/(Tr * ws) * eds  \
                              - vqs + Kmrr * vqr);
      
    f_1 = wel / Ls * ( -R1 * ids - ws * Ls * iqs +\
                              wg/ws * eds - 1/(Tr * ws) * eqs \
                              - vds + Kmrr * vdr);
    
    f_2 = wel * ws * ( R2 * ids  + (1 - wg/ws) * eds - \
                              1/(Tr * ws) * eqs  - Kmrr * vdr);
    
    f_3 = wel * ws * (- R2 * iqs  - (1 - wg/ws) * eqs - \
                              1/(Tr * ws) * eds  + Kmrr * vqr);
    
    f_4 = Teref - Te;
    
    f_5 = Kte * (Teref - Te)  + Kte /Tte * phiTe  - iqr;

    f_6 = Qsref - Qs;
    
    f_7 = Kqs * (Qsref - Qs)  + Kqs /Tqs * phiQs  - idr;
    
    f_8 = 1 /(2*Hg ) * (Pm  - wg*Te);
    
    return f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8

 

    