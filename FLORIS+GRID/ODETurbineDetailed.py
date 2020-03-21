# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:58:23 2019

@author: sanja
"""

import numpy as np
import scipy.linalg as la

def MultipleTurbs(t,x,power_aero):
    nT = [0,1,2,3];
    nState = 11; 
    F= np.zeros(len(nT)*nState)
    for i in nT:
        F[nState*i:nState*i+nState] = TurbineEqns(t,x[nState*i:nState*i+nState],power_aero[i]);
        
    return F 



def TurbineEqns(t,x,power_yaw_1):
    
    No_states = 11 
    F = np.zeros(No_states,) # Initialization

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
    N = 10;
    
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
    theta_tw = x[9];
    wt = x[10];
    
    R1 = R1/N;
    Ls = Ls/N;
    Kid  = Kid/N;
    Kiq = Kiq /N;
    R2 = R2/N;
    ksh = ksh*N;
    csh = csh*N;
    Ht = Ht*N;
    Hg = Hg*N;

     

    # Mechanical Torque from FLORIS
    Tm = N*power_yaw_1/(5e6*wt);
    Teref = N*Kopt*wg**2;
    
    # Intermediate variables
    Te =  +  (eqs /ws) * iqs + (eds /ws) * ids;

    Qs =  -vqs * ids + vds * iqs;

    iqr  = - eds /1 - Kmrr * iqs;
    
    idr = + eqs /1 - Kmrr * ids;
    
    vqr =  Kiq* Kte * (Teref - Te) + Kiq * Kte /Tte * phiTe \
            - Kiq * iqr + Kiq/Tiq * phiIq; 
    
    vdr = +  Kid* Kqs * (Qsref - Qs) + Kid * Kqs /Tqs * phiQs \
            - Kid * idr + Kid/Tid * phiId; 
    
    
    # Differential equations:
    F[0]= wel / Ls * ( -R1 * iqs + ws * Ls * ids  \
                              + wg/ws * eqs - 1/(Tr * ws) * eds  \
                              - vqs + Kmrr * vqr);
      
    F[1] = wel / Ls * ( -R1 * ids - ws * Ls * iqs +\
                              wg/ws * eds - 1/(Tr * ws) * eqs \
                              - vds + Kmrr * vdr);
    
    F[2] = wel * ws * ( R2 * ids  + (1 - wg/ws) * eds - \
                              1/(Tr * ws) * eqs  - Kmrr * vdr);
    
    F[3] = wel * ws * (- R2 * iqs  - (1 - wg/ws) * eqs - \
                              1/(Tr * ws) * eds  + Kmrr * vqr);
    
    F[4] = Teref - Te;
    
    F[5] = Kte * (Teref - Te)  + Kte /Tte * phiTe  - iqr;

    F[6] = Qsref - Qs;
    
    F[7] = Kqs * (Qsref - Qs)  + Kqs /Tqs * phiQs  - idr;
    
    F[8] = 1/(2*Hg) * (ksh*theta_tw + \
                             csh*wel* (wt - wg) - Te);
    
    F[9] = wel*(wt - wg);
    
    F[10] = 1/(2*Ht) * (Tm - ksh * theta_tw \
                             - csh*wel*(wt - wg));
     
    
    return F
 
    

    
    