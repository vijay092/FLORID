# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:58:23 2019

@author: sanja
"""

import numpy as np
import scipy.linalg as la

def TurbineEqns(t,x,power_yaw_1,vds,vqs):
    
    No_states = 11 
    F = np.zeros(No_states,) # Initialization
    
    ### Parameters of the turbine 
    
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
    #Xm = 2*np.pi*60*Lm;

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
     

    
    
    # Mechanical Torque from FLORIS
    Tm = power_yaw_1/wt;
    Teref = Kopt*wg**2;
    
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
 
    
def SyncMachineEqns(t,x, omega_slack,omega_slack_1, Id,Iq):
    
    No_states = 4
    F = np.zeros(No_states)
    
    # Parameters:
    T_do_pr = 9.8 ;
    X_d = 0.069;
    X_d_pr = 0.031;
    T_qo_pr = 1;
    X_q = 0.1;
    X_q_pr = X_d_pr;
    D = 6;
    M = 2*8/(2*np.pi*60);
    Pm = 0.001;
    Efd = 1;
    
    # Rename the states for ease of analysis
    Ed = x[0];
    Eq = x[1];
    delta = x[2];
    omega = x[3];
     
    # Electrical Power
    Pe = Eq*Iq + Ed*Id -(X_d_pr - X_q_pr)*Id*Iq;
    # Differential equations
    #print("id =", Id)
    F[0] =  1/T_qo_pr * (-Ed + (X_q - X_q_pr)*Iq); 
    F[1] = 1/T_do_pr * (-Eq - (X_d - X_d_pr)*Id + Efd); 
    F[2] = omega - omega_slack;
    F[3] = 1/M * (Pm  - Pe - D*(omega - omega_slack_1));
    
    return F
    
    