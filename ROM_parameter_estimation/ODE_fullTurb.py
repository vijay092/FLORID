# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 14:41:53 2020

@author: sanja
"""


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
    Xm = 2*np.pi*60*Lm;
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
    ild = x[9];
    ilq = x[10];
    iod = x[11];
    ioq = x[12];
    vod = x[13];
    voq = x[14];
    gammad = x[15];
    gammaq = x[16];
    P_avg = x[17];
    Q_avg = x[18];
    phi_P = x[19];
    phi_Q = x[20];
    vbdf = x[21];
    phi_PLL = x[22];
    delta_turb = x[23];
    
    R1 = R1/N;
    Ls = Ls/N;
    Kid  = Kid/N;
    Kiq = Kiq /N;
    R2 = R2/N;
    ksh = ksh*N;
    csh = csh*N;
    Ht = Ht*N;
    Hg = Hg*N;
    Xm = Xm/N;
    

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

    iqr  = - eds/Xm - Kmrr * iqs;
    
    idr =  eqs/Xm - Kmrr * ids;
    
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
    
    # inverter dynamics

    omega_nom = 1;
    
    omega_c = 50.26; 
    omega_c_pll = 2*np.pi*200; 
    
    
    # rating and base
    V_inv_rating = 24e3;
    S_m_rating = 5e6;
    P_base = S_m_rating;
    V_base = V_inv_rating/np.sqrt(3)*np.sqrt(2);
    omega_base = 2*np.pi*60;
    
    # scalings
    kappa_p = P_base/1000;
    kappa_v = V_inv_rating/(120*np.sqrt(3));
    
    # the number of inverters
    L_base = (V_base**2)/P_base;
    Z_base = (V_base**2)/P_base;
    C_base = P_base/V_base**2;
    
    PIpll_base = omega_base/V_base;
    PIp_base = 1/V_base;
    PIi_base = (V_base**2)/P_base;
    
    
    P_ref = (vdr*idr + vqr*iqr); 
    Q_ref = (-vqr*idr + vdr*iqr);
    

    vg_d = np.cos((- delta_turb)*omega_base)*vds \
    - np.sin(( - delta_turb)*omega_base)*vqs;
    vg_q = np.sin(( - delta_turb)*omega_base)*vds \
    + np.cos(( - delta_turb)*omega_base)*vqs;
    
    
    # power controllers
    kp_q = 0.01/kappa_v/PIp_base;
    ki_q = 0.1/kappa_v/PIp_base;
    kp_p = 0.01/kappa_v/PIp_base;
    ki_p = 0.1/kappa_v/PIp_base;
    
    # PLL
    kp_pll = 0.25/kappa_v/PIpll_base;
    ki_pll = 2/kappa_v/PIpll_base;
    
    # LCL-filter parameters
    Lf = 1*1e-3/kappa_p*kappa_v**2/L_base/N;
    rf = 0.7/kappa_p*kappa_v**2/Z_base/N;
    Cf = 24*1e-6*kappa_p/kappa_v**2/C_base*N;
    Rd = 0.02/kappa_p*kappa_v**2/Z_base/N;
    Lc = 200*1e-6/kappa_p*kappa_v**2/L_base/N;
    rc = 0.12/kappa_p*kappa_v**2/Z_base/N;
    
    # current controllers
    kp_id = 6/kappa_p*kappa_v**2/PIi_base/N;
    ki_id = 350/kappa_p*kappa_v**2/PIi_base/N;
    kp_iq = 6/kappa_p*kappa_v**2/PIi_base/N;
    ki_iq = 350/kappa_p*kappa_v**2/PIi_base/N;
    
    # reference for i_l
    ild_ref = ki_q*phi_Q + kp_q*(Q_ref - Q_avg);
    ilq_ref = ki_p*phi_P + kp_p*(P_ref - P_avg);
    
    # omega_pll = (omega_nom - kp_pll*vbdf + ki_pll*phi_PLL)*omega_base;
    omega_pll = omega_nom - kp_pll*vbdf + ki_pll*phi_PLL ;
    
    # assume no loss in the inverter: v_i = v_i_ref
    vid = kp_id*(ild_ref - ild) \
    + ki_id*gammad - omega_pll*Lf*ilq;
    
    viq = kp_iq*(ilq_ref - ilq) \
    + ki_iq*gammaq + omega_pll*Lf*ild;
    
    f_9 = 1/Lf*(vid - vod - rf*ild) + omega_pll*ilq; # ild
    
    f_10 = 1/Lf*(viq - voq - rf*ilq) - omega_pll*ild; # ilq
    
    f_11 = 1/Lc*(vod - vg_d - rc*iod) + omega_pll*ioq; # iod
    
    f_12 = 1/Lc*(voq - vg_q - rc*ioq) - omega_pll*iod; # ioq
    
    f_13 = Rd*(f_9  - f_11) \
               + 1/Cf*(ild - iod) + omega_pll*voq \
               - omega_pll*Rd*(ild - iod); # vod
               
    f_14 = Rd*(f_10 - f_12) +\
               1/Cf*(ilq - ioq) - omega_pll*vod \
               + omega_pll*Rd*(ilq - ioq); # voq
    
    f_15 = ild_ref - ild; # gamma_d
    
    f_16 = ilq_ref - ilq; # gamma_q
    
    # instantaneous power
    p = (vg_d*iod + vg_q*ioq);
    q = (vg_q*iod - vg_d*ioq);
    
    f_17 = omega_c*(p - P_avg); # P_avg
    
    f_18 = omega_c*(q - Q_avg); # Q_avg
    
    f_19 = P_ref - P_avg; # phi_P
    
    f_20 = Q_ref - Q_avg; # phi_Q
    
    f_21 = omega_c_pll*vg_d - omega_c_pll*vbdf; # vbd,f
    
    f_22 = -vbdf; # phi_PLL
    
    f_23 =  (- kp_pll*vbdf + ki_pll*phi_PLL) ; # delta
    
    
    return f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, \
            f_10, f_11, f_12, f_13, f_14, f_15, f_16, f_17, f_18, f_19,\
            f_20, f_21, f_22, f_23

    
     

    