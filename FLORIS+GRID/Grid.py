# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 09:46:34 2019

@author: sanja
"""


import numpy as np
import scipy.linalg as la
import ODEs

No_mach = 4;
No_turb = 4;
No_states_mach = 4;
No_states_turb = 11;
R_L = 5;
X_L = 0.5;


def MatrixArrange(a,b,c):
    """
    
    More often than not in Power Systems, we need this matrix representation
    (Recall Rotational Matrices)
    
    D = [a  -b;
         c   a]
    
    """

    left = np.vstack((a, c))
    right = np.vstack((-b, a))
    
    return np.hstack((left, right))


###############################################################################
    
def TurbineCurrentsInGlobal(t,x):
    
    # These are the currents available at the turbine
    i_s_q = x[0::No_states_turb];
    i_s_d = x[1::No_states_turb];
    
    #### ? ? ###### Add the rotor currents, although its small ########

    # The idea is to add turbine only where it is located
    Iturb = np.zeros(2*No_mach);
    Iturb[No_mach -1] = np.sum(i_s_d)
    Iturb[-1] = np.sum(i_s_q)
    
    return Iturb
    
    
    
###############################################################################    
def GridInterconnection(Ybus,x_sync,Iturb):
    """
    
    This function computes voltages of various nodes using the admittance 
    matrices following a number of numerical calculations
    
    Inputs:
        Ybus := Admittance matrix
    Outputs:
        Voltages of every bus
        
    
    """
    
    # Details of the network
    # load_idx = [6, 8];   # Indices where loads are present
    # Ybus[load_idx,load_idx] = Ybus[load_idx,load_idx] \
    #                           + complex(R_L,X_L)/(R_L*complex(0,X_L)); 
                             
    Y_gg = Ybus[:No_mach,:No_mach];
    Y_gl = Ybus[:No_mach,No_mach:];
    Y_ll = Ybus[No_mach:,No_mach:];

    Y_kr = Y_gg - np.dot(np.dot(Y_gl,la.inv(Y_ll)),np.transpose(Y_gl))

    Y_kr_pr = MatrixArrange(Y_kr.real, Y_kr.imag, Y_kr.imag)
    # Parameters of the synchronous machine
    X_d_pr = 0.031;
    X_q_pr = X_d_pr;
    R_s = 0;

    # Machine impedance matrix
    Z_S = MatrixArrange(R_s,X_q_pr, X_d_pr)
    Z_S_aug = la.block_diag(Z_S, Z_S, Z_S, Z_S) # Coz there are 4 machines
    
    # Transformation matrix
    T = []
    for i in range (No_mach):
        theta = x_sync[i*No_states_mach + 2]
        T_cur = MatrixArrange(np.sin(theta), np.cos(theta),np.cos(theta))
        T.append(T_cur)
    T_aug = la.block_diag(*T)
    
    # Arrange as DDQQ instead of DQDQ
    even = np.arange(0,len(Z_S_aug),2)
    odd = np.arange(1,len(Z_S_aug),2)
    new_order = np.hstack((even, odd))
    s = Z_S_aug[:,new_order]
    Z_S_DQ = s[new_order,:];
    
    T_bleh = T_aug[:,new_order]     # Transformation matrix rearranged
    T_DQ = T_bleh[new_order,:]
    
    # Terminal voltages of the synchronus machine
    Ed = x_sync[0::4]  # Collect d-voltages of all the machines
    Eq = x_sync[1::4]  # Collect q-voltages of all the machines
    Edq = np.transpose(np.hstack((Ed,Eq)))
    
    # Compute currents with this info:
    bleh = np.dot(np.dot(T_DQ,la.inv(Y_kr_pr)),np.transpose(T_DQ))
    bleh1 = la.inv(Z_S_DQ + bleh)
    bleh2 = Edq - np.dot(np.dot(T_DQ,la.inv(Y_kr_pr)), Iturb)
    I_DQ = np.dot(bleh1,bleh2)
    Id = I_DQ[:No_mach];
    Iq = I_DQ[No_mach:];
    
    
    # Terminal voltages:
    Vdq = Edq - np.dot(Z_S_DQ,I_DQ)
    Vd = Vdq[:No_mach]
    Vq = Vdq[No_mach:]

    
    return Vd, Vq, Id, Iq
    

###############################################################################    
    


    
def NonLinSystem(t,x, power_yaw, Ybus):
    """
    
    Objective: 
        Describes the system of equations for a coupled machine-turbine system
        dictated by a given network. We use the Kundur 2-area system as a 
        test-case
        
    Inputs:
        x := state vector
        t := time
        fi := struct containing flow info
        Ybus := Admittance matrix
        
    Outputs:
        F := d/dts
        
    Note:
        In this function, we will work with 2 systems
        1) Turbine := No_mach*No_states_mach : No_mach*No_states_mach + 11
        2) Synchronous Machines := 0:No_mach*No_states_mach
        
        The network gives us voltage info of the buses and we use this to
        simulate the dynamics of turbine and synch mach
        
    """

    # So now we have things written in different functions 
    # This is sorta the main function where we get the ODEs
    
    # Initialization:
    Nfull = No_mach*No_states_mach + No_states_turb*No_turb
    
    #==========================================================================
    # 0) Get the currents from the turbine to feed into the grid
    
    Iturb = TurbineCurrentsInGlobal(t,x[No_mach*No_states_mach:])
    
    #==========================================================================
    # 1) Get voltage from the network
    
    Vd, Vq, Id, Iq = GridInterconnection(Ybus,x[:No_mach*No_states_mach], Iturb)
    #==========================================================================
    # 2) Use this voltage to simulate SM
    
    omega_slack = x[3]
    DXsyncDT = np.zeros(No_mach*No_states_mach)
    
    # Loop over machines
    l = 0;
    for mach in range(No_mach):
        if mach == 0:
            omega_slack_1 = 377;
        else:
            omega_slack_1 = x[3];
            
        mach_dyn = ODEs.SyncMachineEqns(t, x[:No_mach*No_states_mach],\
                                     omega_slack, omega_slack_1,\
                                     Id[mach], Iq[mach]);
        DXsyncDT[l:l + No_states_mach] = mach_dyn
        l = l + No_states_mach
    #==========================================================================
    # 3) Use this voltage to simulate turbine
    
    DXturbDT =np.zeros(No_states_turb*No_turb)

    # Loop over turbines
    k = 0;
    for turb in range(No_turb):
        turb_dyn = ODEs.TurbineEqns(t, x[No_mach*No_states_mach:],\
                                     power_yaw[turb]/5e6, Vd[-1],Vq[-1]);
        DXturbDT[k:k + No_states_turb] = turb_dyn
        k = k + No_states_turb
    #==========================================================================
    # 4) Stack 'em
    F = np.hstack((DXsyncDT,DXturbDT))
    
    return F