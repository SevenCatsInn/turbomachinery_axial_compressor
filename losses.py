import numpy as np
import matplotlib.pyplot as plt
import sympy
import math
from numpy.linalg import norm
from numpy import sqrt, arctan, tan, pi, cos, sin

def   losses(radius,chord,Rm,b,Va_out,Va_in,beta_in,beta_out,alpha_in,alpha_out,Vmag_in,Vmag_out,Wa_in,Wa_out,Wmag_in,Wmag_out,rho_in,rho_out,stagger,Nrow,blades,mdot,Pt_in,p_in,shrouded,stator):

    pitch = 2 * pi * radius / blades
    solidity = chord / pitch #computation of the solidity at each radius
    print(solidity)

    #profile losses

    #equivalent diffusion parameter
    if stator==0:
        ax_vel_densityratio=(rho_out*Wa_out)/(rho_in*Wa_in)
        Wmax_W_in=1.12 + 0.61*cos(beta_in)**2/solidity*(tan(beta_in)-tan(beta_out)) #relation page 131 aungier
        Deq=Wmax_W_in*cos(beta_out)/cos(beta_in)*Wa_in/Wa_out #aungier 6-36
        omega_profile = 2 * solidity / cos(beta_out) * (cos(beta_out)/cos(beta_in)*Wa_in/Wa_out) ** (-2) * 0.004 * ( 1 + 3.1 * (Deq-1)**2 +    0.4*(Deq-1)**8) #aungier page 132
    else: 
        ax_vel_densityratio=(rho_out*Va_out)/(rho_in*Va_in)
        Vmax_V_in=1.12 + 0.61*cos(alpha_in)**2/solidity*(tan(alpha_in)-tan(alpha_out)) #relation page 131 aungier
        Deq=Vmax_V_in*cos(alpha_out)/cos(alpha_in)*Va_in/Va_out #aungier 6-36
        omega_profile = 2 * solidity / cos(alpha_out) * (cos(alpha_out)/cos(alpha_in)*Va_in/Va_out) ** (-2) * 0.004 * ( 1 + 3.1 * (Deq-1)**2 +    0.4*(Deq-1)**8) #aungier page 132
    



    #tip leakage losses (from reference book)
    delta_c=0.001 #clearance [m]
    rtip=Rm+b/2-delta_c # [m]
    r_in=rtip # [m]
    r_out=r_in
    if shrouded==0:
        rho_med=(rho_in+rho_out)*0.5
        Vt_in=Va_in*tan(alpha_in)
        Vt_out=Va_out*tan(alpha_out)
        tau=pi*delta_c*(rho_in*r_in*Va_in+rho_out*r_out*Va_out)*(r_out*Vt_out-r_in*Vt_in) #in reference book C_theta and C_m are the absolute t tangential and meridional comp
        deltaP=abs(tau/(blades*rtip*delta_c*chord*cos(stagger))) #blades is the blade number
        U_c=0.816*sqrt(2*deltaP/rho_med)/Nrow**(0.2)
        mdot_c=rho_med*U_c*blades*delta_c*chord*cos(stagger)
        deltaP_t=deltaP*mdot_c/mdot
        omega_tip=deltaP_t/(Pt_in-p_in);
    else: 
        omega_tip=0

    omega_overall = omega_profile+omega_tip


    return omega_overall


#omega_ew=(Pt1_p1-p2_p1*Pt2_p2)/(Pt1_p1-1)
