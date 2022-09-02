import numpy as np
import matplotlib.pyplot as plt
import sympy
import math
from numpy.linalg import norm
from numpy import sqrt, arctan, tan, pi, cos, sin

def off_design(Rm,mdot_off,beta,rho_out,Um,alpha,rho_in,gamma,efficiency_TT,cp,Tt_in,b,Pt_in, n_blades, th, chord, theta):

    sigma = 2 * np.pi * Rm / n_blades

    delta0  = 0.01 * sigma  * abs(beta)  + (0.74*sigma**1.9 + 3 * sigma) * (abs(beta) / 90) ** (1.67 + 1.09 * sigma)

    kdeltath  = 6.25 * (th/chord) + 37.5 * (th/chord)**2

    b1  = 0.9625 - 0.17 * abs(beta)/100 - 0.85 * (abs(beta)/100)**3

    m1  = 0.17 - 0.0333 * abs(beta) /100 + 0.333 * (beta/100) **2

    m  = m1  / (sigma  ** b1 )

    delta = delta0 * kdeltath + m * theta / (sigma ** b1)


    S_in=np.pi*2*Rm*b
    S_out=S_in
    alpha = alpha*np.pi/180
    R=gamma/(gamma-1)*cp
    err=10
    
    

    while abs(err) > 0.001:
        Leul_off = Um * (Um + mdot_off * tan(beta) /( rho_out * S_out )) - Um * mdot_off * tan(alpha) / (rho_in * S_in)
        beta_off = (Leul_off * efficiency_TT / (cp * Tt_in) + 1 ) ** (( gamma )/( gamma - 1 ))
        Tt_out = Leul_off/cp + Tt_in
        V_out = sqrt((mdot_off/( rho_out * S_out ))**2 + (Um + mdot_off * tan(beta) /( rho_out * S_out ))**2)
        T_out= Tt_out - V_out**2/(2*cp)
        M_out=sqrt(gamma*R*T_out)
        Pt_out=Pt_in*beta_off
        P_out=Pt_out*(1+(gamma-1)/2*M_out**2)**(gamma/(gamma-1))
        rho_out_new=P_out/(R*T_out)
        err=rho_out_new-rho_out
        rho_out=rho_out_new
        
        Um=Um
        return Leul_off,beta_off
