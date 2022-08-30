import numpy as np
import matplotlib.pyplot as plt
import sympy
import math
from numpy.linalg import norm
from numpy import sqrt, arctan, tan, pi, cos, sin

def off_design(Rm,mdot_off,beta,rho_out,Um,alpha,rho_in,gamma,efficiency_TT,cp,Tt_in,b):
    S_in=np.pi*2*Rm*b
    S_out=S_in
    alpha = alpha*np.pi/180
    Leul_off = Um * (Um + mdot_off * tan(beta) /( rho_out * S_out )) - Um * mdot_off * tan(alpha) / (rho_in * S_in)
    beta_off = (Leul_off * efficiency_TT / (cp * Tt_in) + 1 ) ** (( gamma )/( gamma - 1 ))
    Um=Um
    return Leul_off,beta_off
