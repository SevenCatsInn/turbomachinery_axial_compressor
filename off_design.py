from sympy import Symbol, solve
import numpy as np
from numpy.linalg import norm

exec(open("./turboproject.py").read())

mdot=90 #[kg/s]
#from the mass flow rate and the values of total pressure and density it is possible to derive the Mach number
#mdot = pt1 * ( np.pi * (r_tip ** 2 - r_hub **2 ) / sqrt ( R * Tt1) * g_m_gamma  

g_m_gamma = mdot * np.sqrt ( R * Tt1) / (Pt1 * ( np.pi * (rtip ** 2 - rhub **2 )))


M1off = Symbol('M1off')
eq1 = g_m_gamma - np.sqrt (gamma) * M1off * (1 + ( gamma - 1 ) / 2 * M1off ** 2) ** ((-gamma - 1 ) / ( 2 * (gamma - 1 )))
sol_M1 = solve (eq1 , M1off)
sol_M1_np = np.array(sol_M1[0], dtype=np.float64)
M1_off = sol_M1_np
T1_off = Tt1 / (1 + (gamma-1) / 2 * M1_off ** 2 )
P1_off = Pt1 / ((1 + (gamma-1) / 2 * M1_off ** 2 ) ** ( gamma / ( gamma - 1 )))
rho1_off = P1_off / (R*T1_off)
Va0_off = mdot / ( rho1_off * np.pi * (rtip ** 2 - rhub ** 2 )) #get the axial velocity

#from the new value of inlet axial velocity we can get the new velocity triangles and therefore the new work
#velocity triangles are computed on the mean line

#on deflector
#given geometric deflection of first stator alpha given by atan(vt/va) at design conditions
#axial component remains the same, tangential changes

#rotor inlet velocity
V1_mag_off= Va0_off / cos (alpha1)
Va1_off = Va0_off
Vt1_off = V1_mag_off * sin (alpha1)
V1_off = [ Va1_off, Vt1_off]

Wt1_off = Vt1_off - Um
Wa1_off = Va1_off
W1_off = [ Wa1_off, Wt1_off ]
W1_mag_off = norm(W1_off) 

#rotor outlet velocity
#outlet geometric angle of the rotor blade is needed alpha 2

Wa2_off = Wa1_off
W2_mag_off = Wa2_off / cos (beta2-beta1)
Wt2_off = W2_mag_off * sin (beta2-beta1)
W2_off = [ Wa2_off, Wt2_off ]

Va2_off = Wa2_off
Vt2_off = Wt2_off + Um
V2_off = [ Va2_off, Vt2_off ]

#stator outlet velocity 
# alpha3 is the outlet geometric angle of the stator
Va3_off = Wa2_off
V3_mag_off = Va3_off / cos (alpha3-alpha2)  
Vt3_off = V3_mag_off * sin (alpha3-alpha2)
V3_off = [ Va3_off, Vt3_off ]

Leul_off = Um * (Vt2_off - Vt1_off)
#Leul_off = cp * (Tt2 - Tt1)
Tt2_off = Leul_off / cp + Tt1
beta_off = (Tt2_off / Tt1 ) ** (gamma / (gamma - 1 ))




