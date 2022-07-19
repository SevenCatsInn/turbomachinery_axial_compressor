# exec(open("./MLD.py").read())

from numpy import * 

# Thermophysical properties
c_p = 1005  # Constant pressure specific heat [J/(kg K)]
gamma = 1.4 # Specific heat ratio
c_v = c_p/gamma
R = c_p * (gamma-1)/gamma # Gas constant [J/(kg K)]

# Geometry
R_m = 0.3

T_t1 = 300
p_t1 = 100000
beta = 1.2
gamma = 1.4
m_dot_req=100

#Non dimensional coefficients
chi = 0.5
phi = 0.5
psi = 0.3
lamb = 2 * psi
eta = 0.92


L_is = c_p * T_t1 * ( beta**((gamma-1)/gamma) - 1)
L_eul = L_is / eta

U_m = sqrt(L_eul/psi)
omega=U_m/R_m

V_a1m = phi * U_m

V_tinfm = (1 - chi) * U_m

V_t1m = U_m * (1 - chi - lamb/4 )

alpha_1= arctan(V_t1m/V_a1m)

V_1m = sqrt( V_a1m**2 + V_t1m**2 )

T_1m = T_t1 - V_1m**2 / (2*c_p)
M_1m = V_1m / sqrt(gamma*R*T_1m)
p_1m = p_t1 * ( 1 + (gamma-1) * M_1m**2 / 2 ) ** (- gamma/(gamma-1))
rho_1m= p_1m/(R*T_1m)


b_1 = m_dot_req / (rho_1m * V_a1m * 2 * pi * R_m)


### Outlet station

eta_R = 0.92

V_t2m = L_eul / U_m + V_t1m
T_t2 = L_eul / c_p + T_t1

V_a2m = V_a1m 

err = 1e10
eps=0.01

while err>eps:
    
    V_a2m_old = V_a2m

    V_2m = sqrt(V_a2m**2 + V_t2m**2)
    T_2m = T_t2 - V_2m**2/(2*c_p)
    T_2mis = eta_R * (T_2m - T_1m) + T_1m
    p_2m = p_t1 * (T_2mis/T_t1)**(gamma/(gamma-1))
    rho_2m = p_2m / (R*T_2m)
    
    V_a2m = m_dot_req / (rho_2m * 2 * pi * R_m * b_1)

    err = abs (V_a2m-V_a2m_old) 
    


chi = (V_a1m**2 + V_t1m**2 - V_a2m**2 - V_t2m**2 + 2*U_m*(V_t2m - V_t1m)) / (2*L_eul)

print("V_a1m = " , V_a1m)
print("V_t1m = " , V_t1m)
print("V_a2m = " , V_a2m)
print("V_t2m = " , V_t2m)

print("T_1m = ", T_1m)
print("p_1m = ", p_1m)
print("T_2m = ", T_2m)
print("p_2m = ", p_2m)

