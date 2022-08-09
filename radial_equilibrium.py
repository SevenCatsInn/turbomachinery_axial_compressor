### Radial equilibrium script ###

exec(open("./turboproject.py").read()) # Run mean line design

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"text.usetex": True})

def finDiff(x,deltaX):
    # Finite differences function over a list
    dx = []
    [dx.append( (x[i+1] - x[i-1]) / (2*deltaX) ) for i in range(1,len(x) - 1)]
    dx = [(x[1] - x[0]) / deltaX] + dx + [(x[-1] - x[-2]) / deltaX]

    return dx

arrayLst = lambda x: np.array(list(x))

# Geometry
R_h = R_m - b_1 / 2  # Hub Radius          [m]   
R_t = R_m + b_1 / 2  # Tip Radius          [m]  

# Discretization
pts = 50  # Total number of points across the radius, 
if pts % 2 == 0: pts = pts + 1 # Make pts uneven if it's even

rr = np.linspace(R_h, R_t, pts) # Discrete space of the radii over which we compute our quantities
deltaR = (R_t - R_h)/ (pts - 1) # Radius interval between points
mean_index = pts//2  # Index of the various lists corresponding to mean radius quantities

omega = 2 * pi * rpm/60                          # Angular velocity     [rad/s]
U = arrayLst( omega * rr[t] for t in range(pts)) # Peripheral velocity  [m/s]

# Input data
T_t1 = 300 * np.ones(pts) # [K]     --> f(r)
p_t1 = 1e5 * np.ones(pts) # [Pa] --> f(r)
s_1 = np.zeros(pts)
ds_1 = np.zeros(pts)
m_dot_req = 100 # Required mass flow [kg/s]
T_1 = T_1m * np.ones(pts)

# Computed Quantities
h_t1 = c_p * T_t1 # Total enthalpy [J/kg]
dh_t1 = finDiff(h_t1,deltaR)

for j in range(pts):
    ds_1[j]  = -R * finDiff(p_t1,deltaR)[j] / p_t1[j] + c_p * finDiff(T_t1,deltaR)[j] / T_t1[j] # Derivative over r of entropy [J/(kg K)]





















# Set the design choice for tangential velocity distribution in the radial direction

# First power vortex distribution
# V_t1 = arrayLst( V_t1m * R_m / rr[t] for t in range(pts)) # Free Vortex
n = 1
V_t1 = arrayLst( a * rr[t]**n - b / rr[t] for t in range(pts)) # Power Design


rV_t1 = arrayLst(rr[t] * V_t1[t] for t in range(pts))
drV_t1 = finDiff(rV_t1, deltaR)

print("")
print("########## ROTOR INLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5
iter= 1

# This loop can be eliminated by varying b_1 to accomodate for m_dot_req
while abs(err) > tol: 
    
    V_a1 = np.zeros(pts) # Create the list
    V_a1[mean_index] = V_a1m  # Initial value for forward integration starting from mean radius
    dV_a1 = np.zeros(pts)

    # N.I.S.R.E. 1 numerical integration
    # --> Start from V_1m at R_m and move forwards and backwards up to R_t and R_h
    # j moves from 1 to the mean_index
    # q and k are a subloop to simplify the code, the first values of q,k corresponding to
    # the "forwards" integration, and the second values to the "backwards" integration
    
    for j in list(range(0,mean_index)):
        for q,k in zip([mean_index + j, mean_index - j],[1,-1]):
            dV_a1[q] = 1 / V_a1[q] * ( dh_t1[q] - T_1[q] * ds_1[q] - V_t1[q] / rr[q] * drV_t1[q] )
            V_a1[q + k*1] = V_a1[q] + dV_a1[q] * k * deltaR 
        
        
    # Initiate all the lists
    V_1 , alpha_1, W_t1, W_a1, W_1, beta_1, p_1, rho_1, M_1, M_1r, p_t1, p_t1r, integrand_1, chi, L_eul = (list(np.zeros(pts)) for t in range(15))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        V_1[j] = np.sqrt(V_a1[j]**2 + V_t1[j]**2)
        alpha_1[j] = np.arctan(V_t1[j]/V_a1[j])
        W_t1[j] = V_t1[j] - U[j]
        W_a1[j] = V_a1[j]
        W_1[j] = np.sqrt(W_t1[j]**2 + W_a1[j]**2)
        beta_1[j] = np.arctan(W_t1[j]/W_a1[j])

        # Thermodynamics
        T_1[j] = T_t1[j] - V_1[j]**2 / (2 * c_p)
        p_1[j] = p_1m * (T_1[j] / T_1m)**(gamma/(gamma-1))
        rho_1[j] = p_1[j] / (R*T_1[j])
        M_1[j]  = V_1[j] / np.sqrt(gamma * R * T_1[j])
        M_1r[j] = W_1[j] / np.sqrt(gamma * R * T_1[j])
        p_t1[j]  = p_1[j]*(1 + (gamma-1) / 2 * M_1[j]**2  ) ** (gamma/(gamma-1))
        p_t1r[j] = p_1[j]*(1 + (gamma-1) / 2 * M_1r[j]**2 ) ** (gamma/(gamma-1))
        
        integrand_1[j] = 2 * np.pi * rr[j] * rho_1[j] * V_a1[j] 

    m_dot_trap = np.trapz(integrand_1, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a1m = V_a1m*(1 + err) # New axial velocity
    
    
    print("")
    print("---Iteration no. " + str(iter))
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("err = "+ str(err))
    print("V_a1m = "+ str(V_a1m) + " [m/s]")
    iter += 1























print("")
print("########## ROTOR OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5 # Tolerance of error wrt the desires mass flow value
iter = 1


# Inputs

# Entropy inputs, NOTE: absolute values are meaningless
omega_loss_R = 0.08 # Coefficient of loss

# Need to transform s_2 and ds_2 into lists otherwise numpy will assign the same id to s_1 and s_2, even with s_2 = s_1[:] why??
s_2  = list( s_1)    # Initial radial entropy distribution in 2
ds_2 = list(ds_1) # Dertivative wrt r of entropy

# V_t2 = arrayLst( V_t2m * R_m / rr[t] for t in range(pts)) # Free Vortex 
V_t2 = arrayLst( a * rr[t]**n + b / rr[t] for t in range(pts)) # Power Design


rV_t2  = arrayLst(rr[t] * V_t2[t] for t in range(pts))
drV_t2 = finDiff(rV_t2,deltaR)

h_t2 = arrayLst(h_t1[t] + U[t] * (V_t2[t] - V_t1[t]) for t in range(pts)) # Total enthalpy in 2 
T_t2 = h_t2 / c_p # Total temperature

T_2 = T_2m * np.ones(pts) # Static temperature

dh_t2 = finDiff(h_t2,deltaR)


# This loop can be avoided using flaired blades b_2 != b_1
while abs(err) > tol: # Begin loop to get mass flow convergence

    V_a2 = list(np.zeros(pts)) # Create the list
    V_a2[mean_index] = V_a2m  # Initial value for forward integration starting from mean radius
    dV_a2 = list(np.zeros(pts))

    # N.I.S.R.E. 2 numerical integration 
    for j in list(range(0,mean_index)):
        for q,k in zip([mean_index + j, mean_index - j],[1,-1]):
            dV_a2[q] = 1 / V_a2[q] * ( dh_t2[q] - T_2[q] * ds_2[q] - V_t2[q] / rr[q] * drV_t2[q] )
            V_a2[q + k*1] = V_a2[q] + dV_a2[q] * k * deltaR 
        
        
    # Initiate all the lists
    V_2 , alpha_2, W_t2, W_a2, W_2, beta_2, p_2, rho_2, M_2, M_2r, p_t2, p_t2r, integrand_2, chi, L_eul = (list(np.zeros(pts)) for t in range(15))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        V_2[j] = np.sqrt(V_a2[j]**2 + V_t2[j]**2)
        alpha_2[j] = np.arctan(V_t2[j]/V_a2[j])
        W_t2[j] = V_t2[j] - U[j]
        W_a2[j] = V_a2[j]
        W_2[j] = np.sqrt(W_t2[j]**2 + W_a2[j]**2)
        beta_2[j] = np.arctan(W_t2[j]/W_a2[j])

        # Thermodynamics
        T_2[j] = T_t2[j] - V_2[j]**2 / (2 * c_p)
        p_2[j] = p_2m * (T_2[j] / T_2m)**(gamma/(gamma-1)) * np.exp(- (s_2[j] - s_2[mean_index]) / R)
        rho_2[j] = p_2[j] / (R*T_2[j])
        M_2[j]  = V_2[j] / np.sqrt(gamma * R * T_2[j])
        M_2r[j] = W_2[j] / np.sqrt(gamma * R * T_2[j])
        p_t2[j] = p_2[j]*(1 + (gamma-1) / 2 * M_2[j]**2 ) ** (gamma/(gamma-1))


        L_eul[j] = U[j] * (V_t2[j] - V_t1[j])
        chi[j] = (W_1[j]**2 - W_2[j]**2) / (2 * L_eul[j])

        integrand_2[j] = 2 * np.pi * rr[j] * rho_2[j] * V_a2[j] 
                

        p_t2r[j] = p_t1r[j] - omega_loss_R * (p_t1r[j] - p_1[j])

        # ENTROPY EVALUATION

        s_2[j]  = s_1[j] - R * np.log(p_t2r[j] / p_t1r[j])

    ds_2 = finDiff(s_2,deltaR)

    m_dot_trap = np.trapz(integrand_2, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a2m = V_a2m*(1 + err) # New axial velocity
    
    
    print("")
    print("---Iteration no. " + str(iter))
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a2m = " + str(V_a2m))
    print("err = "+ str(err))
    iter += 1













exec(open("./turboproject_S2.py").read()) # Mean line design for 2nd stage





print("")
print("########## STATOR OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5 # Tolerance of error wrt the desires mass flow value
iter = 1

# Input data
omega_loss_S = 0.08

V_t3 = arrayLst( a22 * rr[t]**n - b22 / rr[t] for t in range(pts))

# V_t3 = arrayLst( V_t3m * R_m / rr[t] for t in range(pts))

rV_t3  = arrayLst(rr[t] * V_t3[t] for t in range(pts))
drV_t3 = finDiff(rV_t3,deltaR)

# Initial assumptions
T_3  = list(T_3m * np.ones(pts))
s_3  = s_2[:]
ds_3 = s_2[:]

# Imposed by thermodynamics
h_t3 = h_t2
dh_t3 = dh_t2
T_t3 = T_t2

# This loop can be avoided using flaired blades b_2 != b_1
while abs(err) > tol: # Begin loop to get mass flow convergence
    print("")
    print("---Iteration no. " + str(iter))

    V_a3 = list(np.zeros(pts)) # Create the list
    V_a3[mean_index] = V_a3m  # Initial value for forward integration starting from mean radius
    dV_a3 = list(np.zeros(pts))
    
    # N.I.S.R.E at stator outlet (3)
    for j in list(range(0,mean_index)):
        for q,k in zip([mean_index + j, mean_index - j],[1,-1]):
            dV_a3[q] = 1 / V_a3[q] * ( dh_t3[q] - T_3[q] * ds_3[q] - V_t3[q] / rr[q] * drV_t3[q] )
            V_a3[q + k*1] = V_a3[q] + dV_a3[q] * k * deltaR 

    # Initiate all the lists
    V_3 , alpha_3, p_3, rho_3, M_3, p_t3, integrand_3, W_t3, W_a3, W_3, beta_3, M_3r, p_t3, p_t3r = (list(np.zeros(pts)) for t in range(14))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        alpha_3[j] = np.arctan(V_t3[j]/V_a3[j])
        V_3[j] = np.sqrt(float(V_a3[j]**2 + V_t3[j]**2))
        W_t3[j] = V_t3[j] - U[j]
        W_a3[j] = V_a3[j]
        W_3[j] = np.sqrt(W_t3[j]**2 + W_a3[j]**2)
        beta_3[j] = np.arctan(W_t3[j]/W_a3[j])
        
        # Thermodynamics
        T_3[j] = T_t3[j] - V_3[j]**2 / (2 * c_p)
        p_3[j] = p_3m * (T_3[j] / T_3m)**(gamma/(gamma-1)) * np.exp(- (s_3[j] - s_3[mean_index]) / R)
        rho_3[j] = p_3[j] / (R*T_3[j])
        M_3[j]   = V_3[j] / np.sqrt(gamma * R * T_3[j])
        M_3r[j]  = W_3[j] / np.sqrt(gamma * R * T_3[j])
        p_t3[j]  = p_3[j]*(1 + (gamma-1) / 2 * M_3[j]**2  ) ** (gamma/(gamma-1))
        p_t3r[j] = p_3[j]*(1 + (gamma-1) / 2 * M_3r[j]**2 ) ** (gamma/(gamma-1))

        integrand_3[j] = 2 * np.pi * rr[j] * rho_3[j] * V_a3[j]
        
        #Evaluate the q.ties in section 1 (np.expressions) at the current radius
        # tmp = overwritten at every iteration, no need for a new array for _1 quantities
        
        p_t3[j] = p_t2[j] - omega_loss_S * (p_t2[j] - p_2[j])

        # ENTROPY EVALUATION

        s_3[j]  = s_2[j] - R * np.log(p_t3[j] / p_t2[j])

    ds_3 = finDiff(s_3,deltaR) # Derivative of s_3

    m_dot_trap = np.trapz(integrand_3, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a3m = V_a3m*(1 + err) # New axial velocity
    Vt3_check = V_t3[mean_index]
    

    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a3m = "+ str(V_a3m) + " [m/s]")
    print("V_t3m = "+ str(Vt3_check) + " [m/s]")
    print("err = "+ str(err))
    iter += 1

























print("")
print("########## STAGE 2 ROTOR OUTLET ##########")

# Geometry
R_h2 = R_m - b_2 / 2  # Hub Radius          [m]   
R_t2 = R_m + b_2 / 2  # Tip Radius          [m]

rr2 = np.linspace(R_h2, R_t2, pts) # Discrete space of the radii over which we compute our quantities
deltaR2 = (R_t2 - R_h2)/ (pts - 1) # Radius interval between points
mean_index = pts//2  # Index of the various lists corresponding to mean radius quantities

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5 # Tolerance of error wrt the desires mass flow value
iter = 1


# Inputs

# Entropy inputs, NOTE: absolute values are meaningless
omega_loss_R = 0.08 # Coefficient of loss

# Need to transform s_4 and ds_4 into lists otherwise numpy will assign the same id to s     and s_4, even with s_4 = s_3[:] why??
s_4  = list( s_3)    # Initial radial entropy distribution in 2
ds_4 = list(ds_3) # Dertivative wrt r of entropy

V_t4 = arrayLst( a22 * rr[t]**n + b22 / rr[t] for t in range(pts))
# V_t4 = arrayLst( V_t4m * R_m / rr[t] for t in range(pts))


rV_t4  = arrayLst(rr2[t] * V_t4[t] for t in range(pts))
drV_t4 = finDiff(rV_t4,deltaR2)

h_t4 = arrayLst(h_t3[t] + U[t] * (V_t4[t] - V_t3[t]) for t in range(pts)) # Total enthalpy in 2 
T_t4 = h_t4 / c_p # Total temperature

T_4 = T_4m * np.ones(pts) # Static temperature

dh_t4 = finDiff(h_t4,deltaR2)


# This loop can be avoided using flaired blades b_4 != b_3
while abs(err) > tol: # Begin loop to get mass flow convergence
    print("")
    print("---Iteration no. " + str(iter))

    V_a4 = list(np.zeros(pts)) # Create the list
    V_a4[mean_index] = V_a4m  # Initial value for forward integration starting from mean radius
    dV_a4 = list(np.zeros(pts))

    # N.I.S.R.E. 4 numerical integration 
    for j in list(range(0,mean_index)):
        for q,k in zip([mean_index + j, mean_index - j],[1,-1]):
            dV_a4[q] = 1 / V_a4[q] * ( dh_t4[q] - T_4[q] * ds_4[q] - V_t4[q] / rr2[q] * drV_t4[q] )
            V_a4[q + k*1] = V_a4[q] + dV_a4[q] * k * deltaR2 

    # Initiate all the lists
    V_4 , alpha_4, W_t4, W_a4, W_4, beta_4, p_4, rho_4, M_4, M_4r, p_t4, p_t4r, integrand_4, chi_2, L_eul = (list(np.zeros(pts)) for t in range(15))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        V_4[j] = np.sqrt(V_a4[j]**2 + V_t4[j]**2)
        alpha_4[j] = np.arctan(V_t4[j]/V_a4[j])
        W_t4[j] = V_t4[j] - U[j]
        W_a4[j] = V_a4[j]
        W_4[j] = np.sqrt(W_t4[j]**2 + W_a4[j]**2)
        beta_4[j] = np.arctan(W_t4[j]/W_a4[j])        

        # Thermodynamics
        T_4[j] = T_t4[j] - V_4[j]**2 / (2 * c_p)
        p_4[j] = p_4m * (T_4[j] / T_4m)**(gamma/(gamma-1)) * np.exp(- (s_4[j] - s_4[mean_index]) / R)
        rho_4[j] = p_4[j] / (R*T_4[j])
        M_4[j]  = V_4[j] / np.sqrt(gamma * R * T_4[j])
        M_4r[j] = W_4[j] / np.sqrt(gamma * R * T_4[j])
        p_t4[j] = p_4[j]*(1 + (gamma-1) / 2 * M_4[j]**2 ) ** (gamma/(gamma-1))


        L_eul[j] = U[j] * (V_t4[j] - V_t3[j])
        chi_2[j] = (W_3[j]**2 - W_4[j]**2) / (2 * L_eul[j])

        integrand_4[j] = 2 * np.pi * rr2[j] * rho_4[j] * V_a4[j] 
                

        p_t4r[j] = p_t3r[j] - omega_loss_R * (p_t3r[j] - p_3[j])

        # ENTROPY EVALUATION

        s_4[j]  = s_3[j] - R * np.log(p_t4r[j] / p_t3r[j])

    ds_4 = finDiff(s_4,deltaR2)

    m_dot_trap = np.trapz(integrand_4, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a4m = V_a4m*(1 + err) # New axial velocity
    
    
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a4m = " + str(V_a4m))
    print("err = "+ str(err))
    iter += 1














print("")
print("########## STATOR 2 OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-3 # Tolerance of error wrt the desires mass flow value
iter = 1

# Input data
omega_loss_S = 0.08

V_t5 = list( V_t5m * R_m / rr2[t]  for t in range(pts))

rV_t5  = arrayLst(rr2[t] * V_t5[t] for t in range(pts))
drV_t5 = finDiff(rV_t5,deltaR)

# Initial assumptions
T_5  = list(T_5m * np.ones(pts))
s_5  = s_4[:]
ds_5 = s_4[:]

# Imposed by thermodynamics
h_t5 = h_t4
dh_t5 = dh_t4
T_t5 = T_t4

# This loop can be avoided using flaired blades b_2 != b_1
while abs(err) > tol: # Begin loop to get mass flow convergence
    print("")
    print("---Iteration no. " + str(iter))

    V_a5 = list(np.zeros(pts)) # Create the list
    V_a5[mean_index] = V_a5m  # Initial value for forward integration starting from mean radius
    dV_a5 = list(np.zeros(pts))
    
    # N.I.S.R.E at stator outlet (5)
    for j in list(range(0,mean_index)):
        for q,k in zip([mean_index + j, mean_index - j],[1,-1]):
            dV_a5[q] = 1 / V_a5[q] * ( dh_t5[q] - T_5[q] * ds_5[q] - V_t5[q] / rr2[q] * drV_t5[q] )
            V_a5[q + k*1] = V_a5[q] + dV_a5[q] * k * deltaR 

    # Initiate all the lists
    V_5 , alpha_5, p_5, rho_5, M_5, p_t5, integrand_5, W_t5, W_a5, W_5, beta_5, M_5r, p_t5, p_t5r = (list(np.zeros(pts)) for t in range(14))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        alpha_5[j] = np.arctan(V_t5[j]/V_a5[j])
        V_5[j] = np.sqrt(float(V_a5[j]**2 + V_t5[j]**2))
        W_t5[j] = V_t5[j] - U[j]
        W_a5[j] = V_a5[j]
        W_5[j] = np.sqrt(W_t5[j]**2 + W_a5[j]**2)
        beta_5[j] = np.arctan(W_t5[j]/W_a5[j])
        
        # Thermodynamics
        T_5[j] = T_t5[j] - V_5[j]**2 / (2 * c_p)
        p_5[j] = p_5m * (T_5[j] / T_5m)**(gamma/(gamma-1)) * np.exp(- (s_5[j] - s_5[mean_index]) / R)
        rho_5[j] = p_5[j] / (R*T_5[j])
        M_5[j]   = V_5[j] / np.sqrt(gamma * R * T_5[j])
        M_5r[j]  = W_5[j] / np.sqrt(gamma * R * T_5[j])
        p_t5[j]  = p_5[j]*(1 + (gamma-1) / 2 * M_5[j]**2  ) ** (gamma/(gamma-1))
        p_t5r[j] = p_5[j]*(1 + (gamma-1) / 2 * M_5r[j]**2 ) ** (gamma/(gamma-1))

        integrand_5[j] = 2 * np.pi * rr2[j] * rho_5[j] * V_a5[j]
        
        #Evaluate the q.ties in section 1 (np.expressions) at the current radius
        # tmp = overwritten at every iteration, no need for a new array for _1 quantities
        
        p_t5[j] = p_t4[j] - omega_loss_S * (p_t4[j] - p_4[j])

        # ENTROPY EVALUATION

        s_5[j]  = s_4[j] - R * np.log(p_t5[j] / p_t4[j])

    ds_5 = finDiff(s_5,deltaR) # Derivative of s_5

    m_dot_trap = np.trapz(integrand_5, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a5m = V_a5m*(1 + err) # New axial velocity
    
    

    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a5m = "+ str(V_a5m) + " [m/s]")
    print("err = "+ str(err))
    iter += 1




















############### Blade design (Stage 1) ##############


# Design Parameters from computations above
# Inlet
beta_1root = beta_1[0]          * 180 / pi #deg
beta_1mid  = beta_1[mean_index] * 180 / pi #deg
beta_1tip  = beta_1[-1]         * 180 / pi #deg

# Outlet
beta_2root = beta_2[0]          * 180 / pi #deg
beta_2mid  = beta_2[mean_index] * 180 / pi #deg
beta_2tip  = beta_2[-1]         * 180 / pi #deg

# Deflection
deltabeta1_root = beta_1root - beta_2root
deltabeta1_mid  = beta_1mid  - beta_2mid 
deltabeta1_tip  = beta_1tip  - beta_2tip 

# Initial hypothesis ----> NOTE: TUNABLE
percent_th = 8                # [%] Max thickness WRT chord of blade profile 
chord      = 0.06             # [m] Starting point from reference procedure
th = chord * percent_th / 100 # [m] Actual thickness
solidity   = 1                # [ ] ! Initial assumption at midspan


s_mid = solidity * chord

n_blade = round( 2 * np.pi * R_m / s_mid) # Number of blades

# Pitch along blade span
s_mid  = 2 * np.pi * R_m / n_blade
s_tip  = 2 * np.pi * R_t / n_blade
s_root = 2 * np.pi * R_h / n_blade

# Solidity along blade span, recomputed after choosing n blades
sigma_mid  = chord / s_mid
sigma_tip  = chord / s_tip
sigma_root = chord / s_root

# Equivalent camber theta: from graphs on slide 9 ppt
# NOTE: TUNABLE
theta_eq_mid  = 13
theta_eq_tip  = 23
theta_eq_root = 39

# compute C_l = theta/25
cl_mid  = theta_eq_mid  / 25
cl_tip  = theta_eq_tip  / 25
cl_root = theta_eq_root / 25


# TODO: Insert semiempirical correlations below instead of graph values

# Zero camber incidence angle

# Semiempirical exponents
p_mid  = 0.914 + sigma_mid **3 / 160
p_tip  = 0.914 + sigma_tip **3 / 160
p_root = 0.914 + sigma_root**3 / 160

# Incidence angles in degrees
i0_10_mid  = abs(beta_1mid)**p_mid   / (5 + 46*np.exp(-2.3 * sigma_mid))  - 0.1 * sigma_mid**3  * np.exp((abs(beta_1mid)-70)/4)
i0_10_tip  = abs(beta_1tip)**p_tip   / (5 + 46*np.exp(-2.3 * sigma_tip))  - 0.1 * sigma_tip**3  * np.exp((abs(beta_1tip)-70)/4) 
i0_10_root = abs(beta_1root)**p_root / (5 + 46*np.exp(-2.3 * sigma_root)) - 0.1 * sigma_root**3 * np.exp((abs(beta_1root)-70)/4)


# Correction factor for profile thickness
q = 0.28 / (0.1 + (th/chord)**0.3) # Empirical exponent

k_ith_mid  = (10*th/chord)**q
k_ith_tip  = (10*th/chord)**q
k_ith_root = (10*th/chord)**q


# Slope factor n
n_mid  = 0.025 * sigma_mid  - 0.06 - ( (abs(beta_1mid)/90)  ** (1+1.2*sigma_mid)  ) / (1.5 + 0.43 * sigma_mid)
n_tip  = 0.025 * sigma_tip  - 0.06 - ( (abs(beta_1tip)/90)  ** (1+1.2*sigma_tip)  ) / (1.5 + 0.43 * sigma_tip)
n_root = 0.025 * sigma_root - 0.06 - ( (abs(beta_1root)/90) ** (1+1.2*sigma_root) ) / (1.5 + 0.43 * sigma_root)


# Optimal Incidence Angle
i_opt_mid  = i0_10_mid  * k_ith_mid  + n_mid  * theta_eq_mid
i_opt_tip  = i0_10_tip  * k_ith_tip  + n_tip  * theta_eq_tip
i_opt_root = i0_10_root * k_ith_root + n_root * theta_eq_root

# Deviation Angle
delta0_mid  = 0.01 * sigma_mid  * abs(beta_1mid)  + (0.74*sigma_mid**1.9 + 3 * sigma_mid) * (abs(beta_1mid) / 90) ** (1.67 + 1.09 * sigma_mid) 
delta0_tip  = 0.01 * sigma_tip  * abs(beta_1tip)  + (0.74*sigma_tip**1.9 + 3 * sigma_tip) * (abs(beta_1tip) / 90) ** (1.67 + 1.09 * sigma_tip)  
delta0_root = 0.01 * sigma_root * abs(beta_1root) + (0.74*sigma_root**1.9 + 3 * sigma_root) * (abs(beta_1root) / 90) ** (1.67 + 1.09 * sigma_root) 


# Correction for thickness effects on deviation
kdeltath_mid  = 6.25 * (th/chord) + 37.5 * (th/chord)**2    
kdeltath_tip  = 6.25 * (th/chord) + 37.5 * (th/chord)**2   
kdeltath_root = 6.25 * (th/chord) + 37.5 * (th/chord)**2  


# Exponent Factor b
b_mid  = 0.9625 - 0.17 * abs(beta_1mid)/100 - 0.85 * (abs(beta_1mid)/100)**3
b_tip  = 0.9625 - 0.17 * abs(beta_1tip)/100 - 0.85 * (abs(beta_1tip)/100)**3
b_root = 0.9625 - 0.17 * abs(beta_1root)/100 - 0.85 * (abs(beta_1root)/100)**3


# Slope factor m - m1 = Reference value for solidity = 1
m1_mid  = 0.17 - 0.0333 * abs(beta_1mid) /100 + 0.333 * (beta_1mid/100) **2
m1_tip  = 0.17 - 0.0333 * abs(beta_1tip) /100 + 0.333 * (beta_1tip/100) **2
m1_root = 0.17 - 0.0333 * abs(beta_1root)/100 + 0.333 * (beta_1root/100)**2


m_mid  = m1_mid  / (sigma_mid  **b_mid )
m_tip  = m1_tip  / (sigma_tip  **b_tip )
m_root = m1_root / (sigma_root **b_root)


#computation of deviation
delta_mid = delta0_mid * kdeltath_mid + m_mid * theta_eq_mid / (sigma_mid ** b_mid)
delta_tip = delta0_tip * kdeltath_tip + m_tip * theta_eq_tip / (sigma_tip ** b_tip)
delta_root = delta0_root * kdeltath_root + m_root * theta_eq_root / (sigma_root ** b_root)


#final delta beta
deltabetafinal_mid = theta_eq_mid - delta_mid + i_opt_mid
deltabetafinal_tip = theta_eq_tip - delta_tip + i_opt_tip
deltabetafinal_root = theta_eq_root - delta_root + i_opt_root

print(deltabetafinal_mid)
print(abs(deltabeta1_mid))

input()


















############################ PLOTS BELOW ############################


plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,W_1,"b")
plt.plot(rr,W_2,"g")
plt.plot(rr,W_3,"r")
plt.plot(rr,W_4,"c")
plt.ylabel(r" $W$ $[m/s]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Rotor 2 In", "Rotor 2 Out"])
plt.title("Relative Velocity")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,V_a1,"b")
plt.plot(rr,V_a2,"g")
plt.plot(rr,V_a3,"r")
plt.plot(rr,V_a4,"c")
plt.plot(rr,V_a5,"m")
plt.ylabel(r"$V_a$ $[m/s]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Axial Absolute Velocity")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,V_t1,"b")
plt.plot(rr,V_t2,"g")
plt.plot(rr,V_t3,"r")
plt.plot(rr,V_t4,"c")
plt.plot(rr,V_t5,"m")
plt.ylabel(r"$V_t$ $[m/s]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Tangential Absolute Velocity")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,p_1,"b")
plt.plot(rr,p_2,"g")
plt.plot(rr,p_3,"r")
plt.plot(rr,p_4,"c")
plt.plot(rr,p_5,"m")
plt.ylabel(r"$p$ $[Pa]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Static Pressure")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,p_t1,"b")
plt.plot(rr,p_t2,"g")
plt.plot(rr,p_t3,"r")
plt.plot(rr,p_t4,"c")
plt.plot(rr,p_t5,"m")
plt.ylabel(r"$p_t$ $[Pa]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Total Pressure")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,h_t1,"b")
plt.plot(rr,h_t2,"g")
plt.plot(rr,h_t3,"r")
plt.plot(rr,h_t4,"c")
plt.plot(rr,h_t5,"m")
plt.ylabel(r"$h_t$ $[Pa]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Total Enthalpy")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr, T_1 ,"b")
plt.plot(rr, T_2 ,"g")
plt.plot(rr, T_3 ,"r")
plt.plot(rr, T_4,"c")
plt.plot(rr, T_5,"m")
plt.ylabel(r"$T$ $[K]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Static Temperature")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,rho_1,"b")
plt.plot(rr,rho_2,"g")
plt.plot(rr,rho_3,"r")
plt.plot(rr,rho_4,"c")
plt.plot(rr,rho_5,"m")
plt.ylabel(r"$\rho$ $[kg/m^3]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Density")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,s_1,"b")
plt.plot(rr,s_2,"g")
plt.plot(rr,s_3,"r")
plt.plot(rr,s_4,"c")
plt.plot(rr,s_5,"m")
plt.ylabel(r"$s$ $[J/K]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Entropy")
plt.grid(alpha=0.2)

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,180/np.pi * np.array(alpha_1),"b")
plt.plot(rr,180/np.pi * np.array(alpha_2),"g")
plt.plot(rr,180/np.pi * np.array(alpha_3),"r")
plt.plot(rr,180/np.pi * np.array(alpha_4),"c")
plt.plot(rr,180/np.pi * np.array(alpha_5),"m")
plt.ylabel(r"$\alpha$ [deg]")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
plt.title("Absolute Flow Angle")
plt.grid(alpha=0.2)
 
 
plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,chi)
plt.plot(rr2,chi_2)
plt.ylabel(r"$\chi$")
plt.xlabel(r"$r \  [m]$")
plt.title("Reaction Degree")
plt.legend(["Stage 1","Stage 2"])
plt.grid(alpha=0.2)

# This should be constant if a free vortex distribution is used
# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,L_eul)
# plt.title("Euler Work [J/kg]")

# Plot inlet and outlet velocity triangles at hub, mean radius and tip
# P stands for plotting


print("Average Exit Total Pressure = " , np.average(p_t5))










fig, axs = plt.subplots(3,1, sharex = True,  figsize=(4, 7), dpi=65) # Create figure

j = 0 # Index used to move through the subplots
for i, name in zip([R_t, R_m, R_h], ["Tip", "Mean", "Hub"]):
    
    # Find the index of the radius we are considering
    index = np.where(np.isclose(rr, i))
    index = (index[0])[0]
    
    #Evaluate q.ties at that radius
    U_P   = U[index]
    V_a1P = V_a1[index]
    V_t1P = V_t1[index]
    W_a1P = W_a1[index]
    W_t1P = W_t1[index]
    V_a2P = V_a2[index]
    V_t2P = V_t2[index]
    W_a2P = W_a2[index]
    W_t2P = W_t2[index]

    # axs[j].grid(alpha=0.2) #Add grid
    
    #Plot inlet and outlet triangles
    axs[j].quiver([0,U_P - V_t1P, U_P - V_t1P] , [0,V_a1P,V_a1P] , [U_P,V_t1P,W_t1P] , [0,-V_a1P,-W_a1P] , angles='xy',scale_units='xy', scale=1.0, color=["black","blue","blue"])
    axs[j].quiver([0,U_P - V_t2P, U_P - V_t2P] , [0,V_a2P,V_a2P] , [U_P,V_t2P,W_t2P] , [0,-V_a2P,-W_a2P] , angles='xy',scale_units='xy', scale=1.,  color=["black","green","green"])
    
    axs.flat[j].set_xlim(-50, 20 + U[-1]) #Set the limits for the x axis
    axs.flat[j].set_ylim(-5,  20 + max(V_a2[0],V_a1[-1]) )  #Set the limits for the y axis
    
    axs[j].set_aspect('equal') #Equal aspect ratio axes
    axs[j].set_ylabel(r"Axial Component $[m/s]$")
    axs[j].set_title(name)

    j = j+1

axs[2].set_xlabel(r"Tangential Component $[m/s]$")

fig, axs = plt.subplots(3,1, sharex=True, sharey=True, figsize=(4, 7), dpi=65) # Create figure

j = 0 # Index used to move through the subplots
for i, name in zip([R_t, R_m, R_h], ["Tip", "Mean", "Hub"]):

    # Find the index of the radius we are considering
    index = np.where(np.isclose(rr, i))
    index = (index[0])[0]

    #Evaluate q.ties at that radius
    U_P   = U[index]
    V_a2P = V_a2[index]
    V_t2P = V_t2[index]
    V_a3P = V_a3[index]
    V_t3P = V_t3[index]

    # axs[j].grid(alpha=0.2) #Add grid
    
    #Plot inlet and outlet triangles
    axs[j].quiver([0,U_P - V_t2P] , [0,V_a2P] , [U_P,V_t2P] , [0,-V_a2P] , angles='xy',scale_units='xy', scale=1.0, color=["black","green"])
    axs[j].quiver([0,U_P - V_t3P] , [0,V_a3P] , [U_P,V_t3P] , [0,-V_a3P] , angles='xy',scale_units='xy', scale=1.,  color=["black","red"])
    
    
    axs.flat[j].set_xlim(-50, 20 + U[-1]) #Set the limits for the x axis
    axs.flat[j].set_ylim(-5, 20 + max(float(V_a2[0]), float(V_a2[-1])))  #Set the limits for the y axis
    
    axs[j].set_aspect('equal') #Equal aspect ratio axes
    axs[j].set_ylabel(r"Axial Component $[m/s]$")
    axs[j].set_title(name)

    j = j+1

axs[2].set_xlabel(r"Tangential Component $[m/s]$")


plt.show()
