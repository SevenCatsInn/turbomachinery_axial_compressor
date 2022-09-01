
### Radial equilibrium script ###

exec(open("./turboproject.py").read()) # Run mean line design

import numpy as np
import matplotlib.pyplot as plt
from blade_design import *
from losses import *
from off_design import *

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
s_0 = np.zeros(pts)
ds_0 = np.zeros(pts)
m_dot_req = 100 # Required mass flow [kg/s]
T_1 = T_1m * np.ones(pts)

# Computed Quantities
h_t1 = c_p * T_t1 # Total enthalpy [J/kg]
dh_t1 = finDiff(h_t1,deltaR)

# !! WARNING !!
#  Here the values pt1, Tt1 refer to the freestream, NOT ROTOR 1 INLET, I won't change them here yet to a 0 index (correct)
#  These values are nonetheless used to compute ds_0 (which is 0 actually lol) 
for j in range(pts): 
    ds_0[j]  = -R * finDiff(p_t1,deltaR)[j] / p_t1[j] + c_p * finDiff(T_t1,deltaR)[j] / T_t1[j] # Derivative over r of entropy [J/(kg K)]







print("")
print("########## DEFLECTOR INLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-3 # Tolerance of error wrt the desires mass flow value
iter = 1

# Input data

V_t0 = np.zeros(pts)

rV_t0  = arrayLst(rr[t] * V_t0[t] for t in range(pts))
drV_t0 = finDiff(rV_t0,deltaR)

# Initial assumptions
T_0  = list(T_0m * np.ones(pts))
s_0  = np.zeros(pts)
# ds0 TODO

# Imposed by thermodynamics
h_t0 = h_t1
dh_t0 = dh_t1
T_t0 = T_t1

# This loop can be avoided using flaired blades b_2 != b_1
while abs(err) > tol: # Begin loop to get mass flow convergence
    print("")
    print("---Iteration no. " + str(iter))

    V_a0 = list(np.zeros(pts)) # Create the list
    V_a0[mean_index] = V_a0m  # Initial value for forward integration starting from mean radius
    dV_a0 = list(np.zeros(pts))
    
    # N.I.S.R.E at stator outlet (0)
    for j in list(range(0,mean_index)):
        for q,k in zip([mean_index + j, mean_index - j],[1,-1]):
            dV_a0[q] = 1 / V_a0[q] * ( dh_t0[q] - T_0[q] * ds_0[q] - V_t0[q] / rr[q] * drV_t0[q] )
            V_a0[q + k*1] = V_a0[q] + dV_a0[q] * k * deltaR

    # Initiate all the lists
    V_0 , alpha_0, p_0, rho_0, M_0, p_t0, integrand_0, W_t0, W_a0, W_0, beta_0, M_0r, p_t0, p_t0r = (np.zeros(pts) for t in range(14))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        alpha_0[j] = np.arctan(V_t0[j]/V_a0[j])
        V_0[j] = np.sqrt(float(V_a0[j]**2 + V_t0[j]**2))
        W_t0[j] = V_t0[j] - U[j]
        W_a0[j] = V_a0[j]
        W_0[j] = np.sqrt(W_t0[j]**2 + W_a0[j]**2)
        beta_0[j] = np.arctan(W_t0[j]/W_a0[j])
        
        # Thermodynamics
        T_0[j] = T_t0[j] - V_0[j]**2 / (2 * c_p)
        p_0[j] = p_0m * (T_0[j] / T_0m)**(gamma/(gamma-1)) * np.exp(- (s_0[j] - s_0[mean_index]) / R)
        rho_0[j] = p_0[j] / (R*T_0[j])
        M_0[j]   = V_0[j] / np.sqrt(gamma * R * T_0[j])
        M_0r[j]  = W_0[j] / np.sqrt(gamma * R * T_0[j])
        p_t0[j]  = p_0[j]*(1 + (gamma-1) / 2 * M_0[j]**2  ) ** (gamma/(gamma-1))
        p_t0r[j] = p_0[j]*(1 + (gamma-1) / 2 * M_0r[j]**2 ) ** (gamma/(gamma-1))

        integrand_0[j] = 2 * np.pi * rr[j] * rho_0[j] * V_a0[j]
        
        #Evaluate the q.ties in section 1 (np.expressions) at the current radius
        # tmp = overwritten at every iteration, no need for a new array for _1 quantities
        

    m_dot_trap = np.trapz(integrand_0, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a0m = V_a0m*(1 + err) # New axial velocity
    
    

    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a0m = "+ str(V_a0m) + " [m/s]")
    print("err = "+ str(err))
    iter += 1













# Set the design choice for tangential velocity distribution in the radial direction

# First power vortex distribution
# V_t1 = arrayLst( V_t1m * R_m / rr[t] for t in range(pts)) # Free Vortex
V_t1 = arrayLst( a * rr[t]**n - b / rr[t] for t in range(pts)) # Power Design

rV_t1 = arrayLst(rr[t] * V_t1[t] for t in range(pts))
drV_t1 = finDiff(rV_t1, deltaR)

s_1  = arrayLst( s_0)    # Initial radial entropy distribution in 2
ds_1 = arrayLst(ds_0) # Dertivative wrt r of entropy

omega_loss_D = 0.02

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
    V_1 , alpha_1, W_t1, W_a1, W_1, beta_1, p_1, rho_1, M_1, M_1r, p_t1, p_t1r, integrand_1, chi, L_eul = (np.zeros(pts) for t in range(15))

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
        
        #Losses across the deflector
        chordD1=0.1
        solidityD1=1.0

        tmp_staggerD1 = 25 * np.array([0.4146902837280073, 8.562112644741685, 15.346296015211571]) # Times 25 is just because we are reusing a function for computing the C_l distribution that divides by 25
        tmp2_staggerD1 = compute_C_l(tmp_staggerD1,pts)
        staggerD1 = np.average(tmp2_staggerD1)

        NrowD1=1
        bladesD1=20
        shrouded_D1=0
        statorD1=1 #flag to say if the stage is a stator or a rotor: necessary in losses function for the profile losses, check there

        omega_overall_D1= losses(rr[j],chordD1,R_m,b_1,V_a1[j],V_a0[j],beta_0[j],beta_1[j],alpha_0[j],alpha_1[j],V_0[j],V_1[j],W_a0[j],W_a1[j],W_0[j],W_1[j],rho_0[j],rho_1[j],staggerD1,NrowD1,bladesD1,mdot,p_t0[j],p_0[j],shrouded_D1,statorD1) # Coefficient of loss
        p_t1[j] = p_t0[j] - omega_overall_D1 * (p_t0[j] - p_0[j])
        
        integrand_1[j] = 2 * np.pi * rr[j] * rho_1[j] * V_a1[j] 

        # ENTROPY EVALUATION

        s_1[j]  = s_0[j] - R * np.log(p_t1[j] / p_t0[j])
    
    ds_1 = finDiff(s_1,deltaR)
    
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
omega_loss_R = 0.02 # Coefficient of loss

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
    V_2 , alpha_2, W_t2, W_a2, W_2, beta_2, p_2, rho_2, M_2, M_2r, p_t2, p_t2r, integrand_2, chi, L_eul = (np.zeros(pts) for t in range(15))

    Beta = np.zeros(pts)
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
        
        Beta[j]=(T_t2[j]/T_t1[j])**(gamma/(gamma-1))

        L_eul[j] = U[j] * (V_t2[j] - V_t1[j])
        chi[j] = (W_1[j]**2 - W_2[j]**2) / (2 * L_eul[j])

        integrand_2[j] = 2 * np.pi * rr[j] * rho_2[j] * V_a2[j] 

        #LOSSES across first rotor: between section 1 and 2        
        chordR1=0.1
        
        tmp_staggerR1 = 25 * np.array([-11.676659681774582, -24.956847050581516, -47.56992271320476]) # Times 25 is just because we are reusing a function for computing the C_l distribution that divides by 25
        tmp2_staggerR1 = compute_C_l(tmp_staggerR1,pts)
        staggerR1 = np.average(tmp2_staggerR1)
        
        NrowR1 = 2
        bladesR1 = 26
        shrouded_R1 = 0
        statorR1 = 1 #flag to say if the stage is a stator or a rotor: necessary in losses function for the profile losses, check there
        
        omega_overall_R1= losses(rr[j],chordR1,R_m,b_1,V_a2[j],V_a1[j],beta_1[j],beta_2[j],alpha_1[j],alpha_2[j],V_1[j],V_2[j],W_a1[j],W_a2[j],W_1[j],W_2[j],rho_1[j],rho_2[j],staggerR1,NrowR1,bladesR1,mdot,p_t1r[j],p_1[j],shrouded_R1,statorR1) # Coefficient of loss # Coefficient of loss
        p_t2r[j] = p_t1r[j] - omega_overall_R1 * (p_t1r[j] - p_1[j])
        
        
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

# Geometry
R_h2 = R_m - b_2 / 2   # Hub Radius          [m]   
R_t2 = R_m + b_2 / 2  # Tip Radius          [m]

rr2 = np.linspace(R_h2, R_t2, pts) # Discrete space of the radii over which we compute our quantities
deltaR2 = (R_t2 - R_h2)/ (pts - 1) # Radius interval between points
mean_index = pts//2  # Index of the various lists corresponding to mean radius quantities
U2 = arrayLst( omega * rr2[t] for t in range(pts)) # Peripheral velocity  [m/s]

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5 # Tolerance of error wrt the desires mass flow value
iter = 1

# Input data
omega_loss_S = 0.02

V_t3 = arrayLst( a22 * rr2[t]**n - b22 / rr2[t] for t in range(pts))

# V_t3 = arrayLst( V_t3m * R_m / rr[t] for t in range(pts))

rV_t3  = arrayLst(rr2[t] * V_t3[t] for t in range(pts))
drV_t3 = finDiff(rV_t3,deltaR2)

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
            dV_a3[q] = 1 / V_a3[q] * ( dh_t3[q] - T_3[q] * ds_3[q] - V_t3[q] / rr2[q] * drV_t3[q] )
            V_a3[q + k*1] = V_a3[q] + dV_a3[q] * k * deltaR2 

    # Initiate all the lists
    V_3 , alpha_3, p_3, rho_3, M_3, p_t3, integrand_3, W_t3, W_a3, W_3, beta_3, M_3r, p_t3, p_t3r = (np.zeros(pts) for t in range(14))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        alpha_3[j] = np.arctan(V_t3[j]/V_a3[j])
        V_3[j] = np.sqrt(float(V_a3[j]**2 + V_t3[j]**2))
        W_t3[j] = V_t3[j] - U2[j]
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

        integrand_3[j] = 2 * np.pi * rr2[j] * rho_3[j] * V_a3[j]
        
        #Evaluate the q.ties in section 1 (np.expressions) at the current radius
        # tmp = overwritten at every iteration, no need for a new array for _1 quantities
        
        #LOSSES across stator of the first stage: between section 2 and 3        
        chordS1=0.1
        
        tmp_staggerS1 = 25 * np.array([17.465834987449313, 29.176544391755733, 46.28802243060806]) # Times 25 is just because we are reusing a function for computing the C_l distribution that divides by 25
        tmp2_staggerS1 = compute_C_l(tmp_staggerS1,pts)
        staggerS1= np.average(tmp2_staggerS1)

        NrowS1=3
        bladesS1=30
        shrouded_S1=0
        statorS1=1 #flag to say if the stage is a stator or a rotor: necessary in losses function for the profile losses, check there
        
        omega_overall_S1= losses(rr[j],chordS1,R_m,b_1,V_a3[j],V_a2[j],beta_2[j],beta_3[j],alpha_2[j],alpha_3[j],V_2[j],V_3[j],W_a2[j],W_a3[j],W_2[j],W_3[j],rho_2[j],rho_3[j],staggerS1,NrowS1,bladesS1,mdot,p_t2[j],p_2[j],shrouded_S1,statorS1) # Coefficient of loss # Coefficient of loss
    
        
        p_t3[j] = p_t2[j] - omega_overall_S1 * (p_t2[j] - p_2[j])

        # ENTROPY EVALUATION

        s_3[j]  = s_2[j] - R * np.log(p_t3[j] / p_t2[j])

    ds_3 = finDiff(s_3,deltaR2) # Derivative of s_3

    m_dot_trap = np.trapz(integrand_3, rr2)

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


err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5 # Tolerance of error wrt the desires mass flow value
iter = 1


# Inputs

# Entropy inputs, NOTE: absolute values are meaningless
omega_loss_R = 0.02 # Coefficient of loss

# Need to transform s_4 and ds_4 into lists otherwise numpy will assign the same id to s_3 and s_4, even with s_4 = s_3[:] why??
s_4  = list( s_3)    # Initial radial entropy distribution in 2
ds_4 = list(ds_3) # Dertivative wrt r of entropy

V_t4 = arrayLst( a22 * rr2[t]**n + b22 / rr2[t] for t in range(pts))
# V_t4 = arrayLst( V_t4m * R_m / rr[t] for t in range(pts))


rV_t4  = arrayLst(rr2[t] * V_t4[t] for t in range(pts))
drV_t4 = finDiff(rV_t4,deltaR2)

h_t4 = arrayLst(h_t3[t] + U2[t]* (V_t4[t] - V_t3[t]) for t in range(pts)) # Total enthalpy in 2 
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
    V_4 , alpha_4, W_t4, W_a4, W_4, beta_4, p_4, rho_4, M_4, M_4r, p_t4, p_t4r, integrand_4, chi_2, L_eul = (np.zeros(pts) for t in range(15))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        V_4[j] = np.sqrt(V_a4[j]**2 + V_t4[j]**2)
        alpha_4[j] = np.arctan(V_t4[j]/V_a4[j])
        W_t4[j] = V_t4[j] - U2[j]
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


        L_eul[j] = U2[j] * (V_t4[j] - V_t3[j])
        chi_2[j] = (W_3[j]**2 - W_4[j]**2) / (2 * L_eul[j])

        integrand_4[j] = 2 * np.pi * rr2[j] * rho_4[j] * V_a4[j] 
                
        #LOSSES across rotor of the second stage: between section 3 and 4        
        chordR2=0.1
        
        tmp_staggerR2 = 25 * np.array([-3.731720378208422, -20.285718712071258, -49.148251438192744]) # Times 25 is just because we are reusing a function for computing the C_l distribution that divides by 25
        tmp2_staggerR2 = compute_C_l(tmp_staggerR2,pts)
        staggerR2 = np.average(tmp2_staggerR2)

        NrowR2=4
        bladesR2=26
        shrouded_R2=0
        b_R2=0.3
        statorR2=0 #flag to say if the stage is a stator or a rotor: necessary in losses function for the profile losses, check there
        omega_overall_R2= losses(rr2[j],chordR2,R_m,b_2,V_a4[j],V_a3[j],beta_3[j],beta_4[j],alpha_3[j],alpha_4[j],V_3[j],V_4[j],W_a3[j],W_a4[j],W_3[j],W_4[j],rho_3[j],rho_4[j],staggerR2,NrowR2,bladesR2,mdot,p_t3r[j],p_3[j],shrouded_R2,statorR2) # Coefficient of loss 
        p_t4r[j] = p_t3r[j] - omega_overall_R2 * (p_t3r[j] - p_3[j])

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
omega_loss_S = 0.02

V_t5 = list( V_t5m / R_m * rr2[t]  for t in range(pts))

rV_t5  = arrayLst(rr2[t] * V_t5[t] for t in range(pts))
drV_t5 = finDiff(rV_t5,deltaR2)

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
            V_a5[q + k*1] = V_a5[q] + dV_a5[q] * k * deltaR2 

    # Initiate all the lists
    V_5 , alpha_5, p_5, rho_5, M_5, p_t5, integrand_5, W_t5, W_a5, W_5, beta_5, M_5r, p_t5, p_t5r = (np.zeros(pts) for t in range(14))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        alpha_5[j] = np.arctan(V_t5[j]/V_a5[j])
        V_5[j] = np.sqrt(float(V_a5[j]**2 + V_t5[j]**2))
        W_t5[j] = V_t5[j] - U2[j]
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
        #LOSSES across stator of the second stage: between section 4 and 5        
        chordS2=0.1
        
        tmp_staggerS2 = 25 * np.array([24.160523251671574, 32.05712108570312, 46.50612563659101]) # Times 25 is just because we are reusing a function for computing the C_l distribution that divides by 25
        tmp2_staggerS2 = compute_C_l(tmp_staggerS2,pts)
        staggerS2 = np.average(tmp2_staggerS2)

        NrowS2=4
        bladesS2=36
        shrouded_S2=0
        statorS2=1 #flag to say if the stage is a stator or a rotor: necessary in losses function for the profile losses, check there
        omega_overall_S2= losses(rr2[j],chordS2,R_m,b_2,V_a5[j],V_a4[j],beta_4[j],beta_5[j],alpha_4[j],alpha_5[j],V_4[j],V_5[j],W_a4[j],W_a5[j],W_4[j],W_5[j],rho_4[j],rho_5[j],staggerS2,NrowS2,bladesS2,mdot,p_t4[j],p_4[j],shrouded_S2,statorS2) # Coefficient of loss # Coefficient of loss        
        #Evaluate the q.ties in section 1 (np.expressions) at the current radius
        # tmp = overwritten at every iteration, no need for a new array for _1 quantities
        
        p_t5[j] = p_t4[j] - omega_loss_S * (p_t4[j] - p_4[j])

        # ENTROPY EVALUATION

        s_5[j]  = s_4[j] - R * np.log(p_t5[j] / p_t4[j])

    ds_5 = finDiff(s_5,deltaR2) # Derivative of s_5

    m_dot_trap = np.trapz(integrand_5, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a5m = V_a5m*(1 + err) # New axial velocity
    
    

    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a5m = "+ str(V_a5m) + " [m/s]")
    print("err = "+ str(err))
    iter += 1









print("")
print("--------------- BLADE DESIGN ---------------" )


############### Blade design Deflector ##############
percent_th0 = 10               # [%] Max thickness WRT chord of blade profile 
chord0      = 0.1             # [m] Starting point from reference procedure
solidity0   = 1.0              # [ ] ! Initial assumption at midspan
theta0 = [-0.8, -16, -28.3]

alpha_0 = np.zeros(pts)

print("")
print("###### DEFLECTOR BLADE DESIGN ######")


inc0, theta0, dev0, deltaAlpha0 = lieblein_design(alpha_0,alpha_1,percent_th0,chord0,solidity0, theta0, rr)

geom, prof_names = printPlot_blade(alpha_0,alpha_1, deltaAlpha0, inc0, theta0, percent_th0, chord0,pts)
plt.title("Deflector")








############### Blade design (Stage 1 Rotor) ##############

# Mechanical Properties
rho_b = 7850                   # [kg/m^3] Blade material density
stress_Y = 472e6               # [Pa]   Material Yield stress (Annealed 4340 steel)

percent_th1 = 10               # [%] Max thickness WRT chord of blade profile 
chord1      = 0.1              # [m] Starting point from reference procedure
solidity1   = 1.3              # [ ] ! Initial assumption at midspan
theta1 = [26.7, 15.8, -3.4]

print("")
print("###### STAGE 1 ROTOR BLADE DESIGN ######")

inc1, theta1, dev1, deltaBeta1 = lieblein_design(beta_1,beta_2,percent_th1,chord1,solidity1, theta1, rr)

geom1, prof_names1 = printPlot_blade(beta_1,beta_2, deltaBeta1, inc1, theta1, percent_th1, chord1, pts)
plt.title("Rotor Stage 1")


# Stress Calculations on First Rotor
stress_Ax1 = np.zeros(pts)

for j in range(pts):
    stress_Ax1[j] = rho_b*omega**2 * (R_t**2 - rr[j]**2) / 2


C_l1 = compute_C_l(theta1,pts)

integrand_tmp1 = np.zeros(pts)
for j in range(pts):
    integrand_tmp1[j] = ( 0.5 * (p_1[j]/(R*T_1[j])) * W_1[j]**2 * chord1 * C_l1[j] * (rr[j] - R_h)  ) # [N * m / m]


M_f_hub1 = np.trapz(integrand_tmp1, rr)

I_x1 = geom1[0][2] # First index is the Hub (1 mid, 2 tip) and second index is the geometrical information to extract

stress_Bend1 = M_f_hub1 / I_x1 * (percent_th1*chord1/100)/2 

stress_tot1 = abs(stress_Ax1[0]) + abs(stress_Bend1)

print("")
print("Total Stress at Hub = ", stress_tot1, "[Pa]")
print("Safety Factor =", stress_Y/stress_tot1)




############### Blade design (Stage 1 Stator) ##############

percent_th2 = 10               # [%] Max thickness WRT chord of blade profile 
chord2      = 0.1             # [m] Starting point from reference procedure
solidity2   = 1.5              # [ ] ! Initial assumption at midspan
theta2 = [9, 7, 35]

print("")
print("###### STAGE 1 STATOR BLADE DESIGN ######")

inc2, theta2, dev2, deltaAlpha2 = lieblein_design(alpha_2,alpha_3,percent_th2,chord2,solidity2, theta2, (rr+rr2)/2)

geom, prof_names = printPlot_blade(alpha_2,alpha_3, deltaAlpha2, inc2, theta2, percent_th2, chord2, pts)
plt.title("Stator Stage 1 ")






############### Blade design (Stage 2 Rotor) ##############

percent_th3 = 10               # [%] Max thickness WRT chord of blade profile 
chord3      = 0.1              # [m] Starting point from reference procedure
solidity3   = 1.3              # [ ] ! Initial assumption at midspan
theta3 = [27.5, 17.3, -7.5]

print("")
print("###### STAGE 2 ROTOR BLADE DESIGN ######")

inc3, theta3, dev3, deltaBeta3 = lieblein_design(beta_3,beta_4,percent_th3,chord3,solidity3, theta3, rr2)

geom3, prof_names = printPlot_blade(beta_3,beta_4, deltaBeta3, inc3, theta3, percent_th3, chord3, pts)
plt.title("Rotor Stage 2")



# Stress Calculations on Second Rotor
stress_Ax3 = np.zeros(pts)

for j in range(pts):
    stress_Ax3[j] = rho_b*omega**2 * (R_t2**2 - rr2[j]**2) / 2


C_l3 = compute_C_l(theta3,pts)

integrand_tmp3 = np.zeros(pts)
for j in range(pts):
    integrand_tmp3[j] = ( 0.5 * (p_3[j]/(R*T_3[j])) * W_3[j]**2 * chord3 * C_l3[j] * (rr2[j] - R_h2)  ) # [N * m / m]


M_f_hub3 = np.trapz(integrand_tmp3, rr2)

I_x3 = geom3[0][2] # First index is the Hub (1 mid, 2 tip) and second index is the geometrical information to extract

stress_Bend3 = M_f_hub3 / I_x3 * (percent_th3*chord3/100)/2 

stress_tot3 = abs(stress_Ax3[0]) + abs(stress_Bend3)

print("")
print("Total Stress at Hub = ", stress_tot3, "[Pa]")
print("Safety Factor =", stress_Y/stress_tot3)










############### Blade design (Stage 2 Stator) ##############

percent_th4 = 10               # [%] Max thickness WRT chord of blade profile 
chord4      = 0.1             # [m] Starting point from reference procedure
solidity4   = 1.8              # [ ] ! Initial assumption at midspan
theta4 = [14.5, 19.5, 48 ]

print("")
print("###### STAGE 2 STATOR BLADE DESIGN ######")

inc4, theta4, dev4, deltaAlpha4 = lieblein_design(alpha_4,alpha_5,percent_th4,chord4,solidity4, theta4, rr2)

geom, prof_names = printPlot_blade(alpha_4, alpha_5, deltaAlpha4, inc2, theta4, percent_th4, chord4, pts)
plt.title("Stator Stage 2")



# plt.show()

print("")
print("--------------- OFF-DESIGN ---------------")

mdot_off = 90
Leul1_off, beta1_off= off_design(R_m,mdot_off,beta2,rho2,Um,alpha1,rho1,gamma,efficiency_TT,cp,T_t1[mean_index],b1)
#Leul1_off, beta1_off= off_design(R_m,mdot_off,beta_2[mean_index],rho2[mean_index],U[mean_index],alpha_1[mean_index],rho1_[mean_index],gamma,efficiency_TT,cp,T_t1[mean_index],b1)

Leul2_off, beta2_off = off_design(R_m,mdot_off,beta4,rho4,Um,alpha3,rho3,gamma,efficiency_TT,cp,T_t3[mean_index],b2)
#Leul2_off, beta2_off= off_design(rr[mean_index],mdot_off,beta_4[mean_index],rho_4[mean_index],U[mean_index],alpha_3[mean_index],rho_3[mean_index],gamma,efficiency_TT,cp,T_t3[mean_index],b2)

print("Leul1_off=", Leul1_off)
print("beta1_off=", beta1_off)

print("Leul2_off=", Leul2_off)
print("beta2_off=", beta2_off)

beta_off_overall= beta1_off*beta2_off

print("beta_off_overall=",beta_off_overall)















print("")


############################ PLOTS BELOW ############################


# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,W_1,"b")
# plt.plot(rr,W_2,"g")
# plt.plot(rr2,W_3,"r")
# plt.plot(rr2,W_4,"c")
# plt.ylabel(r" $W$ $[m/s]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Rotor 2 In", "Rotor 2 Out"])
# plt.title("Relative Velocity")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,V_a1,"b")
# plt.plot(rr,V_a2,"g")
# plt.plot(rr2,V_a3,"r")
# plt.plot(rr2,V_a4,"c")
# plt.plot(rr2,V_a5,"m")
# plt.ylabel(r"$V_a$ $[m/s]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Axial Absolute Velocity")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,V_t1,"b")
# plt.plot(rr,V_t2,"g")
# plt.plot(rr2,V_t3,"r")
# plt.plot(rr2,V_t4,"c")
# plt.plot(rr2,V_t5,"m")
# plt.ylabel(r"$V_t$ $[m/s]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Tangential Absolute Velocity")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,p_0,"k")
# plt.plot(rr,p_1,"b")
# plt.plot(rr,p_2,"g")
# plt.plot(rr2,p_3,"r")
# plt.plot(rr2,p_4,"c")
# plt.plot(rr2,p_5,"m")
# plt.ylabel(r"$p$ $[Pa]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Deflector In","Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Static Pressure")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,p_t1,"b")
# plt.plot(rr,p_t2,"g")
# plt.plot(rr2,p_t3,"r")
# plt.plot(rr2,p_t4,"c")
# plt.plot(rr2,p_t5,"m")
# plt.ylabel(r"$p_t$ $[Pa]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Total Pressure")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,h_t1,"b")
# plt.plot(rr,h_t2,"g")
# plt.plot(rr2,h_t3,"r")
# plt.plot(rr2,h_t4,"c")
# plt.plot(rr2,h_t5,"m")
# plt.ylabel(r"$h_t$ $[Pa]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Total Enthalpy")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr, T_1 ,"b")
# plt.plot(rr, T_2 ,"g")
# plt.plot(rr2, T_3 ,"r")
# plt.plot(rr2, T_4,"c")
# plt.plot(rr2, T_5,"m")
# plt.ylabel(r"$T$ $[K]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Static Temperature")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,rho_1,"b")
# plt.plot(rr,rho_2,"g")
# plt.plot(rr2,rho_3,"r")
# plt.plot(rr2,rho_4,"c")
# plt.plot(rr2,rho_5,"m")
# plt.ylabel(r"$\rho$ $[kg/m^3]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Density")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,s_0,"k")
# plt.plot(rr,s_1,"b")
# plt.plot(rr,s_2,"g")
# plt.plot(rr2,s_3,"r")
# plt.plot(rr2,s_4,"c")
# plt.plot(rr2,s_5,"m")
# plt.ylabel(r"$s$ $[J/K]$")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Deflector In","Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Entropy")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,180/np.pi * np.array(alpha_0),"k")
# plt.plot(rr,180/np.pi * np.array(alpha_1),"b")
# plt.plot(rr,180/np.pi * np.array(alpha_2),"g")
# plt.plot(rr2,180/np.pi * np.array(alpha_3),"r")
# plt.plot(rr2,180/np.pi * np.array(alpha_4),"c")
# plt.plot(rr2,180/np.pi * np.array(alpha_5),"m")
# plt.ylabel(r"$\alpha$ [deg]")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Deflector In","Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Absolute Flow Angle")
# plt.grid(alpha=0.2)

# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,180/np.pi * np.array(beta_1),"b")
# plt.plot(rr,180/np.pi * np.array(beta_2),"g")
# plt.plot(rr2,180/np.pi * np.array(beta_3),"r")
# plt.plot(rr2,180/np.pi * np.array(beta_4),"c")
# plt.plot(rr2,180/np.pi * np.array(beta_5),"m")
# plt.ylabel(r"$\beta$ [deg]")
# plt.xlabel(r"$r \  [m]$")
# plt.legend(["Rotor In","Rotor Out","Stator Out","Rotor 2 Out", "Stator 2 Out"])
# plt.title("Relative Flow Angle")
# plt.grid(alpha=0.2)
 
 
# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,chi)
# plt.plot(rr2,chi_2)
# plt.ylabel(r"$\chi$")
# plt.xlabel(r"$r \  [m]$")
# plt.title("Reaction Degree")
# plt.legend(["Stage 1","Stage 2"])
# plt.grid(alpha=0.2)

# This should be constant if a free vortex distribution is used
# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,L_eul)
# plt.title("Euler Work [J/kg]")

# Plot inlet and outlet velocity triangles at hub, mean radius and tip
# P stands for plotting


print("Average Exit Total Pressure S1 = " , np.average(p_t3))
print("Average Exit Total Pressure S2 = " , np.average(p_t5))



# plt.show()






# VELOCITY TRIANGLES ACROSS ROTOR 1
fig, axs = plt.subplots(1,3, sharey = False,  figsize=(13, 5), dpi=70) # Create figure


j = 0 # Index used to move through the subplots
for i, name in zip([R_h, R_m, R_t], ["Hub", "Mean", "Tip"]):
    
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
    axs[j].quiver([0,U_P - V_t2P, U_P - V_t2P] , [0,V_a2P,V_a2P] , [U_P,V_t2P,W_t2P] , [0,-V_a2P,-W_a2P] , angles='xy',scale_units='xy', scale=1.0,  color=["black","green","green"])
    
    axs.flat[j].set_xlim(-50, 20 + U[-1]) #Set the limits for the x axis
    axs.flat[j].set_ylim(-5,  20 + max(V_a2[0],V_a1[-1]) )  #Set the limits for the y axis
    
    axs[j].set_aspect('equal') #Equal aspect ratio axes
    axs[j].set_ylabel(r"Axial Component $[m/s]$")
    axs[j].set_xlabel(r"Tangential Component $[m/s]$")
    axs[j].set_title(name)

    j = j+1




# VELOCITY TRIANGLES ACROSS STATOR 1

fig, axs = plt.subplots(1,3, sharey = False,  figsize=(13, 5), dpi=70) # Create figure

j = 0 # Index used to move through the subplots
for i, name in zip([R_h, R_m, R_t], ["Hub", "Mean", "Tip"]):

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
    axs[j].quiver([0,U_P - V_t3P] , [0,V_a3P] , [U_P,V_t3P] , [0,-V_a3P] , angles='xy',scale_units='xy', scale=1.0,  color=["black","red"])
    
    
    axs.flat[j].set_xlim(-50, 20 + U[-1]) #Set the limits for the x axis
    axs.flat[j].set_ylim(-5, 20 + max(float(V_a2[0]), float(V_a2[-1])))  #Set the limits for the y axis
    
    axs[j].set_aspect('equal') #Equal aspect ratio axes
    axs[j].set_ylabel(r"Axial Component $[m/s]$")
    axs[j].set_xlabel(r"Tangential Component $[m/s]$")
    axs[j].set_title(name)

    j = j+1











# VELOCITY TRIANGLES ACROSS ROTOR 2











# VELOCITY TRIANGLES ACROSS STATOR 2


plt.show()