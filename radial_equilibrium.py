### Radial equilibrium script ###

exec(open("./turboproject.py").read()) # Run mean line design

import numpy as np
import matplotlib.pyplot as plt

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


# Thermophysical properties
c_p = 1005  # Constant pressure specific heat [J/(kg K)]
gamma = 1.4 # Specific heat ratio
c_v = c_p/gamma
R = c_p * (gamma-1)/gamma # Gas constant [J/(kg K)]

# Discretization
pts = 101  # Total number of points across the radius, UNEVEN!
rr = np.linspace(R_h, R_t, pts) # Discrete space of the radii over which we compute our quantities
deltaR = (R_t - R_h)/ (pts - 1) # Radius interval between points
mean_index = pts//2  # Index of the various lists corresponding to mean radius quantities

omega = 2 * pi * rpm/60           # Angular velocity     [rad/s]
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
V_t1 = arrayLst(R_m * V_t1m / rr[t] for t in range(pts)) # e.g. Free vortex distribution r * V_t = const

rV_t1 = arrayLst(rr[t] * V_t1[t] for t in range(pts))

drV_t1 = finDiff(rV_t1, deltaR)

print("")
print("########## INLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 0.001
rel = 1
iter= 1

# This loop can be eliminated by varying b_1 to accomodate for m_dot_req
while abs(err) > tol: 
    
    V_a1 = np.zeros(pts) # Create the list
    V_a1[mean_index] = V_a1m  # Initial value for forward integration starting from mean radius
    dV_a1 = np.zeros(pts)

    # N.I.S.R.E. 1 numerical integration 
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
    V_a1m = V_a1m*(1 + err * rel) # New axial velocity
    
    
    print("")
    print("---Iteration no. " + str(iter))
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("err = "+ str(err))
    print("V_a1m = "+ str(V_a1m) + " [m/s]")
    iter += 1


print("")
print("########## OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 0.001 # Tolerance of error wrt the desires mass flow value
rel = 0.8 # Relaxation factor
iter = 1


# Inputs

# Entropy inputs, NOTE: absolute values are meaningless
omega_loss_R = 0.5 # Coefficient of loss

# Need to transform s_2 and ds_2 into lists otherwise numpy will assign the same id to s_1 and s_2, even with s_2 = s_1[:] why??
s_2  = list( s_1)    # Initial radial entropy distribution in 2
ds_2 = list(ds_1) # Dertivative wrt r of entropy

V_t2 = arrayLst(V_t2m * R_m / rr[t] for t in range(pts)) # Outlet tangential velocity distribution (e.g. free vortex)

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
    # --> Start from V_2m at R_m and move forwards and backwards up to R_t and R_h
    # j moves from 1 to the mean_index
    # q and k are a subloop to simplify the code, the first values of q,k corresponding to
    # the "forwards" integration, and the second values to the "backwards" integration
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
    V_a2m = V_a2m*(1 + err * rel) # New axial velocity
    
    
    print("")
    print("---Iteration no. " + str(iter))
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("err = "+ str(err))
    iter += 1

print("V_a2m = " + str(V_a2m))


# Plot inlet and outlet velocity triangles at hub, mean radius and tip
# P stands for plotting

fig, axs = plt.subplots(3,1, sharex=True, sharey=True, figsize=(4, 7), dpi=65) # Create figure

j = 0 # Index used to move through the subplots
for i in [R_t, R_m, R_h]:
    # Evaluate the quantities to plot on the desired radius
    
    index = np.where(np.isclose(rr, i))
    index = (index[0])[0]

    U_P   = U[index]
    V_a1P = V_a1[index]
    V_t1P = V_t1[index]
    W_a1P = W_a1[index]
    W_t1P = W_t1[index]
    V_a2P = V_a2[index]
    V_t2P = V_t2[index]
    W_a2P = W_a2[index]
    W_t2P = W_t2[index]

    # axs[j].grid() #Add grid
    
    #Plot inlet and outlet triangles
    axs[j].quiver([0,U_P - V_t1P, U_P - V_t1P] , [0,V_a1P,V_a1P] , [U_P,V_t1P,W_t1P] , [0,-V_a1P,-W_a1P] , angles='xy',scale_units='xy', scale=1.0, color=["black","blue","blue"])
    axs[j].quiver([0,U_P - V_t2P, U_P - V_t2P] , [0,V_a2P,V_a2P] , [U_P,V_t2P,W_t2P] , [0,-V_a2P,-W_a2P] , angles='xy',scale_units='xy', scale=1.,  color=["black","red","red"])
    
    axs.flat[j].set_xlim(-50, 20 + U[-1]) #Set the limits for the x axis
    axs.flat[j].set_ylim(-5,  20 + max(V_a2[0],V_a1[-1]) )  #Set the limits for the y axis
    
    axs[j].set_aspect('equal') #Equal aspect ratio axes

    j = j+1



print("")
print("########## STATOR ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 0.001 # Tolerance of error wrt the desires mass flow value
rel = 0.8 # Relaxation factor
iter = 1

# Input data
omega_loss_S = 0.5

V_t3 = list(R_m * V_t3m / radius for radius in rr)

drV_t3 = list(np.zeros(pts)) # Free vortex distribution

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

    V_a3 = list(np.zeros(pts)) # Create the list
    V_a3[mean_index] = V_a3m  # Initial value for forward integration starting from mean radius
    dV_a3 = list(np.zeros(pts))
    
    # N.I.S.R.E at stator outlet (3)
    for j in list(range(0,mean_index)):
        for q,k in zip([mean_index + j, mean_index - j],[1,-1]):
            dV_a3[q] = 1 / V_a3[q] * ( dh_t3[q] - T_3[q] * ds_3[q] - V_t3[q] / rr[q] * drV_t3[q] )
            V_a3[q + k*1] = V_a3[q] + dV_a3[q] * k * deltaR 

    # Initiate all the lists
    V_3 , alpha_3, p_3, rho_3, M_3, p_t3, integrand_3 = (list(np.zeros(pts)) for t in range(7))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        alpha_3[j] = np.arctan(V_t3[j]/V_a3[j])
        V_3[j] = np.sqrt(float(V_a3[j]**2 + V_t3[j]**2))

        # Thermodynamics
        T_3[j] = T_t3[j] - V_3[j]**2 / (2 * c_p)
        p_3[j] = p_3m * (T_3[j] / T_3m)**(gamma/(gamma-1)) * np.exp(- (s_3[j] - s_3[mean_index]) / R)
        rho_3[j] = p_3[j] / (R*T_3[j])
        M_3[j]  = V_3[j] / np.sqrt(gamma * R * T_3[j])

        integrand_3[j] = 2 * np.pi * rr[j] * rho_3[j] * V_a3[j]
        
        #Evaluate the q.ties in section 1 (np.expressions) at the current radius
        # tmp = overwritten at every iteration, no need for a new array for _1 quantities
        
        p_t3[j] = p_t2[j] - omega_loss_S * (p_t2[j] - p_2[j])

        # ENTROPY EVALUATION

        s_3[j]  = s_2[j] - R * np.log(p_t3[j] / p_t2[j])

    ds_3 = finDiff(s_3,deltaR) # Derivative of s_3

    m_dot_trap = np.trapz(integrand_3, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a3m = V_a3m*(1 + err * rel) # New axial velocity
    
    
    print("")
    print("---Iteration no. " + str(iter))
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a3m = "+ str(V_a3m) + " [m/s]")
    print("err = "+ str(err))
    iter += 1


fig, axs = plt.subplots(3,1, sharex=True, sharey=True, figsize=(4, 7), dpi=65) # Create figure

j = 0 # Index used to move through the subplots
for i in [R_t, R_m, R_h]:

    index = np.where(np.isclose(rr, i))
    index = (index[0])[0]

    U_P   = U[index]
    V_a2P = V_a2[index]
    V_t2P = V_t2[index]
    V_a3P = V_a3[index]
    V_t3P = V_t3[index]

    # axs[j].grid() #Add grid
    
    #Plot inlet and outlet triangles
    axs[j].quiver([0,U_P - V_t2P] , [0,V_a2P] , [U_P,V_t2P] , [0,-V_a2P] , angles='xy',scale_units='xy', scale=1.0, color=["black","blue"])
    axs[j].quiver([0,U_P - V_t3P] , [0,V_a3P] , [U_P,V_t3P] , [0,-V_a3P] , angles='xy',scale_units='xy', scale=1.,  color=["black","red"])
    
    axs.flat[j].set_xlim(-50, 20 + U[-1]) #Set the limits for the x axis
    axs.flat[j].set_ylim(-5, 20 + max(float(V_a2[0]), float(V_a2[-1])))  #Set the limits for the y axis
    
    axs[j].set_aspect('equal') #Equal aspect ratio axes

    j = j+1



# Some more plots
plt.figure(figsize=(5, 5), dpi=65)
plt.plot(rr,W_1)
plt.plot(rr,W_2)
plt.title("Relative Velocity [m/s]")
plt.legend(["1","2","3"])

plt.figure(figsize=(5, 5), dpi=65)
plt.plot(rr,V_1)
plt.plot(rr,V_2)
plt.plot(rr,V_3)
plt.title("Absolute Velocity [m/s]")
plt.legend(["1","2","3"])

plt.figure(figsize=(5, 5), dpi=65)
plt.plot(rr,p_1)
plt.plot(rr,p_2)
plt.plot(rr,p_3)
plt.title("Pressure [Pa]")
plt.legend(["1","2","3"])

plt.figure(figsize=(5, 5), dpi=65)
plt.plot(rr,chi)
plt.title("Degree of reaction")

# This should be constant if a free vortex distribution is used
# plt.figure(figsize=(5, 5), dpi=65)
# plt.plot(rr,L_eul)
# plt.title("Euler Work")

plt.show()
