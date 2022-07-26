### Radial equilibrium script ###

exec(open("./turboproject_S2.py").read()) # Run mean line design

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
pts = 200  # Total number of points across the radius, 
if pts % 2 == 0: pts = pts + 1 # Make pts uneven if it's even

rr = np.linspace(R_h, R_t, pts) # Discrete space of the radii over which we compute our quantities
deltaR = (R_t - R_h)/ (pts - 1) # Radius interval between points
mean_index = pts//2  # Index of the various lists corresponding to mean radius quantities

omega = 2 * pi * rpm/60                          # Angular velocity     [rad/s]
U = arrayLst( omega * rr[t] for t in range(pts)) # Peripheral velocity  [m/s]

# Input data
T_t1 = T_t3 # [K]     --> f(r)
p_t1 = p_t3 # [Pa] --> f(r)
s_1 = s_3
ds_1 = ds_3
m_dot_req = 100 # Required mass flow [kg/s]
T_1 = T_3

# Computed Quantities
h_t1 = h_t3 # Total enthalpy [J/kg]
dh_t1 = dh_t3

# Tangential Velocity
V_t1 = V_t3

rV_t1 = rV_t3
drV_t1 = drV_t3

print("")
print("########## ROTOR INLET ##########")

V_a1 = V_a3
V_1 = V_3
alpha_1 = alpha_3
W_t1 = W_t3 
W_a1 = V_a3
beta_1 = beta_3

# Thermodynamics
T_1 = T_3 
p_1 = p_3 
rho_1 = rho_3 
M_1  = M_3  
M_1r = M_3r 
p_t1  = p_t3  
p_t1r = p_t3r 

print("")
print("########## ROTOR OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5 # Tolerance of error wrt the desires mass flow value
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


# Plot inlet and outlet velocity triangles at hub, mean radius and tip
# P stands for plotting

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

    # axs[j].grid() #Add grid
    
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


print("")
print("########## STATOR OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 1e-5 # Tolerance of error wrt the desires mass flow value
iter = 1

# Input data
omega_loss_S = 0.5

V_t3 = list(R_m * V_t3m / rr[t] for t in range(pts))

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
    V_a3m = V_a3m*(1 + err) # New axial velocity
    
    
    print("")
    print("---Iteration no. " + str(iter))
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("V_a3m = "+ str(V_a3m) + " [m/s]")
    print("err = "+ str(err))
    iter += 1


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

    # axs[j].grid() #Add grid
    
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

# Some more plots
plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,W_1,"b")
plt.plot(rr,W_2,"g")
plt.ylabel(r" $W$ $[m/s]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out"])
plt.title("Relative Velocity")
plt.grid()

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,V_1,"b")
plt.plot(rr,V_2,"g")
plt.plot(rr,V_3,"r")
plt.ylabel(r"$V$ $[m/s]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out"])
plt.title("Absolute Velocity")
plt.grid()

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,p_1,"b")
plt.plot(rr,p_2,"g")
plt.plot(rr,p_3,"r")
plt.ylabel(r"$p$ $[Pa]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out"])
plt.title("Static Pressure")
plt.grid()

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,rho_1,"b")
plt.plot(rr,rho_2,"g")
plt.plot(rr,rho_3,"r")
plt.ylabel(r"$\rho$ $[kg/m^3]$")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out"])
plt.title("Density")
plt.grid()

plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,180/np.pi * np.array(alpha_1),"b")
plt.plot(rr,180/np.pi * np.array(alpha_2),"g")
plt.plot(rr,180/np.pi * np.array(alpha_3),"r")
plt.ylabel(r"$\alpha$ [deg]")
plt.xlabel(r"$r \  [m]$")
plt.legend(["Rotor In","Rotor Out","Stator Out"])
plt.title("Absolute Flow Angle")
plt.grid()


plt.figure(figsize=(6, 5), dpi=80)
plt.plot(rr,chi)
plt.ylabel(r"$\chi$")
plt.xlabel(r"$r \  [m]$")
plt.title("Reaction Degree")
plt.grid()

# # This should be constant if a free vortex distribution is used
# plt.figure(figsize=(6, 5), dpi=80)
# plt.plot(rr,L_eul)
# plt.title("Euler Work [J/kg]")

plt.show()
