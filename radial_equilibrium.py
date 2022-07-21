### Radial equilibrium script ###

from sympy import *
from sympy import init_printing
import numpy as np
import matplotlib.pyplot as plt

init_printing() 


r = Symbol("r", positive=True) # Radius variable declaration


## Data from mean line design

# Geometry
R_m = 0.3            # Mean radius         [m]
b_1 = 0.2914         # Inlet blade height  [m]
R_h = R_m - b_1 / 2  # Hub Radius          [m]   
R_t = R_m + b_1 / 2  # Tip Radius          [m]  

#Velocities
V_a1m =  157.46  # Mean radius inlet axial velocity [m/s]
V_t1m =  0       # Mean radius inlet tang. velocity [m/s]
V_a2m =  161.87  # Mean radius outl. axial velocity [m/s]
V_t2m =  68.9    # Mean radius outl. tang. velocity [m/s]

rpm = 6265              # Rotations per minute [giri/mins]
omega = 2 * pi * rpm/60 # Angular velocity     [rad/s]
U = omega * r           # Peripheral velocity  [m/s]

# Thermodynamics
T_1m =  287.6526 # Mean radius inlet static temperature [K]
p_1m =  86383    # Mean radius inlet static pressure    [Pa]
T_2m =  298.09   # Mean radius outl. static temperature [K]
p_2m =  96659    # Mean radius outl. static pressure    [Pa]


# Thermophysical properties
c_p = 1005  # Constant pressure specific heat [J/(kg K)]
gamma = 1.4 # Specific heat ratio
c_v = c_p/gamma
R = c_p * (gamma-1)/gamma # Gas constant [J/(kg K)]

# Input data
T_t1 = 300 # [K]     --> f(r)
p_t1 = 100000 # [Pa] --> f(r)
m_dot_req = 100 # Required mass flow [kg/s]

# Computed Quantities
h_t1 = c_p * T_t1 # Total enthalpy [J/kg]

ds_1  = -R * diff(p_t1,r) / p_t1 + c_p * diff(T_t1,r) / T_t1 # Derivative over r of entropy [J/(kg K)]

# Set the design choice for tangential velocity distribution in the radial direction
V_t1 = R_m * V_t1m / r # e.g. Free vortex distribution r * V_t = const

print("")
print("########## INLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 0.001
iter= 1

# This loop can be eliminated by varying b_1 to accomodate for m_dot_req
while abs(err) > tol: 
    
    V_a1 = Function('V_a1') # Declare the function to solve in the O.D.E.
    
    # Non isentropic radial equilibrium in station 1
    nisre_1 = Eq(V_a1(r)*V_a1(r).diff(r) + V_t1 / r * diff(r*V_t1,r) + T_1m * ds_1, diff(h_t1,r) ) 
    sol_nisre_1 = dsolve(nisre_1, ics = {V_a1(R_m):V_a1m})
    V_a1 = sol_nisre_1.rhs # Extract the value of the NISRE solution, assign to V_a1(r)

    # Prints
    print("")
    print("---Iteration no. " + str(iter))
    print("N.I.S.R.E. in Station 1")
    pprint(nisre_1)
    print("Solution:")
    pprint(sol_nisre_1)
    print("")

    ## Compute quantities in station 1 --> f(r)

    # Kinematics
    V_1 = sqrt(V_a1**2 + V_t1**2)
    alpha_1 = atan(V_t1/V_a1)
    W_t1 = V_t1 - U
    W_a1 = V_a1
    W_1 = sqrt(W_t1**2 + W_a1**2)
    beta_1 = atan(W_t1/W_a1)

    # Thermodynamics
    T_1 = T_t1 - V_1**2 / (2 * c_p)
    p_1 = p_1m * (T_1 / T_1m)**(gamma/(gamma-1))
    rho_1 = p_1 / (R*T_1)
    M_1  = V_1 / sqrt(gamma * R * T_1)
    M_1r = W_1 / sqrt(gamma * R * T_1)
    p_t1  = p_1*(1 + (gamma-1) / 2 * M_1**2  ) ** (gamma/(gamma-1))
    p_t1r = p_1*(1 + (gamma-1) / 2 * M_1r**2 ) ** (gamma/(gamma-1))

    ## Compute the mass flow at the inlet

    rr = np.linspace(R_h, R_t, 100) # Extremes of integration

    integrand_1 = [] # 2 * pi * r * V_a1 * rho_1
    for radius in rr:
        integrand_1.append(2 * np.pi * radius * (V_a1.subs(r,radius)).evalf() * (rho_1.subs(r,radius)).evalf() ) 

    m_dot_trap = np.trapz(integrand_1,rr)
    
    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a1m = V_a1m * (1 + err) # New axial velocity
    iter += 1
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("err = "+ str(err))


print("")
print("########## OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 0.001
iter= 1


# Inputs
# Discretization
pts = 101  # Total number of points across the radius, uneven!
rr = np.linspace(R_h, R_t, pts ) # Extremes of integration
deltaR = (R_t - R_h)/ (pts - 1)
mean_index = pts//2  # Index of the various list corresponding to mean radius quantities


# Entropy inputs, NOTE: absolute values are meaningless
omega_loss = 0.0
s_1 = 0
s_2 = list(s_1 * ones(1,pts)) # Initial radial entropy distribution in 2
ds_2 = list(ds_1 * ones(1,pts)) # Initial assumption, negligible s variation

V_t2 = V_t2m * R_m / r # Outlet tangential velocity distribution (e.g. free vortex)
h_t2 = h_t1 + U * (V_t2 - V_t1)
T_t2 = h_t2 / c_p

T_2 = list(T_2m * ones(1,pts) )

dh_t2_lst = []
T_t2_lst = []
V_t2_lst = []
drV_t2_lst = []

for radius in rr:
    dh_t2_lst.append(diff(h_t2,r).subs(r,radius).evalf())
    T_t2_lst.append(T_t2.subs(r,radius).evalf())
    V_t2_lst.append(V_t2.subs(r,radius).evalf())
    drV_t2_lst.append(diff(r*V_t2,r).subs(r,radius).evalf())

dh_t2 = dh_t2_lst
T_t2 = T_t2_lst
V_t2 = V_t2_lst
drV_t2 = drV_t2_lst

# This loop can be avoided using flaired blades b_2 != b_1
while abs(err) > tol:

    V_a2 = list(zeros(1,pts))
    V_a2[mean_index] = V_a2m # Initial value for forward integration
    dV_a2 = list(zeros(1,pts))

    # Variable index to generate the lists 
    # TODO: Rework the upper and lower portion into one line with for q in [m_i +- j]
    for j in list(range(0,mean_index)):
        # Upper portion
        
        dV_a2[mean_index + j] = 1 / V_a2[mean_index + j] * ( dh_t2[mean_index + j] - T_2[mean_index + j] * ds_2[mean_index + j] - V_t2[mean_index + j] / rr[mean_index + j] * drV_t2[mean_index + j] )
        
        dV_a2[mean_index - j] = 1 / V_a2[mean_index - j] * ( dh_t2[mean_index - j] - T_2[mean_index - j] * ds_2[mean_index - j] - V_t2[mean_index - j] / rr[mean_index - j] * drV_t2[mean_index - j] )
        
        V_a2[mean_index + j + 1] = V_a2[mean_index + j] + dV_a2[mean_index + j] * deltaR 
        
        V_a2[mean_index - j - 1] = V_a2[mean_index - j] - dV_a2[mean_index - j] * deltaR 

    V_2 = list(zeros(1,pts))
    alpha_2 = list(zeros(1,pts))
    W_t2 = list(zeros(1,pts))
    W_a2 = list(zeros(1,pts))
    W_2 = list(zeros(1,pts))
    beta_2 = list(zeros(1,pts))
    p_2 = list(zeros(1,pts))
    rho_2 = list(zeros(1,pts))
    M_2 = list(zeros(1,pts))
    M_2r = list(zeros(1,pts))
    p_t2 = list(zeros(1,pts))
    p_t2r = list(zeros(1,pts))
    integrand_2 = list(zeros(1,pts))

    for j in list(range(pts)): # Compute quantities along the radius
        # Kinematics
        V_2[j] = sqrt(V_a2[j]**2 + V_t2[j]**2)
        alpha_2[j] = atan(V_t2[j]/V_a2[j])
        W_t2[j] = V_t2[j] - U.subs(r,rr[j])
        W_a2[j] = V_a2[j]
        W_2[j] = sqrt(W_t2[j]**2 + W_a2[j]**2)
        beta_2[j] = atan(W_t2[j]/W_a2[j])

        # Thermodynamics
        T_2[j] = T_t2[j] - V_2[j]**2 / (2 * c_p)
        p_2[j] = p_2m * (T_2[j] / T_2m)**(gamma/(gamma-1)) * exp(- (s_2[j] - s_2[mean_index]) / R)
        rho_2[j] = p_2[j] / (R*T_2[j])
        M_2[j]  = V_2[j] / sqrt(gamma * R * T_2[j])
        M_2r[j] = W_2[j] / sqrt(gamma * R * T_2[j])
        p_t2[j] = p_2[j]*(1 + (gamma-1) / 2 * M_2[j]**2 ) ** (gamma/(gamma-1))

        integrand_2[j] = 2 * np.pi * rr[j] * rho_2[j] * V_a2[j] 
        
        #Evaluate the q.ties in section 1 (expressions) at the current radius
        # tmp = overwritten at every iteration, no need for a new array for _1 quantities
        p_1_tmp = p_1.subs(r,rr[j]).evalf()
        p_t1r_tmp = p_t1r.subs(r,rr[j]).evalf()
        
        p_t2r[j] = p_t1r_tmp - omega_loss * (p_t1r_tmp - p_1_tmp)

        # ENTROPY EVALUATION

        s_2[j]  = s_1 - R * ln(p_t2r[j] / p_t1r_tmp)

    ds_2 = []
    [ds_2.append( (s_2[i-1] - s_2[i+1]) / (2*deltaR) ) for i in range(1,len(s_2) - 1)]

    ds_2 = [(s_2[1] - s_2[0]) / deltaR] + ds_2 + [(s_2[-2] - s_2[-1]) / deltaR]

    # print(ds_2)
    # plot(p_t1r, (r,R_h,R_t))

    m_dot_trap = np.trapz(integrand_2, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a2m = V_a2m*(1 + err) # New axial velocity
    iter += 1
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("err = "+ str(err))

# L_eul = U * (V_t2 - V_t1)
# chi = (W_1**2 - W_2**2)/(2*L_eul)


# # Plot inlet and outlet velocity triangles at hub, mean radius and tip
# # P stands for plotting

# fig, axs = plt.subplots(3,1, sharex=True, sharey=True, figsize=(3, 6), dpi=65) # Create figure

# j = 0 # Index used to move through the subplots
# for i in [R_h, R_m, R_t]:
#     # Evaluate the quantities to plot on the desired radius
#     U_P   = float(U.subs(r,   i).evalf())
#     V_a1P = float(V_a1.subs(r,i).evalf())
#     V_t1P = float(V_t1.subs(r,i).evalf())
#     W_a1P = float(W_a1.subs(r,i).evalf())
#     W_t1P = float(W_t1.subs(r,i).evalf())
#     V_a2P = float(V_a2.subs(r,i).evalf())
#     V_t2P = float(V_t2.subs(r,i).evalf())
#     W_a2P = float(W_a2.subs(r,i).evalf())
#     W_t2P = float(W_t2.subs(r,i).evalf())

#     # axs[j].grid() #Add grid
    
#     #Plot inlet and outlet triangles
#     axs[j].quiver([0,U_P - V_t1P, U_P - V_t1P] , [0,V_a1P,V_a1P] , [U_P,V_t1P,W_t1P] , [0,-V_a1P,-W_a1P] , angles='xy',scale_units='xy', scale=1.0, color=["black","blue","blue"])
#     axs[j].quiver([0,U_P - V_t2P, U_P - V_t2P] , [0,V_a2P,V_a2P] , [U_P,V_t2P,W_t2P] , [0,-V_a2P,-W_a2P] , angles='xy',scale_units='xy', scale=1.,  color=["black","red","red"])
    
#     axs.flat[j].set_xlim(-50, 300) #Set the limits for the x axis
#     axs.flat[j].set_ylim(-5, 200)  #Set the limits for the y axis
    
#     axs[j].set_aspect('equal') #Equal aspect ratio axes

#     j = j+1
# plt.show()

