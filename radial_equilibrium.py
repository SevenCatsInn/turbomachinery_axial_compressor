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

ds_1dr  = -R * diff(p_t1,r) / p_t1 + c_p * diff(T_t1,r) / T_t1 # Derivative over r of entropy [J/(kg K)]

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
    nisre_1 = Eq(V_a1(r)*V_a1(r).diff(r) + V_t1 / r * diff(r*V_t1,r) + T_1m * ds_1dr, diff(h_t1,r) ) 
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
    p_t1 = p_1*(1 + (gamma-1) / 2 * M_1**2 ) ** (gamma/(gamma-1))

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

# quit()

print("")
print("########## OUTLET ##########")

err = 1e10 # Inital value to enter the loop, meaningless
tol = 0.001
iter= 1


# Inputs
ds_2dr = ds_1dr # Initial assumption, negligible s variation
s_2m = 0   # Reference arbitrary value at mean radius
s_2 = s_2m # Initial radial entropy distribution in 2
V_t2 = V_t2m * R_m / r # Outlet tangential velocity distribution (e.g. free vortex)

h_t2 = h_t1 + U * (V_t2 - V_t1)
T_t2 = h_t2 / c_p

# This loop can be avoided using flaired blades b_2 != b_1
while abs(err) > tol:
    V_a2 = Function('V_a2') # Declare the function to solve in the O.D.E.

    # Non isentropic radial equilibrium in station 2
    nisre_2 = Eq(V_a2(r)*V_a2(r).diff(r) + V_t2 / r * diff(r*V_t2,r) + T_2m * ds_2dr, diff(h_t2,r) ) 
    sol_nisre_2 = dsolve(nisre_2, ics = {V_a2(R_m):V_a2m})
    V_a2 = sol_nisre_2.rhs # Extract the value of the NISRE solution, assign to V_a2(r)

    # Prints
    print("")
    print("---Iteration no. " + str(iter))
    print("N.I.S.R.E. in Station 2")
    # pprint(nisre_2)
    print("Solution:")
    pprint(sol_nisre_2)
    print("")

    ## Compute quantities in station 2 --> f(r)

    # Kinematics
    V_2 = sqrt(V_a2**2 + V_t2**2)
    alpha_2 = atan(V_t2/V_a2)
    W_t2 = V_t2 - U
    W_a2 = V_a2
    W_2 = sqrt(W_t2**2 + W_a2**2)
    beta_2 = atan(W_t2/W_a2)

    # Thermodynamics
    T_2 = T_t2 - V_2**2 / (2 * c_p)
    p_2 = p_2m * (T_2 / T_2m)**(gamma/(gamma-1)) * exp(- (s_2 - s_2m) / R)
    rho_2 = p_2 / (R*T_2)
    M_2  = V_2 / sqrt(gamma * R * T_2)
    M_2r = W_2 / sqrt(gamma * R * T_2)
    p_t2 = p_2*(1 + (gamma-1) / 2 * M_2**2 ) ** (gamma/(gamma-1))
    
    # ENTROPY EVALUATION

    ds2_dr = -R * diff(p_t2,r) / p_t2 + c_p * diff(T_t2,r) / T_t2 # New entropy derivative
    
    points = 8 // 2  * 2 # Number of points in which we compute s_2, rounded to the nearest even integer
    s_2l = [ s_2m ] # From R_m to tip
    s_2u = [ s_2m ] # From hub to R_m
    
    deltaR = (R_t - R_h)/ points
    i = 0
    for radius in np.linspace(R_m, R_t, points // 2):
        s_2u.append(s_2u[i] + deltaR * (ds2_dr.subs(r,radius)).evalf())
        i=i+1
    i = 0
    for radius in np.linspace(R_m, R_h, points // 2):
        s_2l.append(s_2l[i] - deltaR * (ds2_dr.subs(r,radius)).evalf())
        i=i+1
    s_2l = list(reversed(s_2l)) # The list was built in reverse, invert it

    s_2_points = s_2l + s_2u[1:] # Sum the two lists into one
    
    # Generate a spline f(r) to keep on using symbolic computations below
    s_2 = interpolating_spline(1, r, np.linspace(R_h, R_t, points + 1), s_2_points)

    ## Compute the mass flow at the inlet
    rr = np.linspace(R_h, R_t, 100) # Extremes of integration
    
    integrand_2 = [] # 2 * pi * r * V_a2 * rho_2
    for radius in rr:
        integrand_2.append(2 * np.pi * radius * (V_a2.subs(r,radius)).evalf() * (rho_2.subs(r,radius)).evalf() ) 


    m_dot_trap = np.trapz(integrand_2, rr)

    err  = 1 - m_dot_trap/m_dot_req # Error
    V_a2m = V_a2m*(1 + err) # New axial velocity
    iter += 1
    print("mass flow = "+ str(m_dot_trap) + " [kg/s]")
    print("err = "+ str(err))

L_eul = U * (V_t2 - V_t1)
chi = (W_1**2 - W_2**2)/(2*L_eul)


# Plot input and output velocity triangles at hub, mean radius and tip
fig, axs = plt.subplots(3,1, sharex=True, sharey=True)
# P stands for plotting

j = 0
for i in [R_h, R_m, R_t]:
    U_P   = float(U.subs(r,   i).evalf())
    V_a1P = float(V_a1.subs(r,i).evalf())
    V_t1P = float(V_t1.subs(r,i).evalf())
    W_a1P = float(W_a1.subs(r,i).evalf())
    W_t1P = float(W_t1.subs(r,i).evalf())
    V_a2P = float(V_a2.subs(r,i).evalf())
    V_t2P = float(V_t2.subs(r,i).evalf())
    W_a2P = float(W_a2.subs(r,i).evalf())
    W_t2P = float(W_t2.subs(r,i).evalf())

    
    axs[j].grid()
    
    axs[j].quiver([0,U_P - V_t1P, U_P - V_t1P] , [0,V_a1P,V_a1P] , [U_P,V_t1P,W_t1P] , [0,-V_a1P,-W_a1P] , angles='xy',scale_units='xy', scale=1.0)
    axs[j].quiver([0,U_P - V_t2P, U_P - V_t2P] , [0,V_a2P,V_a2P] , [U_P,V_t2P,W_t2P] , [0,-V_a2P,-W_a2P] , angles='xy',scale_units='xy', scale=1.0)
    
    axs.flat[j].set_xlim(-50, 300)
    
    axs[j].set_aspect('equal')

    j = j+1

plt.show()