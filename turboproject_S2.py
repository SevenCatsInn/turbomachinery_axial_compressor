from numpy import sqrt, arctan, tan, pi, cos

Norm = lambda x : sqrt(x[0]**2 + x[1]**2)

# Thermophysical properties
c_p = cp = 1005  # Constant pressure specific heat [J/(kg K)]
gamma = 1.4 # Specific heat ratio
c_v = c_p/gamma
R = c_p * (gamma-1)/gamma # Gas constant [J/(kg K)]

#Input
mdot=100 #mass flow rate
Pt3=p_t3[mean_index] # pressure [bar]
Tt3=T_t3[mean_index] #inlet temperature [K]
beta=1.45/beta #compression ratio
Rm=0.30 #mean line radius
rho3=Pt3/(R*Tt3)
Q=mdot/rho1


## Non dimensional quantities <3 <3

## axial compressor
#vavra: get reaction degree and flow coefficient to get maximum efficiency
xi=0.5 #reaction degree
efficiency_TT=0.92
eta_S = 0.92
eta_R = 0.92
Um = U[mean_index]

#determine loading
L_is=cp*Tt3*(beta**((gamma-1)/gamma)-1)
L_eul=L_is/efficiency_TT

psi = L_eul / Um**2
lamda=psi*2
omega=rpm * pi / 30

V3a=V_a3[mean_index]
phi = V3a / Um
V3t=V_t3[mean_index]

V3=[V3a, V3t]
V1_mag=Norm(V3)
W3a=V3a
W3t=V3t-Um
W3=[W3a, W3t]
W3_mag=Norm(W3)
beta3=arctan(W3t/W3a)

T3=Tt3-V3_mag**2/(2*cp)
M3=V3_mag/sqrt(gamma*R*T3)
p3=Pt3*(1+(gamma-1)/2*M3**2)**((-gamma)/(gamma-1))

#quantities at station 4 (after rotor)
V4t=L_eul/Um + V3t
Tt4=L_eul/cp + Tt3
W4_mag=sqrt(W3_mag**2-xi*L_eul*2)
W4t=V4t-Um
W4a=sqrt(W4_mag**2-W4t**2)
W4=[W4a, W4t]

V4a=W4a
err=20
i=0

while abs(err) > 10**(-4):
    V4_mag=sqrt(V4a**2+V4t**2)
    T4=Tt4-V4_mag**2/(2*cp)
    T4is=T4+eta_R*(T4-T3)
    p4=(T4is/Tt3)**(gamma/(gamma-1))*Pt3
    rho4=p4/(R*T4)
    V4a_new=mdot/(rho4*2*pi*b1*Rm)
    err=abs(V4a_new-V4a)
    V4a=V4a_new
    M4 = V4_mag / sqrt(gamma * R * T4)
    Pt4 = p4 * (1 + (gamma-1)/2 * M4**2)**(gamma/(gamma-1))
    i=i+1


W4a=V4a
W4=[W4a, W4t]
W4_mag_new=Norm(W4)
beta4=arctan(W4t/W4a)
V4=[V4a, V4t]
# xi_new=(V1a**2-V2a**2+V1t**2-V2t**2+2*Um*(V2t-V1t))/(2*L_eul)
xi=(W3_mag**2-W4_mag_new**2)/(2*L_eul) 


# Mean line design for the stator

alpha4 = arctan(V4t/V4a)
alpha5 = 0 * pi/180# Design choice

print(cos(alpha5) / cos(alpha4), "> 0.72 ?") 

V5a = V4a
V5t = V5a * tan(alpha5)
Tt5 = Tt4 # Imposed by thermodynamics, no work in stator

err = 1e10
tol = 0.0001
iter = 0

while abs(err)>tol:
    V5_mag=sqrt(V5a**2+V5t**2)
    T5=Tt5-V5_mag**2/(2*cp)
    T5is=T4 + eta_S*(T5-T4)
    p5=(T5is/Tt4)**(gamma/(gamma-1))*Pt4
    rho5=p5/(R*T5)
    V5a_new=mdot/(rho5*2*pi*b1*Rm)
    err=abs(V5a_new-V5a)
    V5a=V5a_new
    iter=iter+1
V5 = [V5a, V5t]

# Change variables names for radial equilibrium script


# print(beta)
# print("phi,psi,chi = ", phi,psi,xi_new)

T_4m  = T4
p_4m  = p4
V_a4m = V4a
V_t4m = V4t
T_5m  = T5
p_5m  = p5
V_a5m = V5a
V_t5m = V5t




# print(T_4m ,p_4m ,V_a4m,V_t4m)