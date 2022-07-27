from numpy import sqrt, arctan, tan, pi, cos

Norm = lambda x : sqrt(x[0]**2 + x[1]**2)

exec(open("./radial_equilibrium.py").read())
# Thermophysical properties
c_p = cp = 1005  # Constant pressure specific heat [J/(kg K)]
gamma = 1.4 # Specific heat ratio
c_v = c_p/gamma
R = c_p * (gamma-1)/gamma # Gas constant [J/(kg K)]

#Input
mdot=100 #mass flow rate
Pt1=p_t3[mean_index] # pressure [bar]
Tt1=T_t3[mean_index] #inlet temperature [K]
beta=1.45/beta #compression ratio
Rm=0.30 #mean line radius
rho1=Pt1/(R*Tt1)
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
L_is=cp*Tt1*(beta**((gamma-1)/gamma)-1)
L_eul=L_is/efficiency_TT

psi = L_eul / Um**2
lamda=psi*2
omega=rpm * pi / 30
V1a=V_a3[mean_index]
phi = V1a / Um
V1t=V_t3[mean_index]

V1=[V1a, V1t]
V1_mag=Norm(V1)
W1a=V1a
W1t=V1t-Um
W1=[W1a, W1t]
W1_mag=Norm(W1)
beta1=arctan(W1[1]/W1[0])

T1=Tt1-V1_mag**2/(2*cp)
M1=V1_mag/sqrt(gamma*R*T1)
p1=Pt1*(1+(gamma-1)/2*M1**2)**((-gamma)/(gamma-1))
b1=mdot/(rho1*V1a*2*pi*Rm)

#quantities at station 2 (after rotor)
V2t=L_eul/Um+V1t
Tt2=L_eul/cp+Tt1
W2_mag=sqrt(W1_mag**2-xi*L_eul*2)
W2t=V2t-Um
W2a=sqrt(W2_mag**2-W2t**2)
W2=[W2a, W2t]

V2a=W2a
err=20
i=0

while abs(err) > 10**(-4):
    V2_mag=sqrt(V2a**2+V2t**2)
    T2=Tt2-V2_mag**2/(2*cp)
    T2is=T1+eta_R*(T2-T1)
    p2=(T2is/Tt1)**(gamma/(gamma-1))*Pt1
    rho2=p2/(R*T2)
    V2a_new=mdot/(rho2*2*pi*b1*Rm)
    err=abs(V2a_new-V2a)
    V2a=V2a_new
    M2 = V2_mag / sqrt(gamma * R * T2)
    Pt2 = p2 * (1 + (gamma-1)/2 * M2**2)**(gamma/(gamma-1))
    i=i+1


W2a=V2a
W2=[W2a, W2t]
W2_mag_new=sqrt(W2a**2+W2t**2)
beta2=arctan(W2t/W2a)
V2=[V2a, V2t]
xi_new=(V1a**2-V2a**2+V1t**2-V2t**2+2*Um*(V2t-V1t))/(2*L_eul)
xi=(W1_mag**2-W2_mag_new**2)/(2*L_eul) 


# Mean line design for the stator

alpha_2 = arctan(V2t/V2a)
alpha_3 = 0 * pi/180# Design choice

print(cos(alpha_3) / cos(alpha_2), "> 0.72 ?") 

V3a = V2a
V3t = V3a * tan(alpha_3)
Tt3 = Tt2 # Imposed by thermodynamics, no work in stator

err = 1e10
tol = 0.0001
iter = 0

while abs(err)>tol:
    V3_mag=sqrt(V3a**2+V3t**2)
    T3=Tt3-V3_mag**2/(2*cp)
    T3is=T2 + eta_S*(T3-T2)
    p3=(T3is/Tt2)**(gamma/(gamma-1))*Pt2
    rho3=p3/(R*T3)
    V3a_new=mdot/(rho3*2*pi*b1*Rm)
    err=abs(V3a_new-V3a)
    V3a=V3a_new
    iter=iter+1
V3 = [V3a, V3t]

# Change variables names for radial equilibrium script
T_1m  = T1
p_1m  = p1
V_a1m = V1a
V_t1m = V1t
T_2m  = T2
p_2m  = p2
V_a2m = V2a
V_t2m = V2t
T_3m  = T3
p_3m  = p3
V_a3m = V3a
V_t3m = V3t
b_1 = b1
R_m = Rm
rpm = RPM

# print(beta)
# print("phi,psi,chi = ", phi,psi,xi_new)