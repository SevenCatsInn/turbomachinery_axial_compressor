from numpy import sqrt, arctan, tan, pi, cos, sin
import numpy as np

Norm = lambda x : sqrt(x[0]**2 + x[1]**2)


print("--------------- MEAN LINE DESIGN STAGE 1 ---------------")
# Thermophysical properties
c_p = cp = 1005  # Constant pressure specific heat [J/(kg K)]
gamma = 1.4 # Specific heat ratio
c_v = c_p/gamma
R = c_p * (gamma-1)/gamma # Gas constant [J/(kg K)]

#Input
mdot=100 #mass flow rate
Pt1=99439 # pressure [bar]
Tt1=300 #inlet temperature [K]
beta=1.22 #compression ratio


## Non dimensional quantities

## axial compressor
#vavra: get reaction degree and flow coefficient to get machimum efficiency
phi=0.88 #from slide 10 achial compressors
chi=0.5 #reaction degree
psi=0.34 #from first graph slide 12
Rm=0.32 #mean line radius
efficiency_TT=0.905
eta_S = 0.92
eta_R = 0.92

#determine loading
L_is=cp*Tt1*(beta**((gamma-1)/gamma)-1)
L_eul=L_is/efficiency_TT
lamda=psi*2
Um=sqrt(L_eul/psi)
omega=Um/Rm
RPM=omega*30/pi
V1a=phi*Um
V1t = Um*(1-chi-lamda/4)
V1=[V1a, V1t]
alpha1= arctan(V1t/V1a) * 180/pi
V1_mag=Norm(V1)
W1a=V1a
W1t=V1t-Um
W1=[W1a, W1t]
W1_mag=Norm(W1)
beta1=arctan(W1[1]/W1[0])

T1=Tt1-V1_mag**2/(2*cp)
M1=V1_mag/sqrt(gamma*R*T1)
p1=Pt1*(1+(gamma-1)/2*M1**2)**((-gamma)/(gamma-1))
rho1=p1/(R*T1)

b1=mdot/(rho1*V1a*2*pi*Rm)
b2=0.9*b1

#quantities at station 2 (after rotor)
V2t=L_eul/Um+V1t
Tt2=L_eul/cp+Tt1
W2_mag=sqrt(W1_mag**2-chi*L_eul*2)
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
chi_new=(V1a**2-V2a**2+V1t**2-V2t**2+2*Um*(V2t-V1t))/(2*L_eul)
chi=(W1_mag**2-W2_mag_new**2)/(2*L_eul) 

print("")
print("Relative deflection in rotor")
print(cos(beta1) / cos(beta2), "> 0.72 ?") 

# Mean line design for the stator
alpha2 = arctan(V2t/V2a)
alpha3 = 30 * pi/180 # Design choice, stage 2 inlet

print("")
print("Absolute deflection in stator")
print(cos(alpha3) / cos(alpha2), "> 0.72 ?") 

V0a = V1a
V0t = 0
Tt0 = Tt1 # Imposed by thermodynamics, no work in stator


err = 1e10
tol = 0.0001
iter = 0

while abs(err)>tol:
    V0_mag=sqrt(V0a**2+V0t**2)
    T0 = Tt0 - V0_mag**2/(2*cp)
    T0is= T1 - 1/eta_S * (T1-T0)
    p0=(T0is/Tt1)**(gamma/(gamma-1))*Pt1
    T1is=T0 + eta_S*(T1-T0)
    rho0=p0/(R*T0)
    V0a_new=mdot/(rho0*2*pi*b1*Rm)
    err=abs(V0a_new-V0a)
    V0a=V0a_new
    M0 = V0_mag / sqrt(gamma * R * T0)
    Pt0 = p0 * (1 + (gamma-1)/2 * M0**2)**(gamma/(gamma-1))
    iter=iter+1

print(Pt0, "Pt0 = 100000 Pa ?")


V3a = V2a
V3t = V3a * tan(alpha3)
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
    V3a_new=mdot/(rho3*2*pi*b2*Rm)
    err=abs(V3a_new-V3a)
    V3a=V3a_new
    V3t = V3a * tan(alpha3)
    M3 = V3_mag / sqrt(gamma * R * T3)
    Pt3 = p3 * (1 + (gamma-1)/2 * M3**2)**(gamma/(gamma-1))
    iter=iter+1
V3 = [V3a, V3t]

# Change variables names for radial equilibrium script
T_0m  = T0
p_0m  = p0
V_a0m = V0a
V_t0m = V0t

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
p_t3  = Pt3
T_t3  = Tt3
V_a3m = V3a
V_t3m = V3t

b_1 = b1
R_m = Rm
rpm = RPM

print("")
print("\u03C7, \u03A6, \u03A8 = ", chi, phi, psi)

# General Whirl Design
# a * R_m - b / R_m = V_t1m
# a * R_m + b / R_m = V_t2m

n = 1
matA = np.array([[R_m**n, -1 / R_m], 
                 [R_m**n,  1 / R_m]])

vecB = np.array([[V_t1m],
                 [V_t2m]])

x = np.linalg.solve(matA,vecB)

a = (x[0])[0]
b = (x[1])[0]

rtip = Rm + b_1/2
rhub = Rm - b_1/2

