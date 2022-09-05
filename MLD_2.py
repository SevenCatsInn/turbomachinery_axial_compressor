from numpy import sqrt, arctan, tan, pi, cos

print("")
print("--------------- MEAN LINE DESIGN STAGE 2 ---------------")

Norm = lambda x : sqrt(x[0]**2 + x[1]**2)


# Thermophysical properties
c_p = cp = 1005  # Constant pressure specific heat [J/(kg K)]
gamma = 1.4 # Specific heat ratio
c_v = c_p/gamma
R = c_p * (gamma-1)/gamma # Gas constant [J/(kg K)]

#Input
mdot=100 #mass flow rate
Pt3=p_t3 # pressure [bar]
Tt3=T_t3 #inlet temperature [K]
beta=1.457/beta #compression ratio



## Non dimensional quantities <3 <3

## achi2al compressor
#vavra: get reaction degree and flow coefficient to get machi2mum efficiency
chi2=0.5 #reaction degree
efficiency_TT=0.905
eta_S = 0.92
eta_R = 0.92
Um = U[mean_index]

#determine loading
L_is=cp*Tt3*(beta**((gamma-1)/gamma)-1)
# L_is = 9000
L_eul=L_is/efficiency_TT

psi = L_eul / Um**2
lamda=psi*2
omega=rpm * pi / 30

V3a=V_a3m
phi = V3a / Um
V3t=V_t3m

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
rho3=p3/(R*T3)



#quantities at station 4 (after rotor)
V4t=L_eul/Um + V3t
Tt4=L_eul/cp + Tt3
W4_mag=sqrt(W3_mag**2-chi2*L_eul*2)
W4t=V4t-Um
W4a=sqrt(W4_mag**2-W4t**2)
W4=[W4a, W4t]

V4a=W4a
err=20
i=0

while abs(err) > 10**(-4):
    V4_mag=sqrt(V4a**2+V4t**2)
    T4=Tt4-V4_mag**2/(2*cp)
    T4is=T3 + eta_R*(T4-T3)
    p4=(T4is/Tt3)**(gamma/(gamma-1))*Pt3
    rho4=p4/(R*T4)
    V4a_new=mdot/(rho4*2*pi*b2*Rm)
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
# chi2_new=(V1a**2-V2a**2+V1t**2-V2t**2+2*Um*(V2t-V1t))/(2*L_eul)
chi2=(W3_mag**2-W4_mag_new**2)/(2*L_eul) 


# Mean line design for the stator

alpha4 = arctan(V4t/V4a)
alpha5 = 30 * pi/180# Design choice

print("")
print("Relative deflection in rotor")
print(cos(beta3) / cos(beta4), "> 0.72 ?") 



print("")
print("Absolute deflection in stator")
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
    V5a_new=mdot/(rho5*2*pi*b2*Rm)
    V5t = V5a * tan(alpha5)
    err=abs(V5a_new-V5a)
    V5a=V5a_new
    iter=iter+1
V5 = [V5a, V5t]

# Change variables names for radial equilibrium script


b_2 = b2
T_4m  = T4
p_4m  = p4
V_a4m = V4a
V_t4m = V4t
T_5m  = T5
p_5m  = p5
V_a5m = V5a
V_t5m = V5t


print("")
print("\u03C7, \u03A6, \u03A8 = ", chi2, phi, psi)
# General Whirl Design
# a * R_m - b / R_m = V_t1m
# a * R_m + b / R_m = V_t2m


n = 0.3
matA2 = np.array([[R_m**n, -1 / R_m], 
                  [R_m**n,  1 / R_m]])

vecB2 = np.array([[V_t3m],[V_t4m]])


x2 = np.linalg.solve(matA2,vecB2)

a22 = (x2[0])[0]
b22 = (x2[1])[0]
