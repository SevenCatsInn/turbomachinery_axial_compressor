from numpy import sqrt, arctan, tan, pi, cos

Norm = lambda x : sqrt(x[0]**2 + x[1]**2)

#compressor project
mdot=100 #mass flow rate
Pt1=100000 # pressure [bar]
Tt1=300 #inlet temperature [K]
beta=1.15 #compression ratio
cp=1004
gamma=1.4
R=8314.46/28.84
rho1=Pt1/(R*Tt1)
Q=mdot/rho1


## Non dimensional quantities
#centrifugal compressor
#chose specific rotational speed and specific diameter from Balje Diagram
omega_s_centr=2
D_s_centr=2
efficiency_TT_centr=0.85

L_is_centr=cp*Tt1*(beta**((gamma-1)/gamma)-1) #(total enthalpy jump(increase))
L_eul_centr=L_is_centr/efficiency_TT_centr
#compute rotational speed
omega_centr=omega_s_centr*L_is_centr**(0.75)/(sqrt(Q)) #[rad/s]
RPM_centr=omega_centr*30/pi
D_centr=D_s_centr*sqrt(Q)/(L_is_centr**0.25)  #2 stages?

## axial compressor
#vavra: get reaction degree and flow coefficient to get maximum efficiency
phi_ax=0.8 #from slide 10 axial compressors
xi_ax=0.5 #reaction degree
efficiency_TT_ax=0.905
#determine loading
psi_ax=0.35 #from first graph slide 12
L_is_ax=cp*Tt1*(beta**((gamma-1)/gamma)-1)
L_eul_ax=L_is_ax/efficiency_TT_ax
lamda_ax=psi_ax*2
#psi=L_eul/U**2
Um_ax=sqrt(L_eul_ax/psi_ax)
Rm_ax=0.30 #mean line radius
omega_ax=Um_ax/Rm_ax
RPM_ax=omega_ax*30/pi
V1a=phi_ax*Um_ax
V1t=0
V1=[V1a, V1t]
V1_mag=Norm(V1)
W1a=V1a
W1t=V1t-Um_ax
W1=[W1a, W1t]
W1_mag=Norm(W1)
beta=arctan(W1[1]/W1[0])

T1=Tt1-V1_mag**2/(2*cp)
M1=V1_mag/sqrt(gamma*R*T1)
p1=Pt1*(1+(gamma-1)/2*M1**2)**((-gamma)/(gamma-1))
b1=mdot/(rho1*V1a*2*pi*Rm_ax)

#quantities at station 2 (after rotor)
V2t=L_eul_ax/Um_ax+V1t
Tt2=L_eul_ax/cp+Tt1
W2_mag=sqrt(W1_mag**2-xi_ax*L_eul_ax*2)
W2t=V2t-Um_ax
W2a=sqrt(W2_mag**2-W2t**2)
W2=[W2a, W2t]

V2a=W2a
err=20
i=0

while abs(err) > 10**(-4):
    V2_mag=sqrt(V2a**2+V2t**2)
    T2=Tt2-V2_mag**2/(2*cp)
    T2is=T1+efficiency_TT_ax*(T2-T1)
    p2=(T2is/Tt1)**(gamma/(gamma-1))*Pt1
    rho2=p2/(R*T2)
    V2a_new=mdot/(rho2*2*pi*b1*Rm_ax)
    err=abs(V2a_new-V2a)
    V2a=V2a_new
    M2 = V2_mag / sqrt(gamma * R * T2)
    Pt2 = p2 * (1 + (gamma-1)/2 * M2**2)**(gamma/(gamma-1))
    i=i+1


W2a=V2a
W2=[W2a, W2t]
W2_mag_new=sqrt(W2a**2+W2t**2)
beta2=(W2t/W2a)
V2=[V2a, V2t]
xi_new=(V1a**2-V2a**2+V1t**2-V2t**2+2*Um_ax*(V2t-V1t))/(2*L_eul_ax)
xi=(W1_mag**2-W2_mag_new**2)/(2*L_eul_ax) 


# Mean line design for the stator

eta_S = 0.92
alpha_2 = arctan(V2t/V2a)
alpha_3 = 10 * pi/180# Design choice

print(cos(alpha_2) / cos(alpha_3), "> 0.72 ?") 

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
    V3a_new=mdot/(rho3*2*pi*b1*Rm_ax)
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
rpm = RPM_ax