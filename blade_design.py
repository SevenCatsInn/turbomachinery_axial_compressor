import numpy as np
import matplotlib.pyplot as plt
from naca65_plotter import * 

def lieblein_design(beta_in, beta_out, percent_th, chord, solidity, theta, rr):
    # Design Parameters from computations above
    # Inlet
    mean_index = len(beta_in)//2
    R_m = rr[mean_index]
    R_t = rr[-1]
    R_h = rr[0]
    th = chord * percent_th / 100 # [m] Actual thickness
    
    beta_1 = beta_in
    beta_1root = beta_1[0]          * 180 / np.pi #deg
    beta_1mid  = beta_1[mean_index] * 180 / np.pi #deg
    beta_1tip  = beta_1[-1]         * 180 / np.pi #deg

    # Outlet
    beta_2 = beta_out
    beta_2root = beta_2[0]          * 180 / np.pi #deg
    beta_2mid  = beta_2[mean_index] * 180 / np.pi #deg
    beta_2tip  = beta_2[-1]         * 180 / np.pi #deg

    # Deflection
    deltabeta1_root = beta_1root - beta_2root
    deltabeta1_mid  = beta_1mid  - beta_2mid 
    deltabeta1_tip  = beta_1tip  - beta_2tip 
    
    deltaBeta = [deltabeta1_root, deltabeta1_mid, deltabeta1_tip]

    s_mid = solidity * chord

    n_blade = round( 2 * np.pi * R_m / s_mid) # Number of blades

    # Pitch along blade span
    s_mid  = 2 * np.pi * R_m / n_blade
    s_tip  = 2 * np.pi * R_t / n_blade
    s_root = 2 * np.pi * R_h / n_blade

    # Solidity along blade span, recomputed after choosing n blades
    sigma_mid  = chord / s_mid
    sigma_tip  = chord / s_tip
    sigma_root = chord / s_root

    # Equivalent camber theta: from graphs on slide 9 ppt
    # NOTE: TUNABLE
    theta_eq_root = abs(theta[0])
    theta_eq_mid  = abs(theta[1])
    theta_eq_tip  = abs(theta[2])

    # Zero camber incidence angle

    # Semiempirical exponents
    p_mid  = 0.914 + sigma_mid **3 / 160
    p_tip  = 0.914 + sigma_tip **3 / 160
    p_root = 0.914 + sigma_root**3 / 160

    # Incidence angles in degrees
    i0_10_mid  = abs(beta_1mid)**p_mid   / (5 + 46*np.exp(-2.3 * sigma_mid))  - 0.1 * sigma_mid**3  * np.exp((abs(beta_1mid)-70)/4)
    i0_10_tip  = abs(beta_1tip)**p_tip   / (5 + 46*np.exp(-2.3 * sigma_tip))  - 0.1 * sigma_tip**3  * np.exp((abs(beta_1tip)-70)/4) 
    i0_10_root = abs(beta_1root)**p_root / (5 + 46*np.exp(-2.3 * sigma_root)) - 0.1 * sigma_root**3 * np.exp((abs(beta_1root)-70)/4)


    # Correction factor for profile thickness
    q = 0.28 / (0.1 + (th/chord)**0.3) # Empirical exponent

    k_ith_mid  = (10*th/chord)**q
    k_ith_tip  = (10*th/chord)**q
    k_ith_root = (10*th/chord)**q


    # Slope factor n
    n_mid  = 0.025 * sigma_mid  - 0.06 - ( (abs(beta_1mid)/90)  ** (1+1.2*sigma_mid)  ) / (1.5 + 0.43 * sigma_mid)
    n_tip  = 0.025 * sigma_tip  - 0.06 - ( (abs(beta_1tip)/90)  ** (1+1.2*sigma_tip)  ) / (1.5 + 0.43 * sigma_tip)
    n_root = 0.025 * sigma_root - 0.06 - ( (abs(beta_1root)/90) ** (1+1.2*sigma_root) ) / (1.5 + 0.43 * sigma_root)


    # Optimal Incidence Angle
    i_opt_mid  = i0_10_mid  * k_ith_mid  + n_mid  * theta_eq_mid
    i_opt_tip  = i0_10_tip  * k_ith_tip  + n_tip  * theta_eq_tip
    i_opt_root = i0_10_root * k_ith_root + n_root * theta_eq_root

    inc = [i_opt_root, i_opt_mid, i_opt_tip ]

    # Deviation Angle
    delta0_mid  = 0.01 * sigma_mid  * abs(beta_1mid)  + (0.74*sigma_mid**1.9 + 3 * sigma_mid) * (abs(beta_1mid) / 90) ** (1.67 + 1.09 * sigma_mid) 
    delta0_tip  = 0.01 * sigma_tip  * abs(beta_1tip)  + (0.74*sigma_tip**1.9 + 3 * sigma_tip) * (abs(beta_1tip) / 90) ** (1.67 + 1.09 * sigma_tip)  
    delta0_root = 0.01 * sigma_root * abs(beta_1root) + (0.74*sigma_root**1.9 + 3 * sigma_root) * (abs(beta_1root) / 90) ** (1.67 + 1.09 * sigma_root) 


    # Correction for thickness effects on deviation
    kdeltath_mid  = 6.25 * (th/chord) + 37.5 * (th/chord)**2    
    kdeltath_tip  = 6.25 * (th/chord) + 37.5 * (th/chord)**2   
    kdeltath_root = 6.25 * (th/chord) + 37.5 * (th/chord)**2  


    # Exponent Factor b
    b_mid  = 0.9625 - 0.17 * abs(beta_1mid)/100 - 0.85 * (abs(beta_1mid)/100)**3
    b_tip  = 0.9625 - 0.17 * abs(beta_1tip)/100 - 0.85 * (abs(beta_1tip)/100)**3
    b_root = 0.9625 - 0.17 * abs(beta_1root)/100 - 0.85 * (abs(beta_1root)/100)**3


    # Slope factor m - m1 = Reference value for solidity = 1
    m1_mid  = 0.17 - 0.0333 * abs(beta_1mid) /100 + 0.333 * (beta_1mid/100) **2
    m1_tip  = 0.17 - 0.0333 * abs(beta_1tip) /100 + 0.333 * (beta_1tip/100) **2
    m1_root = 0.17 - 0.0333 * abs(beta_1root)/100 + 0.333 * (beta_1root/100)**2


    m_mid  = m1_mid  / (sigma_mid  **b_mid )
    m_tip  = m1_tip  / (sigma_tip  **b_tip )
    m_root = m1_root / (sigma_root **b_root)


    #computation of deviation
    delta_mid = delta0_mid * kdeltath_mid + m_mid * theta_eq_mid / (sigma_mid ** b_mid)
    delta_tip = delta0_tip * kdeltath_tip + m_tip * theta_eq_tip / (sigma_tip ** b_tip)
    delta_root = delta0_root * kdeltath_root + m_root * theta_eq_root / (sigma_root ** b_root)

    dev = [delta_root, delta_mid, delta_tip]

    #final delta beta
    deltabetafinal_mid = theta_eq_mid - delta_mid + i_opt_mid
    deltabetafinal_tip = theta_eq_tip - delta_tip + i_opt_tip
    deltabetafinal_root = theta_eq_root - delta_root + i_opt_root
    
    # Fix the signs
    Beta = [beta_1root, beta_1mid, beta_1tip]
    
    # Sign corrected quantities
    Theta = [] 
    Inc = []
    Dev = []
    DeltaBeta = []

    j = 0
    for bb in Beta:
        if bb < 0 :
            Theta.append(-theta[j])
            Inc.append(-inc[j])
            Dev.append(-dev[j])
        else:
            Theta.append(theta[j])
            Inc.append(inc[j])
            Dev.append(dev[j])
        DeltaBeta.append(Theta[j] - Dev[j] + Inc[j])
    
        

        j += 1

            
            
    return Inc, Theta, Dev, DeltaBeta







def compute_C_l(theta, pts):
    # Simple function to compute an interpolation of C_l starting from the equivalent camber at 3 sections
    C_l_tmp0 = (list(np.linspace( theta[0]/25, theta[1]/25, pts//2 + 1) ))[0:-1]
    C_l_tmp1 =  list(np.linspace( theta[1]/25, theta[2]/25, pts//2 + 1) )

    C_l = C_l_tmp0 + C_l_tmp1

    return C_l








def printPlot_blade (alpha_0,alpha_1, deltaAlpha0, inc0, theta0, percent_th0, chord0, pts) :
    # Prints the relevant angles and plots the 3 blades at hub, mid and tip
    mean_index = pts // 2

    print("")
    print("Lieblein deflection ROOT = ", deltaAlpha0[0])
    print("Design deflection   ROOT = ", (180/np.pi*(alpha_0[0]-alpha_1[0])))

    print("")
    print("Lieblein deflection MID = ", deltaAlpha0[1])
    print("Design deflection   MID = ", (180/np.pi*(alpha_0[mean_index]-alpha_1[mean_index])))

    print("")
    print("Lieblein deflection TIP = ", deltaAlpha0[2])
    print("Design deflection   TIP = ", (180/np.pi*(alpha_0[-1]-alpha_1[-1])))



    Alpha0 = np.array([alpha_0[0],alpha_0[mean_index],alpha_0[-1]]) * 180/np.pi


    plt.figure()
    
    Geom = []
    Prof_names = []

    for theta, beta, inc, color in zip(theta0, Alpha0, inc0, ['c','b','y']):
        
        stagger =  beta - inc

        Xc,Yc,Ux,Uy,Lx,Ly, profile_name, geom = naca65(theta, percent_th0/100 , chord0, "False", stagger )

        Geom.append(geom)
        Prof_names.append(profile_name)
        
        # plt.plot(Xc,Yc,'-.',color='r', linewidth=1)
        plt.plot(Ux,Uy, color)
        plt.plot(Lx,Ly, color)

        plt.axis('square')
        plt.grid(alpha=0.2)
        # plt.title(profile_name)
    plt.legend(["Hub","","Mean","","","Tip"])

    

    return Geom, Prof_names

