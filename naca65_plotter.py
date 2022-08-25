# Naca 65 series plotter function
import numpy as np
import matplotlib.pyplot as plt

def naca65(theta, maxTh, chord, origin, angle_rot):
    # theta = equivalent camber angle
    # maxTh = Profile maximum thickness
    # origin = Origin point for the profile rotation
        # If origin == "False" the rotation is done around the center of the area 
        # After the Rotation the center of rotation is moved to the origin (0,0)
    # angle_rot = angle of rotation

    import numpy as np
    import matplotlib.pyplot as plt
    # Starting Point Naca 65-(10)10 
    # C_l = 1.0
    # max thickness / chord = 0.10 

    xc = np.array([0, 0.5, 0.75, 1.25, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]) # Coordinate along chord x/c %
    yc = np.array([0, 0.25, 0.35, 0.535, 0.93, 1.58, 2.12, 2.585, 3.365, 3.98, 4.475, 4.86, 5.15, 5.355, 5.475, 5.515, 5.475, 5.355, 5.15, 4.86, 4.475, 3.98, 3.365, 2.585, 1.58, 0]) # Camber line coordinates y/c %
    th = np.array([0, 1.544, 1.864, 2.338, 3.480, 4.354, 5.294, 6.08 , 7.332, 8.286, 9.006, 9.52 , 9.848, 9.992, 9.926, 9.624, 9.06 , 8.292, 7.364, 6.312, 5.168, 3.974, 2.77 , 1.62 , 0.612, 0]) # Profile thickness t/c %

    # Modified thickness distribution for blunter T.E.
    # i=0
    # for j in range(17,25):
    #     add = np.linspace(0,0.986, 8)
    #     th[j] = th[j]+add[i]
    #     i=i+1

    
    C_l = theta / 25

    profile_name = "NACA 65-(" + str(int(round(abs(10*C_l)))) + ")" + str(int(maxTh * 100))

    # Scaling and correction for chosen theta, thickness and chord
    yc = yc * chord * C_l / 100        # Multiply camber line by C_l to obtain a different lift
    xc = xc * chord / 100
    t =  th * chord * maxTh * 10 / 100 # Multiply by a different max thickness 

    # Upper and Lower surface x and y coordinates
    ux = xc
    uy = yc + t / 2

    lx = xc 
    ly = yc - t / 2

    # Profile Area Calculation
    profile_area = np.trapz(uy,xc) - np.trapz(ly,xc)
    
    xAc = (np.trapz( uy * xc , xc) - np.trapz(ly * xc,xc) )* 1/profile_area  # Center of the Area position in x
    index_Ac = np.where(abs(xc-xAc) < 0.05*xAc ) # Approximate the position of A.C. w/ a point on the camber line

    # Moment calculation
    S_x = - np.trapz(xc * uy ,uy) + np.trapz(xc * ly ,ly) # First order (static) moment
    I_x = - np.trapz(xc * uy**2 ,uy) + np.trapz(xc * ly**2 ,ly) # Moment of inertia
    yAc = S_x / profile_area

    ## Rotation
    xi = angle_rot * np.pi/180  # Rotation angle
    
    if origin == "False" : # Rotate around center of the area 
        ox = xc[index_Ac]
        oy = yc[index_Ac]
    else:
        ox, oy = origin # Origin of the rotation
    
    R_mat= np.array([[np.cos(xi),-np.sin(xi)], 
                    [np.sin(xi), np.cos(xi)]]) # Rotation Matrix

    Xc,Yc,Ux,Uy,Lx,Ly = (np.zeros(len(xc)) for t in range(6)) # Initiate empty arrays

    # Rotation loop, remove ox oy for centering at the origin of rotation
    for j in range(len(xc)):
        Xc[j] =  (xc[j]-ox) * np.cos(xi) - (yc[j]-oy) * np.sin(xi) #+ ox
        Yc[j] =  (xc[j]-ox) * np.sin(xi) + (yc[j]-oy) * np.cos(xi) #+ oy
        Ux[j] =  (ux[j]-ox) * np.cos(xi) - (uy[j]-oy) * np.sin(xi) #+ ox
        Uy[j] =  (ux[j]-ox) * np.sin(xi) + (uy[j]-oy) * np.cos(xi) #+ oy
        Lx[j] =  (lx[j]-ox) * np.cos(xi) - (ly[j]-oy) * np.sin(xi) #+ ox
        Ly[j] =  (lx[j]-ox) * np.sin(xi) + (ly[j]-oy) * np.cos(xi) #+ oy
    
    geom = [profile_area, index_Ac[0][0], I_x]
    return Xc,Yc,Ux,Uy,Lx,Ly, profile_name, geom

# Test

# chord = 0.08
# for angle in [0]:
#     Xc,Yc,Ux,Uy,Lx,Ly, profile_name, geom = naca65(33, 0.1, chord, "False", angle)
#     # Plot
    
#     plt.plot(Xc,Yc,'-.',color='r', linewidth=1)
#     plt.plot(Ux,Uy, 'k')
#     plt.plot(Lx,Ly, 'k')
#     #plt.axis('scaled')
#     plt.axis('square')
#     plt.grid(alpha=0.2)
#     #plt.xlim(-0.1*chord*np.cos(angle * np.pi/180), 1.1*chord*np.cos(angle*np.pi/180))
#     #plt.ylim(-0.1*chord*np.sin(angle * np.pi/180), 1.1*chord*np.sin(angle*np.pi/180))
#     plt.title(profile_name)

# print(geom[2])
# plt.show()
