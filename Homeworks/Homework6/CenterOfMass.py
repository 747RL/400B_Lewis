
# coding: utf-8

# In[1]:


# Homework 4
# Center of Mass Position and Velocity
# Ryan Lewis


# ### Keep in mind this is just a template; you don't need to follow every step and feel free to change anything.
# ### We also strongly encourage you to try to develop your own method to solve the homework.

# In[2]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[3]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        # self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.m = self.data['m'][self.index] #*1e10*u.Msun
        self.x = self.data['x'][self.index] #*u.kpc
        self.y = self.data['y'][self.index] #*u.kpc
        self.z = self.data['z'][self.index] #*u.kpc
        self.vx = self.data['vx'][self.index] #*u.km/u.s
        self.vy = self.data['vy'][self.index] #*u.km/u.s
        self.vz = self.data['vz'][self.index] #*u.km/u.s

    
    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = np.sum(a*m) / np.sum(m)
        # ycomponent Center of mass
        Bcom = np.sum(b*m) / np.sum(m)
        # zcomponent Center of mass
        Ccom = np.sum(c*m) / np.sum(m)
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        
        # compute the magnitude of the COM position vector.
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW<RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            # print ("maxR", RMAX)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units usint astropy and round all values
        # and then return the COM positon vector
        
        COMvector = np.around(COMP, 2) #*u.kpc
        return COMvector
        
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x - COMX
        yV = self.y - COMY
        zV = self.z - COMZ
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV<RVMAX)

        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units using astropy
        # round all values
        # write your own code below
        COMV = np.around([VXCOM, VYCOM, VZCOM], 2) #*u.km/u.s

        # return the COM vector                                                                                        
        return COMV


# In[4]:


# Create a Center of mass object for the MW, M31 and M33
MWCOM = CenterOfMass("MW_000.txt", 2)
M31COM = CenterOfMass("M31_000.txt", 2)
M33COM = CenterOfMass("M33_000.txt", 2)


# In[5]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
MW_COMP = MWCOM.COM_P(0.1)
MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])

# M31:   store the position and velocity COM
M31_COMP = M31COM.COM_P(0.1)
M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])

# M33:   store the position and velocity COM
M33_COMP = M33COM.COM_P(0.1)
M33_COMV = M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2])


# In[6]:


# Questions

# 1
# What is the COM position and velocity vector for the MW, M31 and M33 at Snapshot 0 
# using Disk Particles only (use 0.1 kpc as the tolerance so we can have the same answers 
# to compare)? In practice, disk particles work the best for the COM determination. 
# Recall that the MW COM should be close to the origin of the coordinate system (0,0,0).
print("Position vectors: ")
print("The position vector for the MW is ", MW_COMP *u.kpc)
print("The position vector for the M31 is ", M31_COMP *u.kpc)
print("The position vector for the M33 is ", M33_COMP *u.kpc)
print("\n")
print("Velocity vectors: ")
print("The velocity vector for the MW is ", MW_COMV *u.km/u.s)
print("The velocity vector for the M31 is ", M31_COMV *u.km/u.s)
print("The velocity vector for the M33 is ", M33_COMV *u.km/u.s)


# In[7]:


# 2
# What is the magnitude of the current separation and velocity between the MW and M31? 
# From class, you already know what the relative separation and velocity should roughly be
# (Lecture2 Handouts; Jan 16).

SepMW_M31 = MW_COMP - M31_COMP
MagSepMW_M31 = np.linalg.norm(SepMW_M31)
print("The magnitude of the separation between MW and M31 is ", np.around(MagSepMW_M31, 2))

VelMW_M31 = MW_COMV - M31_COMV
MagVelMW_M31 = np.linalg.norm(VelMW_M31)
print("The magnitude of the velocity between MW and M31 is ", np.around(MagVelMW_M31, 2))


# In[8]:


# 3
# What is the magnitude of the current separation and velocity between M33 and M31?

SepM33_M31 = M33_COMP - M31_COMP
MagSepM33_M31 = np.linalg.norm(SepM33_M31)
print("The magnitude of the separation between MW and M31 is ", np.around(MagSepM33_M31, 2))

VelM33_M31 = M33_COMV - M31_COMV
MagVelM33_M31 = np.linalg.norm(VelM33_M31)
print("The magnitude of the velocity between M33 and M31 is ", np.around(MagVelM33_M31, 2))


# In[9]:


# 4
# Given that M31 and the MW are about to merge, why is the iterative process to determine 
# the COM is important?

print(" The iterative process to determine the COM is important because as the galaxies merge, we can determine the COM of all particles within a given radius.")

