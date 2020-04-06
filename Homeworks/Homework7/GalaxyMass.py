
# coding: utf-8

# In[1]:

# Homework 3 Solutions
# Computing the Mass of the Local Group
# G. Besla


# In[2]:

# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# astropy provides unit system for astronomical calculations
import astropy.units as u
# import previous HW functions
from ReadFile import Read


# # Galaxy Mass

# In[3]:


def ComponentMass(filename, PType):
# function to compute the mass of all particles of a given type for a given galaxy 
# input:  filename (str) ,  Particle type (1,2,3)
# output: Mass in units of 1e12 Msun. 
    
     # read in the file 
    time, total, data = Read(filename)
    
      # gather particles with the same type and sum up the mass
    # we can directly round the result and adjust the units
    mass = np.round(data[data['type']==PType]['m'].sum() * 1e10/1e12, 3)
    
    # return the mass 
    return mass


# In[5]:

# MW Mass
MW_halo = ComponentMass("MW_000.txt",1)
MW_disk = ComponentMass("MW_000.txt",2)
MW_bulge = ComponentMass("MW_000.txt",3)


# In[13]:

print(f"MW halo Mass: {MW_halo:.3f} x 10^12 Msun") 
print(f"MW disk Mass: {MW_disk:.3f} x 10^12 Msun")
print(f"MW bulge Mass: {MW_bulge:.3f} x 10^12 Msun")


# In[15]:

# Total MW mass
TotalMW = MW_halo + MW_disk + MW_bulge
print(f"Total MW Mass: {TotalMW:.3f} x 10^12 Msun")


# In[16]:

# M31 Mass
M31_halo = ComponentMass("M31_000.txt",1)
M31_disk = ComponentMass("M31_000.txt",2)
M31_bulge = ComponentMass("M31_000.txt",3)


# In[17]:

print(f"M31 halo Mass: {M31_halo:.3f} x 10^12 Msun") 
print(f"M31 disk Mass: {M31_disk:.3f} x 10^12 Msun")
print(f"M31 bulge Mass: {M31_bulge:.3f} x 10^12 Msun")


# In[19]:

# Total M31 Mass
TotalM31 = M31_halo + M31_disk + M31_bulge
print(f"Total M31 Mass: {TotalM31:.3f} x 10^12 Msun")


# In[20]:

# M33 Mass
M33_halo = ComponentMass("M33_000.txt",1)
M33_disk = ComponentMass("M33_000.txt",2)


# In[21]:

print(f"M33 halo Mass: {M33_halo:.3f} x 10^12 Msun") 
print(f"M33 disk Mass: {M33_disk:.3f} x 10^12 Msun")


# In[22]:

# Total M33 Mass
TotalM33 = M33_halo + M33_disk
print(f"Total M33 Mass: {TotalM33:.3f} x 10^12 Msun")


# In[23]:

# Total Local Group Mass
TotalLG  = TotalMW + TotalM31 + TotalM33
print(f"Total Local Group Mass: {TotalLG:.3f} x 10^12 Msun")


# In[ ]:

