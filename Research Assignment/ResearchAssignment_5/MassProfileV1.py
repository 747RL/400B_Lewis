#!/usr/bin/env python
# coding: utf-8

# In[46]:


# Homework 5 for ASTR 400B
# Determine a galaxy's rotation curve from the mass distribution at SnapNumber 0
# Ryan Lewis

# In collaboration with Cassie Bodin


# In[188]:


# Import modules
import numpy as np
import astropy.units as u
import astropy
import matplotlib.pyplot as plt

from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)


# In[189]:


# Class to determine a galaxy's rotation curve from the mass distribution at SnapNumber 0
class MassProfile:
    
    # Function that initializes the instances of the class MassProfile
    def __init__(self, Galaxy, Snap):
        # Input:
            # Galaxy - the name of the galaxy, e.g. "MW", "M31", or "M33"
            # Snap - the snapshot number, e.g. 0, 1, etc.
        # Outputs:
            # None - (Function is used for initialization purposes only)
        
        ilbl = '000' + str(Snap) # Add a string of the filenumber to the value "000"
        ilbl = ilbl[-3:] # remove all but the last 3 digits
        self.filename = "%s_"%(Galaxy)+ilbl+'.txt'
        
        self.time, self.total, self.data = Read(self.filename) # Read data in the given file using Read

        # Store the mass and positions of only the particles of the given type
        self.m = self.data['m']       # Mass
        self.x = self.data['x']*u.kpc # x-position
        self.y = self.data['y']*u.kpc # y-position
        self.z = self.data['z']*u.kpc # z-position
        
        self.gname = Galaxy # Store the name of the galaxy as a global property
        self.Snap = Snap
        
    # Function that computes the mass (Msun) enclosed within a given radius of the Center of Mass position 
    # for a specified galaxy and specified compnent of that galaxy
    def MassEnclosed(self, ptype, Renc):
        # Input:
            # ptype - the particle type (Halo: 1, Disk: 2, Bulge: 3)
            # Renc - the galaxy radius within which the mass is enclosed
        # Output:
            # Menc - the mass enclosed within the given radius for the specified galaxy
            
        index = np.where(self.data['type'] == ptype) # Create an array to store indexes of particles of desired Ptype
        
        COM = CenterOfMass(self.filename, ptype) # Create a center of mass object
        COMP = COM.COM_P(0.1) # Calculate the center of mass position
            
        # Change the frame of reference to the newly computed COM position
        xNew = self.x[index] - COMP[0]
        yNew = self.y[index] - COMP[1]
        zNew = self.z[index] - COMP[2]
            
        rCOM = np.sqrt(xNew**2 + yNew**2 + zNew**2) # Calculate the distance from the COM position
            
        Menc = np.zeros(len(Renc)) # Initialize an array of 0's to store the total mass within a given radius
            
        for i in range(len(Renc)):
            rad = Renc[i]
            pindex = np.where(rCOM.value <= rad)
            Mpart = self.m[pindex]
            Menc[i] = np.sum(Mpart)
                
        if self.gname == 'M33' and ptype == 1: # To account for M33 having no bulge
            Menc = np.zeros(len(Renc))
            
        return Menc*1e10*u.Msun
        
    # Function that calculates an array of masses (Msun) representing the total enclosed mass (bulge+disk+halo)
    # at each radius of the given array of radii
    def MassEnclosedTotal(self, Renc):
        # Input:
            # Renc - the galaxy radius within which the mass is enclosed
        # Output:
            # MencTot - the total mass enclosed within the given radius for the specified galaxy
            
        mHalo = self.MassEnclosed(1,Renc)
        mDisk = self.MassEnclosed(2,Renc)
        
        if self.gname != "M33":
            mBulge = self.MassEnclosed(3,Renc)
        else:
            mBulge = 0
            
        MencTot = mHalo + mDisk + mBulge
        
        return MencTot
    
    # Function that computes the mass enclosed (Msun) within a given radius using the the Hernquist 1990 dark halo mass 
    # profile
    # M(r) = Mhalo*r**2/(a+r)**2
    def HernquistMass(self, r, a, mHalo):
        # Input:
            # r - the distance from the cetner of the galaxy (kpc)
            # a - the scale radius (kpc)
            # mHalo - the total dark matter halo mass (10^12 Msun)
        # Output:
            # mHernquist - the total dark matter mass enclosed within r (Msun)
        
        mHernquist = mHalo*r**2/(a+r)**2*u.Msun
        
        return mHernquist
    
    # Function that calculates the circular speed of a given particle type at a given radius
    def CircularVelocity(self, ptype, Renc):
        # Input:
            # ptype - the particle type (Halo: 1, Disk: 2, Bulge: 3)
            # Renc - the galaxy radius within which the mass is enclosed
        # Output:
            # Vcirc - the circular velocity at a given radius for the specified galaxy (km/s)
        
        Vcirc = np.zeros(len(Renc)) # Initialize an array of 0's to store the circular velocity of a particle at a given radius
        
        Menc = self.MassEnclosed(ptype, Renc)
        
        Vcirc = np.around(np.sqrt(G*Menc/Renc*u.kpc), 2)
        
        return Vcirc
        
    # Function that calculates the total circular speed at a given radius
    def CircularVelocityTotal(self, Renc):
        # Input:
            # Renc - the galaxy radius within which the mass is enclosed
        # Output:
            # VcircTot - the circular velocity 
    
        Menc = self.MassEnclosedTotal(Renc)
        
        VcircTot = np.around(np.sqrt(G*Menc/Renc*u.kpc), 2)
        
        return VcircTot
        
    # Function that calculates the cirvular velocity using the Hernquist 1990 dark halo mass profile
    def HernquistVCirc(self, r, a, mHalo):
        # Input:
            # r - the distance from the cetner of the galaxy (kpc)
            # a - the scale radius (kpc)
            # mHalo - the total dark matter halo mass (10^12 Msun)
        # Output:
            #VcircHernquist - the circular velocity (km/s)
            
        Menc = self.HernquistMass(r, a, mHalo)
        
        VcircHernquist = np.around(np.sqrt(G*Menc/r*u.kpc))
        
        return VcircHernquist
    
    # Research Assignment Additions--------------------------------------------------------------------
    
    def NFWMass(self, r, a, mHalo):
        # Input:
            # r - the distance from the cetner of the galaxy (kpc)
            # a - the scale radius (kpc)
            # mHalo - the total halo mass (10^12 Msun)
        # Output:
            # mNFW - the total dark matter mass enclosed within r (Msun)
            
        mNFW = np.round(mHalo/((r/a)*(1+(r/a))**2), 2)*1e12*u.Msun
        
        return mNFW