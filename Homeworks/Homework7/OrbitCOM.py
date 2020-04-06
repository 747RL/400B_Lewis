
# coding: utf-8

# In[ ]:


# Homework 6 Solution
# G. Besla & R. Li


# In[3]:


# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass


# for bonus:
# scipy provides powerful functions for scientific computing
#from scipy.signal import find_peaks
#from scipy.optimize import curve_fit


# In[4]:


def OrbitCOM(galaxy,start,end,n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time
    inputs:
           galaxy  the name of the galaxy e.g. "MW"
           start   initial SnapNumber  e.g. 0
           end     final SnapNumber  e.g. 100
           n       integer indicating the frequency with which the snapshots are read in;  n should not be 0
    returns:  file with COM pos/vel of a galaxy at snapshots over the range (start, end+n, n) 
                     columns t, x, y, z, vx, vy, vz for n snapshots
    """
    
    # compose the filename for output
    #fileout = "Orbit_" + galaxy + ".txt"
    fileout = 'Orbit_%s.txt'%(galaxy)   # filename to store output

    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    delta, VolDec = 0.1, 2.0
    # for M33 that is stripped more, use different values for VolDec
    if galaxy is "M33": 
        delta, VolDec = 0.1, 4.0
    
    # generate the snapshot id sequence and check if the input is eligible
    snap_ids = np.arange(start, end+n, n)
    if snap_ids.size == 0:
        raise ValueError("Cannot build a sequence using the input, start = ", 
                         start, "end = ", end, "n = ", n)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([snap_ids.size, 7])
    
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        filename = "VLowRes/" + galaxy + "_VLowRes/" + galaxy +  "_{:03d}".format(snap_id) + ".txt"
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_pos = COM.COM_P(delta, VolDec)
        COM_vel = COM.COM_V(COM_pos[0],COM_pos[1], COM_pos[2])
        
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i] = COM.time.value, *tuple(COM_pos.value), *tuple(COM_vel.value)
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


# In[6]:


# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
print("Starting to Compute MW Orbit")
OrbitCOM("MW",0,800,5)
print("Starting to Compute M31 Orbit")
OrbitCOM("M31",0,800,5)
print("Starting to Compute M33 Orbit")
OrbitCOM("M33",0,800,5)


# In[9]:


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

MWOrbit = np.genfromtxt('Orbit_MW.txt',dtype=None,names=True) 
M31Orbit = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True) 
M33Orbit = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True) 


# In[17]:


# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def relative_mag(orbit1, orbit2): 
    # Input: 3D vector 1 , 3D vector 2
    # returns the magnitude of the relative position and velocities
    pos =  np.sqrt( (orbit1['x']-orbit2['x'])**2 + (orbit1['y']-orbit2['y'])**2 + (orbit1['z']-orbit2['z'])**2)
    vel =  np.sqrt( (orbit1['vx']-orbit2['vx'])**2 + (orbit1['vy']-orbit2['vy'])**2 + (orbit1['vz']-orbit2['vz'])**2)
    return pos,vel


# In[18]:


# Determine the magnitude of the relative position and velocities 

# of MW and M31
M31_MW_R, M31_MW_V = relative_mag(MWOrbit, M31Orbit)

# of M33 and M31
M33_M31_R, M33_M31_V = relative_mag(M31Orbit, M33Orbit)


# In[37]:


# Plot the Orbit of the galaxies 
#################################


fig, ax= plt.subplots(figsize=(10, 8))


# Plot the separtion of M31 and MW
ax.plot(MWOrbit['t']/1e3, M31_MW_R, color='blue', linewidth=5, label='MW-M31')

# Plot the separtion of M33 and M31
ax.plot(M31Orbit['t']/1e3, M33_M31_R, color='red', linewidth=5, label='M31-M33')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)

# Set the title
ax.set_title("Separation vs. Time", fontsize=22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

# save 
fig.savefig('orbit_sep.pdf')

plt.show()



# Plot the orbital velocities of the galaxies 
#################################

fig, ax= plt.subplots(figsize=(10, 8))

# Plot the velocities of M31 and MW
ax.plot(MWOrbit['t']/1e3, M31_MW_V, 'b', lw=5, label='MW-M31')

# Plot the velocities of M33 and M31
ax.plot(MWOrbit['t']/1e3, M33_M31_V, 'r', lw=5, label='M31-M33')

# Add axis labels
ax.set_xlabel('Time (Gyr)', fontsize=22)
ax.set_ylabel('Velocity (km/s)', fontsize=22)
ax.set_title("Velocity vs. Time", fontsize=22)

#set axis limits
ax.set_ylim([1, 600])

#adjust tick label font size
ax.xaxis.set_tick_params(labelsize=22)
ax.yaxis.set_tick_params(labelsize=22)

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize=22)

# tight layout
fig.tight_layout()

# Save to a file
fig.savefig('orbit_vel.pdf')

plt.show()


# In[36]:


# Plot a zoom of the MW-M31 Orbit, using log y axis 
#################################


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot the separtion of M31 and MW
plt.semilogy(MWOrbit['t']/1000, M31_MW_R, color='blue', linewidth=5, label='MW-M31')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)

#set axis limits
plt.ylim(0,50)
plt.xlim(6,6.8)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


# Q1: 2 encounters and then they merge on the 3rd encounter
# 
# Q2:  The velocity peaks when the separation is smallest and is at a minimum when the separation is largest.  Energy is being exchanged - from potential energy to kinetic at pericenter. The opposite occurs at apocenter.
# 
# Q3:  a) MW-M31 Merger will occur in ~ 6.5 Gyr 
# Note that the force resolution (gravitational softening) of the N-body is of order 3 kpc.
# Separations a few times that scale are not going to be well resolved. 
# 
# Q3: b)  M33's orbit is slowly decaying over time, but the orbit does change a bit at the time of the merger
# and starts to decay more quickly. This is because there is more material in the central regions of the 
# combined MW+M31 object than about M31 alone. 

# In[43]:


# BONUS: 

# get orbits after 6 Gyr 
decay_time = MWOrbit['t'][MWOrbit['t']>6e3]/1e3
decay_orbit = M33_M31_R[MWOrbit['t']>6e3]

# find those four local maxima in the figure
indices, _ = find_peaks(decay_orbit)

# get four apocenters
apo_time = decay_time[indices]
apo_sep = decay_orbit[indices]

# print decay rate as defined in the HW instructions
print("decay rate over three orbits:", 
      (apo_sep[1:] - apo_sep[:-1]) / (apo_time[1:] - apo_time[:-1]), 
      'kpc / Gyr')

# now fit an exponential decay function
exp_decay = lambda t, a, b, c : a * np.exp(-b * t) + c
# using curve_fit
popt, pcov = curve_fit(exp_decay, apo_time, apo_sep)

plt.show()

# In[52]:


# Fit the orbital decay of M33-M31
#################################

fig, ax= plt.subplots(figsize=(10, 8));

# Plot the velocities of M33 and M31
ax.plot(MWOrbit['t']/1e3, M33_M31_R, 'b-', lw=5, label='Separations')

# Plot the fitted curve
ax.plot(np.linspace(6, 20, 100), exp_decay(np.linspace(6, 20, 100), *popt), 
        'r-', lw=5, label='Fitted Exp Decay')

# Add axis labels
ax.set_xlabel('Time (Gyr)', fontsize=22)
ax.set_ylabel('Separation (kpc)', fontsize=22)
ax.set_title("M31-M33", fontsize=22)
ax.set_xlim([6, 20]); ax.set_ylim([25, 125])

#adjust tick label font size
ax.xaxis.set_tick_params(labelsize=22)
ax.yaxis.set_tick_params(labelsize=22)

# add a legend with some customizations.
ax.legend(loc='best', fontsize=22)

# tight layout
fig.tight_layout()

# Save to a file
fig.savefig('orbital_decay.pdf')

plt.show()


# Q4 Bonus: 
# 
# The orbital period decays from ~110 to 90 over 1.5 Gyr then from 90 to 75 over 1.2 Gyr 
# if this rate of decay were to continue, you would expect it to take 75 /[(110-90)/1.5] ~ 5.6 Gyr for the orbit of M33 to decay from 75 kpc. 
# 
# BUT the rate of decay is slowing down by the end of the simulation - rate of (75-70)kpc/1 Gyr ~ 4 kpc/Gyr. At that rate it would take more than 20 Gyr to merge!   Galaxies merge because their orbits decay owing to dynamical friction. This friction term scales with the mass of the galaxy. As M33's mass is stripped by tidal effects, the strength of dynamical friction weakens and the decay rate of the orbit slows.  It will take a very long time for the M33 remnant to merge! 
# 
# This fits well with the expectations if you fit an exponential decay rate - like Andy Henrici did. See an adaptation of his code in the plot below.
# 
# NOTE: If you don't decrease RMAX by 4 for M33 the Orbit will look pretty ratty because there is a lot of tidal debris
# Another, more clever way, would be to use a smaller volume centered on the COM from the previous snapshot as the starting point. 