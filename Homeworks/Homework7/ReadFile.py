
# coding: utf-8

# In[3]:


# Homework 2
# Jan 30 2020
# G Besla

# Load Modules
import numpy as np
import astropy.units as u


# In[4]:


# Define a function that will read in a file
# Input:  filename, e.g. "MW_000.txt"
# Returns:  time (in Myr), total number of particles 
#           and an array with data
# USAGE :   time, total, data = Read("filename")
def Read(filename):
    
    
    # open the file 
    file = open(filename,'r')
    
    #read header info line by line (line will be a string)
    # read first two lines FIRST and store as variable
    
    # read and store time
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr

    # read and store total number of particles
    line2 = file.readline()
    label, value = line2.split()
    total = float(value)
    
    # close file
    file.close()

    # read the remainder of the file, 
    # "dtype=None" specifies data type. None is default float
    # default delimiter is line is split using white spaces
    # "skip_header=3"  skipping the first 3 lines 
    # the flag "names=True" creates arrays to store the date
    #       with the column headers given in line 4 like "m", "x"
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    # Note, another option is loadtxt, skipping the first 3 rows.  
    # data = np.loadtxt(filename,skiprows=3)
    # But this loses the information in the headers
    
    # this will return the time of the snapshot, 
    #total number of particles 
    #and an array that stores the remainder of the data. 
    return time, total, data


# Checking to see if the code works:

# In[5]:


time, total, data =Read("MW_000.txt")


# In[6]:


time


# In[7]:


total


# In[8]:


#Mass of first particle, if you used loadgenfromtxt
data['m'][0]*u.Msun*1e10


# In[21]:


# if you used loadtxt, first store the mass of all particles in a new array
# i.e. store the 2nd column
#Mass = data[:,1]
# Print the Mass of the first particle
#Mass[0]*u.Msun*1e10


# In[9]:


# Type of first particle
data['type'][0]

