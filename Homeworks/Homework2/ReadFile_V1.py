#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:31:30 2020

@author: ryan
"""

import numpy as np
import astropy.units as u

def Read(filename):             # Function that reads in data from filename
    file = open(filename, 'r')
    
    line1 = file.readline()     # Reads in first line
    label, value = line1.split()
    time = float(value)*u.Myr
    
    line2 = file.readline()     # Reads in second line
    label, value = line2.split()
    pnum = float(value)
    
    file.close()                # Closes file
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #stores the remiander of the file
    
    return time, pnum, data     # returns the time, particle number and data array