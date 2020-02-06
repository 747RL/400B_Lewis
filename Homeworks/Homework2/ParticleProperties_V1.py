#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 09:47:01 2020

@author: ryan
"""

import numpy as np
import astropy.units as u

def ParticleInfo(filename, ptype, pnum):
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    rx = data['x']
    ry = data['y']
    rz = data['z']
    dmag = np.sqrt(rx**2 + ry**2 + rz**2)*u.kpc
    
    vx = data['vx']
    vy = data['vy']
    vz = data['vz']
    vmag = np.sqrt(vx**2 + vy**2 + vz**2)*u.km/u.second
    
    M = data['m']*u.solMass
    
    #index = np.where(data[_]?)
    #_new = data[_][index]
    
    return dmag, vmag, M