#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 19:15:57 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [6 , 5]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
AEK = '#F1C410'
import mesa_reader as mr

###
# Constants
### 
Rsol = 6.957e10 # [cm]
Rjup = 71.492e8 # [cm]
Rearth = 6.371e8 # [cm]
Msol = 1.98855e33 # [g]
Mjup = 1.8986e30 # [g]
Mjup_in_sol = Mjup / Msol
Mearth = 5.9722e27 # [g]
Bjup = 4170 # [G]
sqrt2 = np.sqrt(2)

path = 'data/12.5-0.2/profile'
p = mr.MesaData(path + '20.data')
Temperature = np.power(10, p.logT) # [K]
Den = np.power(10, p.logRho) # ???
r = np.power(10, p.logR)
Pressure = np.power(10, p.logP)
def Ht(r):
    ''' Temperature Scale Height 
    Ht = cp / (ag)
    
    Cp: Heat Capacity
    a: Thermal Expansivity
    g: Gravity Acceleration
    '''
    
    Ht =1
    return Ht
