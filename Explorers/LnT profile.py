#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:12:18 2024

@author: konstantinos
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.reynolds import rey_mag, profile_sorter

     
name = 'jup17vol3'
lum = False
rs = [0.8, 0.9, 0.97, 0.98, 0.99, 'Rdyn', 'RCB']
colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]

ages = []
p_path = 'data/' + name
profiles = os.popen('ls ' + p_path + '/profile*.data').read()
profiles = list(profiles.split("\n"))
profiles.pop() # Remove last
profiles = profile_sorter(profiles) # Guess what that does

colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.reddish, 'k', 'r']
fig, ax =  plt.subplots(1,2, tight_layout = True, sharex = True,
                           figsize = (8,4))

p1 = mr.MesaData(profiles[5])
r1 = np.power(10, p1.logR)
T1 = np.power(10, p1.logT)
L1_mesa = p1.luminosity
L1_BB = 4 * np.pi * (r1 * c.Rsol)**2 * c.sigma * T1**4 / c.Lsol
age1 = str(np.round(p1.star_age /  1e6, 2)) + ' Myr'


p2 = mr.MesaData(profiles[10])
r2 = np.power(10, p2.logR)
T2 = np.power(10, p2.logT)

L2_mesa = p2.luminosity
L2_BB = 4 * np.pi * (r2 * c.Rsol)**2 * c.sigma * T2**4 / c.Lsol
age2 = str(np.round(p2.star_age /  1e6, 2)) + ' Myr'

ax[0].plot(r1 * c.Rsol/ c.Rearth, T1, c = c.AEK, label = age1)
ax[0].plot(r2 * c.Rsol/ c.Rearth, T2, c = 'k', label = age2)

ax[1].plot(r1 * c.Rsol/ c.Rearth, L1_mesa, c = c.c92, label = age1 + ' - MESA')
ax[1].plot(r2 * c.Rsol/ c.Rearth, L2_mesa, c = c.c97, label = age2 + ' - MESA')

ax[1].plot(r1 * c.Rsol/ c.Rearth, L1_BB, c = c.c92, linestyle = '-.', 
           label = age1 + ' - BB')
ax[1].plot(r2 * c.Rsol/ c.Rearth, L2_BB, c = c.c97,  linestyle = '-.', 
           label = age2 + ' - BB')

ax[0].legend()
ax[1].legend()

ax[0].set_yscale('log')
ax[1].set_yscale('log')

ax[0].set_ylabel('Temperature [K]', fontsize = 15)
ax[1].set_ylabel(r'Luminosity [$L_\oplus$]', fontsize = 15)

ax[0].set_xlabel(r'r $[R_\oplus]$', fontsize = 15)
ax[1].set_xlabel(r'r $[R_\oplus]$', fontsize = 15)
