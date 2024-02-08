#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 13:52:25 2024

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
Ts = np.zeros((len(rs), len(profiles)))
Ls = np.zeros((len(rs), len(profiles)))
Ls2 = np.zeros((len(rs), len(profiles)))
for profile, i in zip(profiles, range(len(profiles))):
    p = mr.MesaData(profile)
    r, reynolds_mag_number, age = rey_mag(p)
    ages.append(age)
    for j in range(len(Ts) - 2):
        idx = np.argmin( np.abs(r - rs[j] * r[0] )) # Find closest r
        Ts[j][i] = np.power(10, p.logT[idx]) 
        Ls[j][i] = p.luminosity[idx]
        Ls2[j][i] = 4 * np.pi * (rs[j] * r[0] * c.Rsol)**2 * c.sigma * Ts[j][i]**4 / c.Lsol
    # Dynamo surface 
    R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
    Rdyn = R_dynamo_active[0] # it's the wrong way round
    idx = np.argmin( np.abs(r - Rdyn )) # Find closest r
    Ts[-2][i] = np.power(10, p.logT[idx]) 
    Ls[-2][i] = p.luminosity[idx]
    Ls2[-2][i] = 4 * np.pi * (Rdyn * c.Rsol)**2 * c.sigma * Ts[-2][i]**4 / c.Lsol


    # RCB | Radiation convective boundary
    rad = p.gradr
    conv = p.grada
    # plt.figure()
    # plt.plot(r, rad)
    # plt.plot(r, conv)
    cmorethanb = r[conv < rad]
    RCB = cmorethanb[0]
    idx2 = np.argmin( np.abs(r - RCB )) # Find closest r
    Ts[-1][i] = np.power(10, p.logT[idx2])
    Ls[-1][i] = p.luminosity[idx]
    Ls2[-1][i] = 4 * np.pi * (RCB * c.Rsol)**2 * c.sigma * Ts[-1][i]**4 / c.Lsol

    
colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.reddish, 'k', 'r']
fig, ax =  plt.subplots(1,3, tight_layout = True, sharex = True,
                           figsize = (8,4))
for i in range(len(Ts)):
    if colors[i] == 'r':
        linestyle = 'dotted'
    elif colors[i] == 'k':
        linestyle = '-.'
    else:
        linestyle = 'solid'
    ax[0].plot(ages, Ts[i], c=colors[i], linestyle = linestyle,
               label = str(rs[i]),)
    ax[1].plot(ages, Ls[i], c=colors[i], linestyle = linestyle)
    ax[2].plot(ages, Ls2[i], c=colors[i], linestyle = linestyle)


fig.legend(fontsize =  12, ncols = 4, alignment = 'center', # Lawful Neutral
                bbox_to_anchor=(0.8, 0.03), bbox_transform = fig.transFigure,)

ax[0].set_xlim(100, 7e3)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[2].set_yscale('log')

ax[0].set_xlabel('Age [Myr]')
ax[1].set_xlabel('Age [Myr]')
ax[2].set_xlabel('Age [Myr]')

ax[0].set_ylabel('Temperature [K]')
ax[1].set_ylabel(r'$\texttt{MESA}$ Luminosity [$L_\odot$]')
ax[2].set_ylabel(r'BB Luminosity [$L_\odot$]')

fig.suptitle('Luminosity \& Temperature Evolution, Jupiter, 50\% Env')