#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 17:19:57 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [6 , 5]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
AEK = '#F1C410'
import mesa_reader as mr
# Constants
Rsol = 6.957e10 # [cm]
Rjup = 71.492e8 # [cm]
Rearth = 6.371e8 # [cm]
Msol = 1.98855e33 # [g]
Mjup = 1.8986e30 # [g]
Mjup_in_sol = Mjup / Msol
Mearth = 5.9722e27 # [g]
c = 2.99792458e10 #[cm/s]
h = 6.62607015e-27 #[gcm^2/s]
Kb = 1.380649e-16 #[gcm^2/s^2K]
Lsol = 3.826e33 # [erg/s]
sigma = 5.670e-5 # [cgs]
sqrt2 = np.sqrt(2)
#%%
mass = '7'
atm = '0.3'
p_path = 'data/' + mass + '-' + atm + '/profile'
h_path = 'data/' + mass + '-' + atm + '/history_7'
fig, axs = plt.subplots(1,2, tight_layout = True, sharex=True, sharey=True,
                       figsize = (10,8))
profiles = np.arange(1, 2+1, step = 1)
profile_interval = 1
for p, ax in zip(profiles, axs.reshape(-1)):
    
    # Load data
    # h = mr.MesaData(h_path + '.data')
    p = mr.MesaData(p_path + str(p) + '.data')
    r = np.power(10, p.logR) * Rsol / Rearth
    Temperature = np.power(10, p.logT) # [K]
    Density = np.power(10, p.logRho) # [g/cm3]
    
    # Check if it cools everywhere
    monotone = np.diff(Temperature) > 0
    colors = np.where(monotone, 'k', 'r')   

    # Plot
    # ax2 = ax.twinx()
    ax.plot(r[1:], Density[1:], c = 'deepskyblue', zorder = 2)
    # ax2.scatter(r[1:], Temperature[1:], c=colors, s = 10, zorder = 3)
    
    # Age text
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.75)
    ax.text(0.55, 0.9, 'Age: ' + str(np.round(p.star_age,2)) + ' yr', 
            transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    ax.grid()
    # ax.set_yscale('log')
    ax.set_xscale('log')
    
fig.suptitle(mass + r' $M_\oplus$ - ' + str(float(atm) * 100) + '$\%$ $f_{at,0}$', 
             fontsize = 18, y = 0.98)
fig.text(0.5, -0.02, r'Radius $[R_\oplus]$', 
         fontsize = 15, transform = fig.transFigure)
# fig.text(0.99, 0.4, 'Temperature [K]', rotation = 270, 
#           fontsize = 15, transform = fig.transFigure)
fig.text(-0.02, 0.4, 'Density [g/cm$^3$]', rotation = 90, 
          fontsize = 15, transform = fig.transFigure)

# Legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='Deepskyblue', lw = 3),
                Line2D([0], [0], color = 'k', lw = 3),
                Line2D([0], [0], color= 'r', lw = 3)]

fig.legend(custom_lines, ['Density', 'Decreasing Temperature', 'Increasing Tempreature'],
           fontsize =  14, ncols = 3, alignment = 'center', # Lawful Neutral
           bbox_to_anchor=(0.9, -0.03), bbox_transform = fig.transFigure,)