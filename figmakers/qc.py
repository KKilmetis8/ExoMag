#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 19:47:07 2024

@author: konstantinos

Figure 1
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.Utilities.profile_sorter import profile_sorter
from src.Bfield.reynolds import dynamo_region
#%%


def convective_flux(cp, T, rho, vel_conv, P, delta):
    return 2 * cp * T * rho**2 * vel_conv**3 / (P * -delta)


def hardq(p, Rdyns):
    r = np.power(10, p.logR)
    start = np.argmin(np.abs(r - Rdyns[-1])) # inner surface
    end = np.argmin(np.abs(r - Rdyns[0])) # outer surface
    # Convective Flux
    T = np.power(10, p.logT)[end:start] # surface -> center
    P = np.power(10, p.logP)[end:start]
    rho = np.power(10, p.logRho)[end:start]
    delta = p.dlnRho_dlnT_const_Pgas[end:start]
    conv_vel = p.conv_vel[end:start]
    cp = p.cp[end:start]
    q = convective_flux(cp, T, rho, conv_vel, P, delta)
    q *= 1e-3 # mW/m^2 to W/m^2

    return q

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.q = []


    def __call__(self, age_h, q_h,):
        self.age.append(age_h)
        self.q.append(q_h)


    
def plotter(names, cols, labels, title, bigfirst = False):
    

    # Makes the calculations
    planets = hardB_doer(names) 
    fig, axs = plt.subplots(1,3, tight_layout = True, sharex = True,
                           figsize = (8,4))
    custom_lines = []
    for planet, color, i in zip(planets, colors, range(len(planets))):
        if i == 0 and bigfirst:
            #axs[0].plot(planet.age, planet.F, color = color, lw=10, zorder = 2)
            axs[0].plot(planet.age, planet.Bdyn, color = color, lw=10, zorder = 2)
            axs[1].plot(planet.age, planet.Bdip, color = color, lw=10, zorder = 2)
            
        if i == 1 and bigfirst:
            axs[2].plot(planet.age, planet.F, color = color, lw=5, zorder = 2)
            axs[0].plot(planet.age, planet.Bdyn, color = color, lw=5, zorder = 2)
            axs[1].plot(planet.age, planet.Bdip, color = color, lw=5, zorder = 2)
        
        axs[2].plot(planet.age, planet.F, color = color)
        axs[0].plot(planet.age, planet.Bdyn, color = color)
        axs[1].plot(planet.age, planet.Bdip, color = color)
    
        # Legend
        custom_lines.append( Line2D([0], [0], color = colors[i], 
                                    label = labels[i]))
        
    # Make nice
    #axs[0].set_ylabel('Efficiency Factor', fontsize = 14)
    axs[0].set_ylabel('Dynamo [G]', fontsize = 14)
    axs[1].set_ylabel('Dipole [G]', fontsize = 14)
    axs[2].set_ylabel('Mean density', fontsize = 14)
    
    #axs[0].grid()
    axs[0].grid()
    axs[1].grid()
        
    axs[1].set_xlim(500, 10_000)
    #axs[0].set_yscale('log')
    if 'm317' in names[0]:
        #axs[0].set_ylim(1_000, 6_000)
        axs[0].set_ylim(0, 400)
        axs[1].set_ylim(0, 250)
        #axs[2].set_ylim(1e2, 1e5)
        #axs[2].set_yscale('log')
    elif 'm17' in names[0]:
        # axs[0].set_ylim(1E-12,1E1)
        axs[0].set_ylim(0, 40)
        axs[1].set_ylim(0, 40)
    
    #axs[0].set_yscale('log')
    #axs[1].set_yscale('log')

    # fig.suptitle(title, 
    #               fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)

    if len(planets) == 2:
        box_x = 0.74
    elif len(planets) == 3:
        box_x = 0.87
    elif len(planets) == 6:
        box_x = 0.83
    elif len(planets) == 8:
        box_x = 0.73
    else:
        box_x = 0.97
        
    fig.legend(custom_lines, labels,
            fontsize =  9, ncols = len(planets), alignment = 'left', # Lawful Neutral
            bbox_to_anchor=(0.94, -0.03), bbox_transform = fig.transFigure,)
#%%    
plt.figure(figsize = (5,5), dpi = 300)

name = 'm317_env0.92_zero_a0.1_s8'
numbers = [17, 37, 57, 77, 97, 112]
colors = [c.c99, c.c97, c.c95, c.AEK, c.c92, c.c91]
texty = [15, 9, 5.5, 4, 3, 2] 
textx = np.array([0.1, 0.095, 0.09, 0.085, 0.08, 0.07]) + 0.9
for num, col, y, x in zip(numbers, colors, texty, textx):
    profile = f'data/{name}/profile{num}.data'
    p = mr.MesaData(profile)
    r = np.power(10, p.logR)
    R_dynamo_active, rmn, age = dynamo_region(p)
    q = hardq(p, R_dynamo_active)
    plt.plot(R_dynamo_active[:-1]/r[0], q, c = col)
    plt.text(x, y, f'{age*1e-3:.1f} Gyr', fontsize = 14, color = col) 
plt.xlabel('Radial Coordinate [$R_\mathrm{p}$]', fontsize = 14)
plt.ylabel('Convective Energy Flux $q_\mathrm{c}$ [erg s$^{-1}$ cm$^{-2}$]', fontsize = 14)
plt.xlim(0.1, 1.2)
plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
plt.axvline(1, c='k', alpha = 0.1, lw = 3)
plt.yscale('log')
plt.savefig('figs/qc.pdf',format='pdf', dpi = 300)

