#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 20:49:45 2024

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.reynolds import profile_sorter, rey_mag

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.mass = []
        self.radius = []
        self.lum = []
        self.lum2 = [] # comparing MESA logT with photosphereL
      
    def __call__(self, age_h, mass_h, radius_h, lum_h, lum2_h):
        self.age.append(age_h)
        self.mass.append(mass_h)
        self.radius.append(radius_h)
        self.lum.append(lum_h)
        self.lum2.append(lum2_h)
      
def doer(names):
    # Count and generate profile lists
    apothikh = []
    for name in names:
        # History data wrangling
        # h_path = 'data/' + name + '/history_7.data'
        # h = mr.MesaData(h_path)
        # print(dir(h))
        # h_age = np.round(10**h.log_star_age /  1e9, 2) # Gyr
        
        # Profile data wrangling and bookeeping
        hold = apothicarios(name)
        p_path = 'data/' + name
        profiles = os.popen('ls ' + p_path + '/profile*.data').read()
        profiles = list(profiles.split("\n"))
        profiles.pop() # Remove last
        profiles = profile_sorter(profiles) # Guess what that does
        
        for profile, i in zip(profiles, range(len(profiles))):
            # Load data
            p = mr.MesaData(profile)
            r = np.power(10, p.logR)
            age = np.round(p.star_age /  1e6, 2)
            
            
            # Find closest age in history file
            # idx = np.argmin(np.abs(h_age - age))
            M = p.star_mass * c.Msol / c.Mearth
            R = r[0] * c.Rsol / c.Rearth
            # print(R * c.Rsol / c.Rjup)
            L = p.photosphere_L
            #
            r_cgs = r[0] * c.Rsol
            
            # Get Rdyn
            _, reynolds_mag_number, _= rey_mag(p)
            R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
            Rdyn = R_dynamo_active[0] * c.Rsol # [cgs]
            idx = np.argmin(np.abs(r - Rdyn)) 
            
            print(np.power(10, p.logT[idx]), Rdyn)
            L2 = 4 * np.pi * Rdyn**2 * c.sigma* np.power(10, p.logT[idx])**4
            L2 /= c.Lsol
            # Save
            hold(age, M, R, L, L2)
        print('--')
        apothikh.append(hold)
    return apothikh


def plotter(names, cols, labels, bigfirst = True):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'maroon' , 'k', 'maroon']
        colors2 = ['b', 'r', 'g'] 

    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
        
    elif cols == 3:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c98, c.c99]

    # Makes the calculations
    planets = doer(names) 

    fig, axs = plt.subplots(1,3, tight_layout = True, sharex = True,
                           figsize = (9,4))
    
    custom_lines = []
    for planet, color, i in zip(planets, colors, range(len(planets))):
        if i == 0 and bigfirst:
            axs[0].plot(planet.age, planet.mass, color = color, lw=5)
            axs[1].plot(planet.age, planet.radius, color = color, lw=5)
            axs[2].plot(planet.age, planet.lum, color = color, lw=5)
            axs[2].plot(planet.age, planet.lum2, color = colors2[i], linestyle = 'dashed', lw = 5)
                
        axs[0].plot(planet.age, planet.mass, color = color)
        axs[1].plot(planet.age, planet.radius, color = color)
        axs[2].plot(planet.age, planet.lum2, color = color)
        #axs[2].plot(planet.age, planet.lum2, color = colors2[i], linestyle = 'dashed')

        # Legend
        custom_lines.append( Line2D([0], [0], color = colors[i], 
                                    label = labels[i]))
        
    # Make nice
    axs[0].set_ylabel('Mass $[M_\oplus]$', fontsize = 14)
    axs[1].set_ylabel('Radius $[R_\oplus]$', fontsize = 14)
    axs[2].set_ylabel('Luminosity [$L_\odot$]', fontsize = 14)
    
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()

    axs[0].set_xlim(200, 10_000)
    axs[1].set_ylim(7.2,14)
    axs[2].set_ylim(0,1.5e-5)

    fig.suptitle('You\'re hot and you\'re cold', 
                  fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)

    if len(names) == 2:
        box_x = 0.87
    elif len(names) == 3:
        box_x = 0.87
    elif len(names) == 6:
        box_x = 0.83
    elif len(names) == 8:
        box_x = 0.73
    else:
        box_x = 0.97
        
    if len(names) == 8:
        fig.legend(custom_lines, labels,
                    fontsize =  10, ncols = 4, alignment = 'center', # Lawful Neutral
                    bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
    else:
        fig.legend(custom_lines, labels,
                fontsize =  14, ncols = 3, alignment = 'center', # Lawful Neutral
                bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
#%%    
# name = 'jup13'
# name2 = 'jup14'
# name3 = 'jup15'
# name4 = 'jup16'
# name5 = 'jup17'
# name6 = 'jup11'
# labels = ['0.025 AU', '0.035 AU', '0.045 AU', '0.05 AU', '0.1 AU', '0.5 AU']
# plotter([name, name2, name3, name4, name5, name6], 2, labels, False)

# name = 'jup17'
# name2 = 'jupnoirr'
# name3 = 'jupfaraway'
# labels = ['Hot Jup', 'Hot Jup - No Irr', 'Cold Jup', ]
# plotter([name, name2, name3,], 1, labels, True)

# name = 'jupenvtest1'
# name2 = 'jupenvtest2'
# name3 = 'jup17'
# name4 = 'jupenvtest3'
# name5 = 'jupenvtest4'
# name6 = 'jupenvtest5'
# name7 = 'jupenvtest6'
# name8 = 'jupenvtest7'
# labels = ['30', '40', '50', '60', '70', '80', '90', '95']
# plotter([name, name2, name3, name4, name5, name6, name7, name8], 3, labels, False)

