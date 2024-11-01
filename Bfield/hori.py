#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 19:13:08 2024

@author: konstantinos

Uses the metallic hydrogen criterion first found Hori 2020 for the dynamo criterion
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.Utilities.profile_sorter import profile_sorter
from src.Bfield.reynolds import rey_mag, dynamo_region
#%%
def hori_met_hyd(p, normalize = False, rearth = False):
    # Extract form MESA.
    r = np.power(10, p.logR) # Rsol
    if normalize:
        r *= 1/r[0]
    if rearth:
        r *= c.Rsol/c.Rearth
    Temperature = np.power(10, p.logT) # [K]  
    Pressure = np.power(10, p.logP) # erg / cm^3
    conv_vel = p.conv_vel # cm/s
    

    # Metallic hydrogen conditions
    critical_T = 2_000 # [K]
    bar_to_cgs = 1_000_000
    critical_P = 2 * 1e6 * bar_to_cgs # 2 MBars
    good_T = np.ma.MaskedArray(Temperature > critical_T)
    good_P = np.ma.MaskedArray(Pressure > critical_P )
    good_vel = np.ma.MaskedArray(conv_vel != 0)
    
    met_hyd_mask = np.multiply(good_T, good_P)
    met_hyd_mask = np.multiply(met_hyd_mask, good_vel)
    # Apply mask
    dynamo_region = r[met_hyd_mask]
    age = np.round(p.star_age /  1e6, 2)
    return r,  met_hyd_mask, age, dynamo_region
          
class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.r = []
        self.met_hyd = []
        self.rdyn = []

    def __call__(self, age, r_h, met_hyd_h, rdyn_h):
        self.age.append(age)
        self.r.append(r_h)
        self.met_hyd.append(met_hyd_h)
        self.rdyn.append(rdyn_h)
        
def hori_doer(names, many = 8):
    apothikh = []
    for name in names:
        # Instanciate Holder Object
        hold = apothicarios(name)
        
        # Get profiles
        p_path = 'data/' + name
        profiles = os.popen('ls ' + p_path + '/profile*.data').read()
        profiles = list(profiles.split("\n"))
        profiles.pop() # Remove last
        profiles = profile_sorter(profiles) # Guess what that does
    
        # Choose how many profiles you want
        step = len(profiles) // many
        indices = np.arange(0, len(profiles), step, dtype=int)
        # Do the thing
        for i in indices:
            p = mr.MesaData(profiles[i])
            r, met_hyd, age, rdyn = hori_met_hyd(p, normalize = True)
            hold(age, r, met_hyd, rdyn)
            
        # Keep em
        apothikh.append(hold)
    return apothikh
       
def plotter(names, cols, labels, title):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'maroon' , 'k', 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = ['dodgerblue', 'forestgreen', 'darkorange']
    elif cols == 4:
        colors = [c.c91, c.c92, c.c93, c.c95, c.c96, c.c98, c.c99,]

    # Makes the calculations
    planets = hori_doer(names) 
    fig, axs = plt.subplots(4,2, tight_layout = True, sharex = True,
                           figsize = (8, 6))
    custom_lines = []
    i = 0
    for ax in axs.reshape(-1):
        plot_age = 0
        j = 1
        for planet, color, label in zip(planets, colors, labels):
            ax.plot(planet.r[i], planet.met_hyd[i]*j, color = color)
            j += 1
            # ax.axvline(planet.rmax[i], color = 'k', linestyle = 'dashed')
            # ax.set_yscale('log')
            plot_age += planet.age[i]
            
            # Legend
            if i == 0:
                custom_lines.append( Line2D([0], [0], color = color, 
                                            label = label))
        # Critical Reymag
        # ax.axhline(50, linestyle = 'dashed', c = 'k')
        
        # Avg age
        plot_age /= len(planets)
        plot_age = np.round(plot_age, decimals = 0)
        ax.set_title(r"Age $\sim$" + str(plot_age) + ' Myrs')
        # ax.set_xlim(0.8, 1.01)
        # Change profile
        i += 1
        
    # # Make nice
    # axs[0].set_ylabel('Dynamo [G]', fontsize = 14)
    # axs[1].set_ylabel('Dipole [G]', fontsize = 14)
    
    # axs[0].grid()
    # axs[1].grid()

    # axs[0].set_xlim(200, 10_000)
    # axs[1].set_xlim(200, 10_000)

    fig.suptitle('Increasingly Puffier Jupiters', 
                  fontsize = 15, y = 0.98)
    fig.text(0.47, -0.02, r' r/Radius', 
              fontsize = 15, transform = fig.transFigure)
    fig.text(-0.02, 0.3, 'Magnetic Reynolds Number', rotation = 'vertical', 
             fontsize = 15, transform = fig.transFigure)

    if len(names) == 2:
        box_x = 0.87
    elif len(names) == 3:
        box_x = 0.87
    elif len(names) == 6:
        box_x = 0.83
    elif len(names) == 8:
        box_x = 0.60
    else:
        box_x = 0.97
        
    if len(names) == 8:
        fig.legend(custom_lines, labels,
                    fontsize =  18, ncols = 4, alignment = 'center', # Lawful Neutral
                    bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
    else:
        fig.legend(custom_lines, labels,
                fontsize =  14, ncols = 5, alignment = 'center', # Lawful Neutral
                bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
    
#%%    
kind = 'nep'
if __name__ == '__main__':
    if kind == 'nep':
        name3 = 'm17_env0.06_zero_a0.1_s8'

        names = [name3,]
        labels = ['85',]
        plotter(names, 4, labels, 'nep ')

