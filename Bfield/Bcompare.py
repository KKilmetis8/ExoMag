#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 20:46:26 2024

@author: konstantinos
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import src.prelude as c
from src.Bfield.B_reynolds import ReyB_doer
from src.Bfield.mediumB import medB_doer
from src.Bfield.hardB import hardB_doer

def plotter(names, names2, names3, cols, labels, title, bigfirst = True):
    
    # Specify Palettes
    if cols == 1:
        colors = ['k', c.AEK, 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c99]

    # Makes the calculations
    planets1 = ReyB_doer(names2)
    planets2 = medB_doer(names) 
    planets3 = hardB_doer(names3)
    planets = planets1 + planets2 + planets3
    fig, axs = plt.subplots(1,2, tight_layout = True, sharex = True, 
                            sharey = True, figsize = (6,4))
    
    custom_lines = []
    for planet, color, i in zip(planets, colors, range(len(planets))):
        if i == 1 and bigfirst:
            axs[0].plot(planet.age, planet.Bdyn, color = color, lw=5)
            axs[1].plot(planet.age, planet.Bdip, color = color, lw=5)
                
        axs[0].loglog(planet.age, planet.Bdyn, color = color)
        axs[1].loglog(planet.age, planet.Bdip, color = color)
        
        # Legend
        custom_lines.append( Line2D([0], [0], color = colors[i], 
                                    label = labels[i]))
        
    # Make nice
    axs[0].set_ylabel('Dynamo [G]', fontsize = 14)
    axs[1].set_ylabel('Dipole [G]', fontsize = 14)
    
    axs[0].grid()
    axs[1].grid()
    axs[0].set_xlim(450, 10_000)
    if 'jup' in names[0]:
        axs[0].set_ylim(1e0, 550)
        axs[1].set_ylim(1e0, 550)
    elif 'nep' in names[0]:
        axs[0].set_ylim(0, 45)
        axs[1].set_ylim(0, 45)
    #axs[1].set_ylim(0, 200)

    # axs[0].set_ylim(300, 600)
    # axs[1].set_ylim(300, 550)
    # axs[0].set_yscale('log')
    #axs[1].set_yscale('log')

    fig.suptitle(title + ' B comparison', 
                  fontsize = 18, y = 0.98)
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
        
    if len(planets) == 8:
        fig.legend(custom_lines, labels,
                    fontsize =  10, ncols = 4, alignment = 'center', # Lawful Neutral
                    bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
    else:
        fig.legend(custom_lines, labels,
                fontsize =  14, ncols = 3, alignment = 'center', # Lawful Neutral
                bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
#%%    
kind = 'jup_e94_zero_a01_s9'
if kind == 'jup_e94_zero_a01_s7':
    name = 'jup_e94_zero_a01_s7'
    labels = ['Scaling', 'F=1', 'Hard']
    plotter([name], [name], [name], 1, labels, 
            r'Jupiter $317 M_\oplus$ 94\% $\alpha$=0.1 AU S=7 kb/bar', False)
if kind == 'jup_e94_zero_a01_s8':
    name = 'jup_e94_zero_a01_s8'
    labels = ['Scaling', 'F=1', 'Hard']
    plotter([name], [name], [name], 1, labels, 
            r'Jupiter $317 M_\oplus$ 94\% $\alpha$=0.1 AU S=8 kb/bar', False)
if kind == 'jup_e94_zero_a01_s9':
    name = 'jup_e94_zero_a01_s9'
    labels = ['Scaling', 'F=1', 'Hard']
    plotter([name], [name], [name], 1, labels, 
            r'Jupiter $317 M_\oplus$ 94\% $\alpha$=0.1 AU S=9 kb/bar', False)
if kind == 'jup_e94_zero_autoS':
    name = 'jup_e94_zero_autoS'
    labels = ['Scaling', 'F=1', 'Hard']
    plotter([name], [name], [name], 1, labels, 
            r'Jupiter $317 M_\oplus$ 94\% $\alpha$=0.1 AU autoS', False)
##
if kind == 'jup_e94_zero_a001_s8':
    name = 'jup_e94_zero_a001_s8'
    labels = ['Scaling', 'F=1', 'Hard']
    plotter([name], [name], [name], 1, labels, 
            r'Jupiter $317 M_\oplus$ 94\% $\alpha$=0.01 AU S=8 kb/bar', False)
if kind == 'nep_e4_zero_7s':
    name = 'nep_e4_zero_7s'
    labels = ['Scaling', 'F=1', 'Hard']
    plotter([name], [name], [name], 1, labels, r'Neptune $17 M_\oplus$ 4\% $\alpha$=0.1 AU ', False)
