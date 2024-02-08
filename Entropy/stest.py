#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 12:04:56 2024

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
# 
import src.prelude as c

def plotter(names, cols, labels):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'maroon' , 'k', 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]*2          
    elif cols == 3:
        colors = ['dodgerblue', 'forestgreen', 'darkorange']
        
    fig, axs = plt.subplots(2,1, tight_layout = True,
                           figsize = (7,4))
    custom_lines = []
    ax = axs[0]
    ax.grid()
    ax.set_title('250 $M_\oplus$, $f_{a,0}$ 50$\%$')
    ax.set_xscale('log')
    ls = 'solid'
    for name, color, label in zip(names, colors, labels):
        print(name)
        # Data load
        path = 'data/' + name + '/'
        h = mr.MesaData(path + 'history_7.data')
        age = h.star_age / 1e6 # Myr
        rp = h.rp_rearth
        
        if name == 'stest21': # single use bullshit
            ls = 'dashdot'
            ax = axs[1]
            
        # Plot
        ax.plot(age, rp, color = color, linestyle = ls)
        
        # For the legend
        custom_lines.append( Line2D([0], [0], color = color, label = label,
                                    linestyle = ls))

    # Make Pretty
    ax.set_title('250 $M_\oplus$, $f_{a,0}$ 25$\%$')
    ax.set_xscale('log')
    ax.grid()

    #ax.set_xlim(10, 10_000)
    #ax.set_ylim(8.5, 10)
    fig.suptitle(r'K$\&$V Fig.4', 
                  fontsize = 15, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)
    fig.text(-0.02, 0.3, 'Planetary Radius $[R_\oplus]$', rotation = 'vertical', 
             fontsize = 15, transform = fig.transFigure)
    
    if len(names) == 2:
        box_x = 0.87
    elif len(names) == 3:
        box_x = 0.97
    elif len(names) == 6:
        box_x = 0.92
    elif len(names) == 12:
        box_x = 0.85
    else:
        box_x = 0.97
        
    # fig.legend(custom_lines, labels,
    #             fontsize = 12, ncols = 3, alignment = 'center', # Lawful Neutral
    #             bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
#%%

kind = 'stest'

if kind == 'planets':
    name = 'neptest'
    name2 = 'supearth'
    name3 = 'juptest'
    labels_p = [r'Superearth 2$M_\oplus$, 1$\%$', 
                              r'Neptune 17$M_\oplus$, 30$\%$', 
                          r'Jupiter 317$M_\oplus$, 50$\%$']
    plotter([name, name2, name3], 3, labels_p)

if kind == 'jups':
    name = 'jup13'
    name2 = 'jup14'
    name3 = 'jup15'
    name4 = 'jup16'
    name5 = 'jup17'
    name6 = 'jup11'
    labels = ['0.025 AU', '0.035 AU', '0.045 AU', '0.05 AU', '0.1 AU', '0.5 AU']
    plotter([name, name2, name3, name4, name5, name6], 2, labels)
    
if kind == 'noirr':
    name = 'jup17'
    name2 = 'jupfaraway'
    name3 = 'jupnoirr'
    labels = ['Hot Jup', 'Cold Jup', 'Hot Jup - No Irr']
    plotter([name, name2, name3], 1, labels)

if kind == 'stest':
    testno = 6
    nos = np.arange(1, testno+1)
    names1 = [ 'stest' + str(no) for no in nos]
    names2 = [ 'stest2' + str(no) for no in nos]
    names = names1 + names2
    
    entropies = np.arange(6.5, 9+0.5, 0.5)
    labels1 = ['50\% $f_{a,0}$ -' + str(s) + r' $k_B$/bar.' for s in entropies]
    labels2 = ['25\% $f_{a,0}$ -' + str(s) + r' $k_B$/bar.' for s in entropies]
    labels = labels1 + labels2
    plotter(names, 2, labels)
