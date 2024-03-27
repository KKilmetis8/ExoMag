#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 20:36:31 2024

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 19:13:08 2024

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import src.prelude as c
from src.Bfield.reynolds import reynolds_doer
from src.Bfield.hori import hori_doer
#%%
       
def plotter(names, title):
    # Makes the calculations
    many = 50
    hori = hori_doer(names, many)[0]
    reynold = reynolds_doer(names, many)[0]
    rey_masks = []
    for h_rs, r_rs in zip(hori.r, reynold.r):
        rey_mask = []
        for r in h_rs:
            if r < np.min(r_rs):
                rey_mask.append(0)
            elif r > np.max(r_rs):
                rey_mask.append(0)
            else:
                rey_mask.append(1)
        rey_masks.append(rey_mask)
        
    plot_age = []
    ylen = len(hori.r[10])
    to_plot_big = np.zeros((ylen, many, 3))
    for i, mask_h, mask_r in zip(range(many), hori.met_hyd, rey_masks):
        plot_age.append(hori.age[i])
        for j, h, r in zip(range(ylen-1, 0, -1), mask_h, mask_r):
            if h and r:
                to_plot_big[j][i] = c.c95_rgb / 255
            elif h and not r:
                to_plot_big[j][i] = c.AEK_rgb / 255
            elif not h and r:
                to_plot_big[j][i] = c.c97_rgb / 255
            elif not h and not r:
                to_plot_big[j][i] = np.array([255,255,255]) / 255
            else:
                to_plot_big[j][i] = np.array([255,255,255]) / 255

    fig, ax = plt.subplots()
    print(to_plot_big.shape)
    ax.pcolormesh(plot_age, np.flip(hori.r[10]), to_plot_big)

    # Minor ticks
    ax.set_xticks(np.linspace(plot_age[0], plot_age[-1], many), minor=True)
    #ax.set_yticks([])
    
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='k', linestyle='--', linewidth=0.5)
    
    #
    labels = ['Neither', 'Only Met Hyd', 'Both' ,'Only Reynolds']
    custom_lines = [ Line2D([0], [0], color = 'k', linewidth = 3),
                    Line2D([0], [0], color = c.AEK, linewidth = 3),
                    Line2D([0], [0], color = c.c95, linewidth = 3),
                    Line2D([0], [0], color = c.c97, linewidth = 3),
        ]
    ax.legend( custom_lines, labels, loc = 'upper right', framealpha = 1,
               edgecolor = 'k',
               bbox_to_anchor=(1.17, 0.65), bbox_transform = fig.transFigure,)
    
    ax.set_xlabel('Age [Myr]', fontsize = 14)
    ax.set_ylabel('r/Radius', fontsize = 14)
    # ax.set_ylim(0, 4500)
    ax.set_xlim(500, 8000)
    ax.set_title('Metallic Hydrogen vs Reynolds ' + title)


    
    #%%    
kind = 'jupenv_zero'
if __name__ == '__main__':
    if kind == 'jupenv_zero':
        name6 = 'jup_e94_zero'
        names = [name6,]
        plotter(names, r'Jupiters 317 $M_\oplus$, Env 94$\%$, 0.1 AU')
    if kind == 'saturn':
        name6 = 'm95_e95_zero_a01_s8'
        names = [name6,]
        plotter(names, r'Saturn 95 $M_\oplus$, Env 95$\%$, Zero Evap, 0.1 AU, S: 8 kb/bar')

