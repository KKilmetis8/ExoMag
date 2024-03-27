#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 20:21:33 2024

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 12:47:45 2024

@author: konstantinos
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.Utilities.profile_sorter import profile_sorter
from src.Bfield.reynolds import dynamo_region
from src.Bfield.hori import hori_met_hyd
#%%

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.r = []
        self.rey_start = []
        self.rey_end = []
        self.hori_start = []
        self.hori_end = []

    def __call__(self, age_h, r, rey_start, rey_end, hori_start, hori_end):
        self.age.append(age_h)
        self.r.append(r)
        self.rey_start.append(rey_start)
        self.rey_end.append(rey_end)
        self.hori_start.append(hori_start)
        self.hori_end.append(hori_end)
        
def hardB_doer(names):
    # Count and generate profile lists
    apothikh = []
    for name in names:
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
            rey_dyn, _, age = dynamo_region(p)
            _, _, _, hori_dyn = hori_met_hyd(p)
            try:
                rey_dyn_end = rey_dyn[0] * c.Rsol / c.Rjup # it's the wrong way round
                hori_dyn_end = hori_dyn[0] # it's the wrong way round
            except IndexError:
                # Save
                hold(age, 0, 0, 0, 0)
                continue  

            # Get Dipole
            # Save
            hold(age, r[0] * c.Rsol / c.Rjup, 
                 rey_dyn[-1] * c.Rsol / c.Rjup, 
                 rey_dyn[0]  * c.Rsol / c.Rjup, 
                 hori_dyn[-1]  * c.Rsol / c.Rjup, 
                 hori_dyn[0] * c.Rsol / c.Rjup)
        apothikh.append(hold)
    return apothikh
    
def plotter(names, cols, title):
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'k', 'maroon',]
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c99, 'k']
    elif cols == 4:
        colors = [c.c91, c.c92, c.c93, c.c95, c.c96, c.c98, c.c99,]
    if cols == 5:
        colors = [c.AEK, 'k', 'r',]

    # Makes the calculations
    planets = hardB_doer(names) 
    fig, axs = plt.subplots(1,1, tight_layout = True, sharex = True,
                           figsize = (4,4))
    custom_lines = []
    for planet, i in zip(planets, range(len(planets))):
        axs.plot(planet.age, planet.rey_start, c = 'k', lw = 3)
        axs.plot(planet.age, planet.rey_end, c = 'k', linestyle = ':', lw = 3)
        axs.plot(planet.age, planet.hori_start, c = 'r')
        axs.plot(planet.age, planet.hori_end, c = 'r', linestyle = ':')
        
    # Make nice
    axs.set_ylabel('Radial Distance from centre [R$_J$]', fontsize = 14)
    axs.grid() 
    #axs.set_xlim(500, 11_000)
    axs.set_xscale('log')
    axs.plot(planet.age, planet.r, c = 'b', linestyle = '-.')
    fig.suptitle(title, fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
             fontsize = 15, transform = fig.transFigure)
    
    # Legend
    custom_lines = [ Line2D([0], [0], color = 'k', linewidth = 3),
                     Line2D([0], [0], color = 'r', linewidth = 3,),
                     ]
    labels = ['Reynolds', 'Met Hyd']
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
            fontsize =  10, ncols = 1, alignment = 'left', # Lawful Neutral
            bbox_to_anchor=(1.27, 0.53), bbox_transform = fig.transFigure,)
    return planets
#%%    
kind = 'saturn'
if __name__ == '__main__':
    if kind == 'jup_e94_zero':
        name = 'jup_e94_zero'       
        names = [name]
        ps = plotter(names, 1, r'Jupiter 317 $M_\oplus$, Env 94$\%$, 0.1 AU')
    if kind == 'saturn':
        name6 = 'm95_e95_zero_a01_s8'
        names = [name6,]
        plotter(names, 1, r'Saturn 95 $M_\oplus$, Env 95$\%$, Zero Evap, 0.1 AU, S: 8 kb/bar')
