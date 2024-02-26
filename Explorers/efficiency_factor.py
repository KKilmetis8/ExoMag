#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 16:51:42 2024

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 14:34:36 2024

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.Bfield.reynolds import profile_sorter, rey_mag
from src.Bfield.hardB import hardB

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.F = []
        self.avgrho = []
        self.q0 = []
      
    def __call__(self, age_h, F_h, avgrho_h, q0_h):#), lum2_h):
        self.age.append(age_h)
        self.F.append(F_h)
        self.avgrho.append(avgrho_h)
        self.q0.append(q0_h)
      
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
            r, reynolds_mag_number, age = rey_mag(p)
            
            # Get Rdyn
            R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
            try:
                Rdyn_end = R_dynamo_active[0] # it's the wrong way round
            except IndexError:
                # Save
                hold(age, 0, 0, 0)
                continue
            F, avgrho, q0 = hardB(p, R_dynamo_active, explore = True)
            hold(age, F, avgrho / 1000, q0)#, L2)
        apothikh.append(hold)
    return apothikh


def plotter(names, cols, labels, title, 
            bigfirst = False, met_hyd_lines = False):
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

    fig, axs = plt.subplots(1,1, tight_layout = True, sharex = True,
                           figsize = (5,4))
    
    custom_lines = []
    for planet, color, i in zip(planets, colors, range(len(planets))):
        axs.plot(planet.age, planet.F, color = color)
        # Legend
        custom_lines.append( Line2D([0], [0], color = colors[i], 
                                    label = labels[i]))
        
    # Make nice
    axs.set_ylabel('Efficiency Factor', fontsize = 14)
    axs.grid()
    axs.set_xlim(500, 10_000)
    axs.set_yscale('log')

    fig.suptitle(title, 
                  fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)
    fig.legend(custom_lines, labels,
            fontsize =  9, ncols = len(planets) // 2, alignment = 'center', # Lawful Neutral
            bbox_to_anchor=(0.8, -0.03), bbox_transform = fig.transFigure,)  
#%%    
kind = 'jupenv'
if __name__ == '__main__':
    if kind == 'jupenv':
        name1 = 'jup_e30'
        name2 = 'jup_e40'
        name3 = 'jup_e50'
        name4 = 'jup_e60'
        name5 = 'jup_e70'
        name6 = 'jup_e80'
        name7 = 'jup_e90'
        name8 = 'jup_e95'
        name9 = 'jup_e98'
        names = [name1, name2, name3, name4, name5, name6, name7, name8, name9]
        labels = ['30', '40', '50', '60', '70', '80', '90', '95', '98']
        plotter(names, 3, labels, 'Jupiter with Diff. Envelopes', 
                met_hyd_lines=True)
    if kind == 'nepenv':
        name1 = 'nep_1'
        name2 = 'nep_2'
        name3 = 'nep_3'
        name4 = 'nep_4'
        name5 = 'nep_5'
        name6 = 'nep_6'
        names = [name1, name2, name3, name4, name5, name6]
        labels = ['10', '20', '30', '40', '50', '60']
        plotter(names, 3, labels, 'Neptunes with Diff. Envelopes')
    if kind == 'se5env':
        name1 = 'se58'
        name2 = 'se59'
        name3 = 'se510'
        name4 = 'se511'
        name5 = 'se512'
        name6 = 'se513'
        name7 = 'se514'
        names = [name1, name2, name3, name4, name5, name6, name7]
        labels = ['1', '2', '6', '12', '18', '24', '30']
        plotter(names, 3, labels, '5xEarth with Diff. Envelopes')
    if kind == 'se5envEL':
        name1 = 'se5_e001_a03_EL'
        name2 = 'se5_e002_a03_EL'
        name3 = 'se5_e005_a03_EL'
        name4 = 'se5_e01_a03_EL'
        name5 = 'se5_e013_a03_EL'
        name6 = 'se5_e017_a03_EL'
        names = [name1, name2, name3, name4, name5, name6]
        labels = ['1', '2', '5', '10', '13', '17',]
        plotter(names, 3, labels, '5xEarth with Diff. Envelopes')
    if kind == 'jup-sep':
        name1 = 'jup_e95_a002'
        name3 = 'jup_e95_a008_2'
        name4 = 'jup_e95_a012'
        name5 = 'jup_e95_a016'
        name6 = 'jup_e95_a020'
        name7 = 'jup_e95_a024'
        name8 = 'jup_e95_a028'
        names = [name1, name3, name4, name5, name6, name7, name8]
        labels = ['0.02', '0.08', '0.12', '0.16', '0.20', '0.24', '0.28',]
        plotter(names, 3, labels, 'Jups far away')