#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 15:27:52 2024

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
# 
import src.prelude as c
from src.Bfield.reynolds import rey_mag, profile_sorter


class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.absolute = []
        self.relative = []
      
    def __call__(self, age, absolute_h, relative_h):
        self.age.append(age)
        self.absolute.append(absolute_h)
        self.relative.append(relative_h)
        
def dynsize_doer(names):
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
    
        # Do the thing
        for profile in profiles:
            p = mr.MesaData(profile)
            r, reynolds_mag_number, age = rey_mag(p)
            r *=  c.Rsol / c.Rearth
            R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
            abs_dyn_size = np.abs(R_dynamo_active[0] - R_dynamo_active[-1]) 
            rel_dyn_size = 100 * abs_dyn_size / r[0]
            hold(age, abs_dyn_size , rel_dyn_size)
            
        # Keep em
        apothikh.append(hold)
        print('---')

    return apothikh

        
def plotter(names, cols, labels, title):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'maroon' , 'k', 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c98, c.c99]

    
    # Makes the calculations
    planets = dynsize_doer(names) 
    fig, axs = plt.subplots(1,2, tight_layout = True, sharex = True,
                           figsize = (7,4))
    custom_lines = []
    for planet, color, label in zip(planets, colors, labels):
        axs[0].plot(planet.age , planet.absolute, color = color)
        axs[1].plot(planet.age , planet.relative, color = color)
        
        custom_lines.append( Line2D([0], [0], color = color, 
                                    label = label))
    # # Make nice
    axs[0].set_ylabel('Absolute Dynamo Region [R$_\oplus$]', fontsize = 14)
    axs[1].set_ylabel('Relative Dynamo Region [$\%$]', fontsize = 14)
    
    axs[0].grid()
    axs[1].grid()
    
    axs[0].set_xlim(300, 10_000)

    fig.suptitle(title,  fontsize = 15, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)


    fig.legend(custom_lines, labels,
            fontsize =  9, ncols = len(planets), alignment = 'center', # Lawful Neutral
            bbox_to_anchor=(0.95, -0.03), bbox_transform = fig.transFigure,)
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
        plotter(names, 3, labels, 'Jupiter with Diff. Envelopes')
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
    
    
