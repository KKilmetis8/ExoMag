#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 15:27:52 2024

@author: konstantinos


Plots the absolute and relative dynamo region
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
# 
import src.prelude as c
from src.Bfield.reynolds import rey_mag, profile_sorter
from src.Bfield.hori import hori_doer

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.start = []
        self.end = []
        self.absolute = []
        self.relative = []
      
    def __call__(self, age, absolute_h, relative_h, start_h, end_h):
        self.age.append(age)
        self.absolute.append(absolute_h)
        self.relative.append(relative_h)
        self.start.append(start_h)
        self.end.append(end_h)
        
def dynsize_doer(names):
    apothikh = []
    for name in names:
        # Instanciate Holder Object
        hold = apothicarios(name)
        
        # Get profiles
        p_path = name # 'data/' + name
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
            try:
                abs_dyn_size = np.abs(R_dynamo_active[0] - R_dynamo_active[-1]) 
            except IndexError:
                abs_dyn_size = 0
                hold(age, 0,0,0,0)
                continue
            rel_dyn_size = 100 * abs_dyn_size / r[0]
            hold(age, abs_dyn_size , rel_dyn_size, 
                 R_dynamo_active[-1], R_dynamo_active[0])
            
        # Keep em
        apothikh.append(hold)
        print('---')
    return apothikh

        
def plotter(names, cols, labels):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'maroon' , 'k', 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c98, c.c99]
    

    
    # Makes the calculations
    planets = dynsize_doer(names) 
    #hori = hori_doer(names, 100)[0]

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
    #axs[0].set_ylim(0,1.2)
    #fig.suptitle(title,  fontsize = 15, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)
    fig.legend(custom_lines, labels,
            fontsize =  9, ncols = len(planets), alignment = 'center', # Lawful Neutral
            bbox_to_anchor=(0.95, -0.03), bbox_transform = fig.transFigure,)
#%%
kind = 'jupenv_orbsep'
if __name__ == '__main__':
    if kind == 'nep_escape':
        pre = '/media/konstantinos/Dunkey/mesadata/lotsaneps/'
        name1 = f'{pre}m5_env0.03_zero_a0.1_s8'
        name2 = f'{pre}m5_env0.06_zero_a0.1_s8'
        name3 = f'{pre}m5_env0.09_zero_a0.1_s8'
        names = [name1, name2, name3,]
        labels = ['3', '6', '9']
        plotter(names, 2, labels)
    if kind == 'jups':
        pre = '/home/konstantinos/Astro-Data/34R-ALM/code/data/'
        name1 = f'{pre}m317_env0.86_zero_a0.1_s8'
        name2 = f'{pre}m317_env0.92_zero_a0.1_s8'
        name3 = f'{pre}m317_env0.96_zero_a0.1_s8'
        names = [name1, name2, name3,]
        labels = ['86', '92', '96']
        plotter(names, 2, labels)
    if kind == 'jupenv_orbsep':
        pre = '/home/konstantinos/Astro-Data/34R-ALM/code/data/'
        name1 = f'{pre}m317_env0.94_zero_a0.1_s8'
        name2 = f'{pre}m317_env0.94_zero_a0.6_s8'
        name3 = f'{pre}m317_env0.94_zero_a0.05_s8'
        names = [name3, name1, name2,]
        labels = ['0.05','0.1', '0.6',]
        plotter(names, 2, labels)
    
    
