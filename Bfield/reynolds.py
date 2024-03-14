#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 05:49:18 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
# 
import src.prelude as c
def profile_sorter(profiles):
    ''' yes this is N^2, yes I could make this better, this is not a problem'''
    zs = np.arange(1, len(profiles)+1) # Zahlen!
    new_profiles = []
    for z in zs: 
        z = str(z)
        z = 'e' + z + '.' # don't mix 10 with 1
        profile = [s for s in profiles if z in s] 
        new_profiles.append(profile[0]) # only 1 match idealy
    return new_profiles

def rey_mag(p):

    # Data Out
    r = np.power(10, p.logR)
    Temperature = np.power(10, p.logT) # [K]    
    conv_len = p.mlt_mixing_length # [cm]
    conv_vel = p.conv_vel # [cm/s]
    
    # Mag Reynolds
    conductivity = 1e7 * Temperature**(3/2)
    mag_diffusivity = c.inv_mu_0_cgs / conductivity
    reynolds_mag =  conv_len * conv_vel / mag_diffusivity
    
    age = np.round(p.star_age /  1e6, 2)
    return r, reynolds_mag, age
#%%
def dynamo_region(p, normalize = False):
    r, rey_mag_num, age = rey_mag(p)
    regions = []
    current_region = None
    for i, rmn in enumerate(rey_mag_num):
        if rmn > c.critical_rey_mag_num:
            if current_region is None:
                current_region = [i] # Make a new region
            else:
                current_region.append(i)
                
        # Finish this region
        elif current_region is not None:
            regions.append(current_region)
            current_region = None
            
    # Save the last one
    if current_region is not None:
       regions.append(current_region)
       
    # Return the largest one
    dyn_region_len = 0 # arbitrary
    for i, region in enumerate(regions):
        temp_len = len(region)
        if temp_len > dyn_region_len:
            dyn_region_idx = i
            dyn_region_len = temp_len
    if dyn_region_len > 0:
        dynamo_region = regions[dyn_region_idx]
    else: 
        return [], [], age
    
    if normalize:
        return r[dynamo_region] / r[0], rey_mag_num[dynamo_region], age
    else:
        return r[dynamo_region], rey_mag_num[dynamo_region], age

#%%

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.r = []
        self.reymag = []
        self.rmax = []
      
    def __call__(self, age, r, raymag, rmax):
        self.age.append(age)
        self.r.append(r)
        self.reymag.append(raymag)
        self.rmax.append(rmax)
        
def doer(names, many = 16):
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
            r, reynolds_mag_number, age = dynamo_region(p, normalize = True)
            hold(age, r, reynolds_mag_number, 0)
            
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
    planets = doer(names) 

    fig, axs = plt.subplots(4,4, tight_layout = True, sharex = True,
                           figsize = (8, 6))
    
    custom_lines = []
    i = 0
    for ax in axs.reshape(-1):
        plot_age = 0
        for planet, color, label in zip(planets, colors, labels):
            ax.plot(planet.r[i], planet.reymag[i], color = color)
            # ax.axvline(planet.rmax[i], color = 'k', linestyle = 'dashed')
            ax.set_yscale('log')
            plot_age += planet.age[i]
            
            # Legend
            if i == 0:
                custom_lines.append( Line2D([0], [0], color = color, 
                                            label = label))
        # Critical Reymag
        ax.axhline(50, linestyle = 'dashed', c = 'k')
        
        # Avg age
        plot_age /= len(planets)
        plot_age = np.round(plot_age, decimals = 0)
        ax.set_title(r"Age $\sim$" + str(plot_age) + ' Myrs')
        ax.set_xlim(0.8, 1.1)
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
kind = 'jups_thekink'
if __name__ == '__main__':
    
    if kind == 'jups':
        name = 'jup13'
        name2 = 'jup14'
        name3 = 'jup15'
        name4 = 'jup16'
        name5 = 'jup17'
        name6 = 'jup11'
        labels = ['0.025 AU', '0.035 AU', '0.045 AU', '0.05 AU', '0.1 AU', '0.5 AU']
        plotter([name, name2, name3, name4, name5, name6], 2, labels)
    if kind == 'jups_thekink':
        name = 'jup_e94_zero_a01_s8'
        name2 = 'jup_e94_zero_a01'
        labels = ['94%', '92%',]
        plotter([name, name2,], 1, labels, 'jups_thekink')
    if kind == 'noirr':
        name = 'jup17'
        name2 = 'jupfaraway'
        name3 = 'jupnoirr'
        labels = ['Hot Jup', 'Cold Jup', 'Hot Jup - No Irr']
        plotter([name, name2, name3], 1, labels)
    if kind == 'jupenv_zero':
        name3 = 'jup_e85_zero'
        name4 = 'jup_e90_zero'
        name5 = 'jup_e92_zero'
        name6 = 'jup_e94_zero'
        name7 = 'jup_e96_zero'
        #name8 = 'jup_e97_zero'
        #name9 = 'jup3_e98_zero'
        names = [name3, name4, name5, name6, name7]#, name4, name9]
        labels = ['85', '90', '92', '94', '96',]#, '95', '98']
        plotter(names, 4, labels, 'Increasingly Puffier Jupiters')
    if kind == 'nepenv':
        name1 = 'nep_1'
        name2 = 'nep_2'
        name3 = 'nep_3'
        name4 = 'nep_4'
        name5 = 'nep_5'
        name6 = 'nep_6'
        names = [name1, name2, name3, name4, name5, name6]
        labels = ['10', '20', '30', '40', '50', '60']
        plotter(names, 2, labels)#, 'Neptunes with Diff. Envelopes')
    
    
