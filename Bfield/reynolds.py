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
        
def doer(names, many = 8):
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
            r, reynolds_mag_number, age = rey_mag(p)
            hold(age, r, reynolds_mag_number, r[0])
            
        # Keep em
        apothikh.append(hold)
    return apothikh

        
def plotter(names, cols, labels):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'maroon' , 'k', 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = ['dodgerblue', 'forestgreen', 'darkorange']
    elif cols == 4:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c98, c.c99]

    
    # Makes the calculations
    planets = doer(names) 

    fig, axs = plt.subplots(2,4, tight_layout = True, sharex = True,
                           figsize = (8, 6))
    
    custom_lines = []
    i = 0
    for ax in axs.reshape(-1):
        plot_age = 0
        for planet, color, label in zip(planets, colors, labels):
            ax.plot(planet.r[i] / planet.rmax[i], planet.reymag[i], color = color)
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
        
        # Change profile
        i += 1
        
    # # Make nice
    # axs[0].set_ylabel('Dynamo [G]', fontsize = 14)
    # axs[1].set_ylabel('Dipole [G]', fontsize = 14)
    
    # axs[0].grid()
    # axs[1].grid()

    # axs[0].set_xlim(200, 10_000)
    # axs[1].set_xlim(200, 10_000)

    fig.suptitle('Increasingly Puffier Neptunes', 
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
                fontsize =  14, ncols = 3, alignment = 'center', # Lawful Neutral
                bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
 #%%

kind = 'nepenv'
if __name__ == '__main__':
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
    
    if kind == 'puffer':
        name = 'jupenvtest1'
        name2 = 'jupenvtest2'
        name3 = 'jup17'
        name4 = 'jupenvtest3'
        name5 = 'jupenvtest4'
        name6 = 'jupenvtest5'
        name7 = 'jupenvtest6'
        name8 = 'jupenvtest7'
        labels = ['30', '40', '50', '60', '70', '80', '90', '95']
        plotter([name, name2, name3, name4, name5, name6, name7, name8], 4, labels)
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
    
    
