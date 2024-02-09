#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 20:49:45 2024

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.Bfield.reynolds import rey_mag, profile_sorter
#%%

def dipole(M, R, L, Rdyn):
    '''
    Parameters
    ----------
    M : Planetary Mass, in Msol
    R : Planetary Radius, in Rsol
    L : Planetary Luminosity, in Lsol
        
    Returns
    -------
    Dipolar Magnetic Field at the pole, in Gauss.

    '''
    B_dyn = 4.8e3 * np.power(M * L**2, 1/6) * np.power(R, -7/6)
    dynamo = (R - Rdyn) / R
    B_dip =  B_dyn * np.power( 1 - dynamo, 3)
    
    return B_dyn, B_dip

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.Bdyn = []
        self.Bdip = []
      
    def __call__(self, age_h, Bdyn_h, Bdip_h):
        self.age.append(age_h)
        self.Bdyn.append(Bdyn_h)
        self.Bdip.append(Bdip_h)
      
def ReyB_doer(names):
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
            Rdyn = R_dynamo_active[0] # it's the wrong way round
            idx = np.argmin(np.abs( r - Rdyn))
            
            # M L R for dynamo surface
            M = p.star_mass
            R = Rdyn
            # print(R * c.Rsol / c.Rjup)
            L = p.luminosity[idx]
            
            # Calc B
            B_dyn, B_dip = dipole(M, R, L, Rdyn)
            
            # Save
            hold(age, B_dyn, B_dip)

        apothikh.append(hold)
    return apothikh


def plotter(names, cols, labels, bigfirst = True):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'maroon' , 'k', 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c99]

    # Makes the calculations
    planets = ReyB_doer(names) 

    fig, axs = plt.subplots(1,2, tight_layout = True, sharex = True,
                           figsize = (6,4))
    
    custom_lines = []
    for planet, color, i in zip(planets, colors, range(len(planets))):
        if i == 0 and bigfirst:
            axs[0].plot(planet.age, planet.Bdyn, color = color, lw=5)
            axs[1].plot(planet.age, planet.Bdip, color = color, lw=5)
                
        axs[0].plot(planet.age, planet.Bdyn, color = color)
        axs[1].plot(planet.age, planet.Bdip, color = color)
        
        # Legend
        custom_lines.append( Line2D([0], [0], color = colors[i], 
                                    label = labels[i]))
        
    # Make nice
    axs[0].set_ylabel('Dynamo [G]', fontsize = 14)
    axs[1].set_ylabel('Dipole [G]', fontsize = 14)
    
    axs[0].grid()
    axs[1].grid()

    axs[0].set_xlim(200, 10_000)
    axs[1].set_xlim(200, 10_000)
    
    axs[0].set_ylim(300, 600)
    axs[1].set_ylim(300, 550)

    fig.suptitle('Increasingly Colder Jupiters', 
                  fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)

    if len(names) == 2:
        box_x = 0.87
    elif len(names) == 3:
        box_x = 0.87
    elif len(names) == 6:
        box_x = 0.83
    elif len(names) == 8:
        box_x = 0.73
    else:
        box_x = 0.97
        
    if len(names) == 8:
        fig.legend(custom_lines, labels,
                    fontsize =  10, ncols = 4, alignment = 'center', # Lawful Neutral
                    bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
    else:
        fig.legend(custom_lines, labels,
                fontsize =  14, ncols = 3, alignment = 'center', # Lawful Neutral
                bbox_to_anchor=(box_x, -0.03), bbox_transform = fig.transFigure,)
#%%    
if __name__ == '__main__':
    # name = 'jup13'
    # name2 = 'jup14'
    # name3 = 'jup15'
    # name4 = 'jup16'
    # name5 = 'jup17'
    # name6 = 'jup11'
    
    # labels = ['0.025 AU', '0.035 AU', '0.045 AU', '0.05 AU', '0.1 AU', '0.5 AU']
    # plotter([name, name2, name3, name4, name5, name6], 2, labels, False)
    
    name = 'jupenvtest1'
    name2 = 'jupenvtest2'
    name3 = 'jup17'
    name4 = 'jupenvtest3'
    name5 = 'jupenvtest4'
    name6 = 'jupenvtest5'
    name7 = 'jupenvtest6'
    name8 = 'jupenvtest7'
    labels = ['30', '40', '50', '60', '70', '80', '90', '95']
    plotter([name, name2, name3, name4, name5, name6, name7, name8], 3, labels, False)
