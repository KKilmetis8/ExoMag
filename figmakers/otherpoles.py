#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 12:42:28 2024

@author: konstantinos

Bonuuuus.

 Okay so this didn't make it in the paper, but this checks out what happens if we dont enforce a dipole geometry. They field strengths for both quadrapoles and octapoles are suprisingly close. This assumes that all the enregy is in one mode. We know (Maus et al., 2006) that the power spectrum of the Earth's dynamo is rather white, that is the same amount of power is contained in the first, like, 32 modes.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.Utilities.profile_sorter import profile_sorter
from src.Bfield.hardB import hardB
from src.Bfield.reynolds import dynamo_region

#%%

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.Bdyn = []
        self.Bdip = []
        self.Bquad = []
        self.Bocto = []
      
    def __call__(self, age_h, Bdyn_h, Bdip_h, Bquad_h, Bocto_h):
        self.age.append(age_h)
        self.Bdyn.append(Bdyn_h)
        self.Bdip.append(Bdip_h)
        self.Bquad.append(Bquad_h)
        self.Bocto.append(Bocto_h)
      
def doer(names):
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
            R_dynamo_active, _, age = dynamo_region(p)
            try:
                Rdyn_end = R_dynamo_active[0] # it's the wrong way round
                # R_dynamo_active *= c.Rsol / c.Rearth
            except IndexError:
                continue 
            B_dyn, F = hardB(p, R_dynamo_active)            
            
            # Calc B
            dynamo = (r[0] - Rdyn_end) / r[0]
            B_dip =  B_dyn * np.power( 1 - dynamo, 3)
            B_quad =  B_dyn * np.power( 1 - dynamo, 4)
            B_octo =  B_dyn * np.power( 1 - dynamo, 5)
            
            # Save
            hold(age, B_dyn, B_dip, B_quad, B_octo)

        apothikh.append(hold)
    return apothikh


def plotter(planets, cols, labels, bigfirst = True):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'k' , 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
        
    # Makes the calculations
    # planets = doer(names) 

    fig, ax = plt.subplots(1,1, tight_layout = True,
                           figsize = (4,4))

    for planet in planets:
        ax.plot(planet.age, planet.Bdip, color = colors[0], label = 'Dipole')
        ax.plot(planet.age, planet.Bquad, color = colors[1], label = 'Quad')
        ax.plot(planet.age, planet.Bocto, color = colors[2], label = 'Octo')
        
    # Make nice
    ax.set_ylabel(r'$B_\mathrm{surface}^\mathrm{(max)}$ [G]', fontsize = 14)
    ax.set_xlabel('Age [Myr]', fontsize = 14)
    ax.set_xlim(300, 2800)
    ax.set_ylim(0, 15)
    plt.savefig('figs/multipole.pdf', format='pdf', dpi = 300)
    # ax.legend()
    #fig.suptitle('Other poles, Hot Jupiter', 
    #              fontsize = 18, y = 0.98)

#%%    
name = 'm17_env0.04_zero_a0.1_s8'
labels = ['j']
planets  = doer([name])
plotter(planets, 1, labels)



