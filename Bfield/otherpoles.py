#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 12:42:28 2024

@author: konstantinos
"""

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
from src.reynolds import rey_mag, profile_sorter
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

def quadpole(M, R, L, Rdyn):
    B_dyn = 4.8e3 * np.power(M * L**2, 1/6) * np.power(R, -7/6)
    dynamo = (R - Rdyn) / R
    B_quad =  B_dyn * np.power( 1 - dynamo, 4)
    
    return  B_quad

def octopole(M, R, L, Rdyn):
    B_dyn = 4.8e3 * np.power(M * L**2, 1/6) * np.power(R, -7/6)
    dynamo = (R - Rdyn) / R
    B_octo =  B_dyn * np.power( 1 - dynamo, 5)
    
    return  B_octo

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
            r, reynolds_mag_number, age = rey_mag(p)
            
            # Get Rdyn
            R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
            Rdyn = R_dynamo_active[0] # it's the wrong way round
            
            # Find closest age in history file
            # idx = np.argmin(np.abs(h_age - age))
            M = p.star_mass
            R = r[0]
            # print(R * c.Rsol / c.Rjup)
            L = p.photosphere_L
            
            # Calc B
            B_dyn, B_dip = dipole(M, R, L, Rdyn)
            B_quad = quadpole(M, R, L, Rdyn)
            B_octo = octopole(M, R, L, Rdyn)
            
            # Save
            hold(age, B_dyn, B_dip, B_quad, B_octo)

        apothikh.append(hold)
    return apothikh


def plotter(names, cols, labels, bigfirst = True):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'k' , 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
        
    # Makes the calculations
    planets = doer(names) 

    fig, ax = plt.subplots(1,1, tight_layout = True,
                           figsize = (4,4))

    for planet in planets:
        ax.plot(planet.age, planet.Bdip, color = colors[0], label = 'Dipole')
        ax.plot(planet.age, planet.Bquad, color = colors[1], label = 'Quad')
        ax.plot(planet.age, planet.Bocto, color = colors[2], label = 'Octo')
        
    # Make nice
    ax.set_ylabel('Magnetic Field [G]', fontsize = 14)
    ax.set_xlabel('Age [Myr]', fontsize = 14)
    ax.set_xlim(200, 10_000)
    ax.legend()
    fig.suptitle('Other poles, Hot Jupiter', 
                  fontsize = 18, y = 0.98)

#%%    
name = 'jup17'
labels = ['Hot Jup']
plotter([name], 1, labels)



