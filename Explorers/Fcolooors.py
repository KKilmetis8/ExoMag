#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:45:21 2024

@author: konstantinos
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
# 
import src.prelude as c
from src.Bfield.reynolds import profile_sorter, rey_mag
from src.Bfield.hardB import convective_flux, temperature_scale_height

#%%
def F_plus_friends(p, Rdyns):
    r = np.power(10, p.logR) * c.Rsol / c.Rearth
    start = np.argmin(np.abs(r - Rdyns[-1])) # inner surface
    end = np.argmin(np.abs(r - Rdyns[0])) # outer surface
    
    # Convective Flux
    T = np.power(10, p.logT)[end:start] # surface -> center
    P = np.power(10, p.logP)[end:start]
    rho = np.power(10, p.logRho)[end:start]
    delta = p.dlnRho_dlnT_const_Pgas[end:start]
    conv_vel = p.conv_vel[end:start]
    cp = p.cp[end:start]
    q = convective_flux(cp, T, rho, conv_vel, P, delta)
    q *= 1e-3 # mW/m^2 to W/m^2

    # Reference flux
    q0 = q[-1] * Rdyns[-1]**2 / Rdyns[0]**2
    
    # Mean Denstity
    mean_density = np.mean(rho) # in the dyn region
    mean_density *= 1000 # g/cm^3 to kg/m^3

    # Temperature Scale Height
    g = p.grav[end:start]
    grada = p.grada[end:start]
    Ht = temperature_scale_height(P, rho, g, grada)
    Ht *= 1e-2 # cm to m
    
    # L
    L = p.mlt_mixing_length[end:start]
    L *= 1e-2 # cm to m
    
    # Convert to SI
    rho *= 1000 # g/cm^3 -> kg/m^3
    r = r[end:start] * c.Rsol * 1e-2 # Rsol to cm to m   
    # Calculate effieciency Factor
    red = (q / q0)**(2/3)
    blue = (L / Ht)**(2/3)
    yellow = (rho/mean_density)**(1/3)
    F = red * blue * yellow
    # F_to_the_2by3 = np.abs(np.trapz(F,r))
    #print(F_to_the_2by3**(3/2))

    return F, red, blue, yellow


class apothicarios:
    def __init__(self, name, age, rmin, rmax):
        self.name = name
        self.age = age
        self.r = []
        self.F = []
        self.col = []
        self.rmin = rmin
        self.rmax = rmax
        self.reds = []
        self.blues = []
        self.yellows = []
        
    def __call__(self, r, F, col, red, blue, yellow):
        self.r.append(r)
        self.F.append(F)
        self.col.append(col)
        self.reds.append(red)
        self.blues.append(blue)
        self.yellows.append(yellow)
        
def doer(name, many = 9):
    apothikh = []
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
        r *= c.Rsol / c.Rearth
        # Get Rdyn
        R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
        try:
            Rdyn_end = R_dynamo_active[-1] # it's the wrong way round
        except IndexError:
            # Save
            hold = apothicarios(name, age, r[-1], r[0])
            apothikh.append(hold)
            continue
        
        hold = apothicarios(name, age, r[-1], r[0]) # Instanciate Holder Object
        Fs, reds, blues, yellows = F_plus_friends(p, R_dynamo_active)
        col = []
        for red, blue, yellow, F, r in zip(reds, blues, yellows, Fs, R_dynamo_active):
            if red > blue and red > yellow:
                if red > 10*blue and red > 10*yellow:
                    col = 'red'
                elif red > 10*blue and red < 10*yellow:
                    col = 'orange'
                elif red < 10* blue and red > 10*yellow:
                    col = 'purple'
                else:
                    col = 'grey'
                    
            if blue > red and blue > yellow:
                if blue > 10*red and blue > 10*yellow:
                    col = 'blue'
                elif blue > 10*red and blue < 10*yellow:
                    col = 'green'
                elif blue < 10* red and blue > 10*yellow:
                    col = 'purple'
                else:
                    col = 'grey'
                
            if yellow > red and yellow > blue:
                if yellow > 10*red and yellow > 10*blue:
                    col = 'yellow'
                elif yellow > 10*red and yellow < 10*blue:
                    col = 'green'
                elif yellow < 10* red and yellow > 10*blue:
                    col = 'orange'  
                else:
                    col = 'grey'
            # Save
            hold(r, F, col, red, blue, yellow)
        # Keep em
        apothikh.append(hold)
    return apothikh
#%%
fig, axs = plt.subplots(2, 4, figsize = (12, 5), tight_layout = True, 
                        sharex = True)
planet = doer('jup_e90_zero')
fig.suptitle('Jupiter 317 $M_\oplus$, 90 $\%$ $f_{env}$, 0.1 AU, No Escape', 
             fontsize = 20)
for ax, prof in zip(axs.reshape(-1), range(1,9)):
    ax.plot(planet[prof].r, planet[prof].F, c ='k')
    ax.axvline(planet[prof].rmin, linestyle = '--', c = 'k')
    ax.axvline(planet[prof].rmax, linestyle = '--', c = 'k')

    ax.set_xlim(0.9 * planet[prof].rmin,)
    ax.set_title(str(planet[prof].age) + ' Myrs')
    for i in range(1, len(planet[prof].r)):

        ax.axvspan(planet[prof].r[i-1], planet[prof].r[i], 
                    facecolor=planet[prof].col[i], alpha = 0.2)
    # Colors
    ax2 = ax.twinx()
    ax2.plot(planet[prof].r, planet[prof].yellows, c = c.AEK)
    ax2.plot(planet[prof].r, planet[prof].blues, c = c.c97)
    ax2.plot(planet[prof].r, planet[prof].reds, c = c.c91)
    
    ax2.set_yscale('log')
    
#fig.text($R_\oplus$]')
fig.text(-0.008, 0.42, r'$F^{2/3}$', fontsize = 20, rotation = 90,
         transform = fig.transFigure)
fig.text(1.001, 0.35, r'Color Terms', fontsize = 20, rotation = 270,
         transform = fig.transFigure)
fig.text(0.45, -0.03, r'$R [R_\oplus]$', fontsize = 20, rotation = 0,
         transform = fig.transFigure)

labels = ['Convection', 'Length Scale', 'Density' ]
lines = [Line2D([0], [0], color = c.c91,),
         Line2D([0], [0], color = c.c97,),
         Line2D([0], [0], color = c.AEK,)]
fig.legend(lines, labels,
        fontsize = 15, ncols = 3,
        bbox_to_anchor=(0.7, -0.03), bbox_transform = fig.transFigure,)