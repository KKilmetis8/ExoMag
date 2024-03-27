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
def temperature_scale_height(P, rho, g, grada):
    return P / (rho * g * grada)


def convective_flux(cp, T, rho, vel_conv, P, delta):
    return 2 * cp * T * rho**2 * vel_conv**3 / (P * -delta)


def hardB(p, Rdyns, explore = False):
    r = np.power(10, p.logR)
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
    F = (q / q0)**(2/3) * (L/Ht)**(2/3) * (rho/mean_density)**(1/3) \
        * 4 * np.pi * r**2
    F /= 4*np.pi/3 * r**3
    F_to_the_2by3 = np.abs(np.trapz(F,r))
    #print(F_to_the_2by3**(3/2))
    
    # Explorer return
    if explore:
        # All in SI
        return F_to_the_2by3**(3/2), mean_density, q0
    
    # Calculate B
    B_dyn_squared_SI = 2  * c.mu_0_SI * c.porp_c * mean_density**(1/3) \
                        * q0**(2/3) * F_to_the_2by3 # T
    B_dyn_SI = np.sqrt(B_dyn_squared_SI)
    B_dyn = B_dyn_SI * 10_000 # Tesla to Gauss
    return B_dyn, F_to_the_2by3**(3/2)

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.Bdyn_rey = []
        self.Bdip_rey = []
        self.Bdyn_hori = []
        self.Bdip_hori = []
        self.F_rey = []
        self.F_hori = []

    def __call__(self, age_h, Bdyn_rey, Bdip_rey, F_rey, 
                 Bdyn_hori, Bdip_hori, F_hori):
        self.age.append(age_h)
        self.Bdyn_rey.append(Bdyn_rey)
        self.Bdip_rey.append(Bdip_rey)
        self.Bdyn_hori.append(Bdyn_hori)
        self.Bdip_hori.append(Bdip_hori)
        self.F_rey.append(F_rey)
        self.F_hori.append(F_hori)
        
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
                rey_dyn_end = rey_dyn[0] # it's the wrong way round
                hori_dyn_end = hori_dyn[0] # it's the wrong way round
            except IndexError:
                # Save
                hold(age, hold.Bdyn[-1], hold.Bdip[-1], 0)
                continue  
            B_dyn_rey, F_rey = hardB(p, rey_dyn)
            B_dyn_hori, F_hori = hardB(p, hori_dyn)

            # Get Dipole
            dynamo_rey = (r[0] - rey_dyn_end) / r[0]
            dynamo_hori = (r[0] - hori_dyn_end) / r[0]
            B_dip_rey =  B_dyn_rey * np.power( 1 - dynamo_rey, 3) / np.sqrt(2)
            B_dip_hori =  B_dyn_hori * np.power( 1 - dynamo_hori, 3) / np.sqrt(2)
            print(B_dip_hori)
            # Save
            hold(age, B_dyn_rey, B_dip_rey, F_rey, 
                      B_dyn_hori, B_dip_hori, F_hori)
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
    fig, axs = plt.subplots(1,3, tight_layout = True, sharex = True,
                           figsize = (8,4))
    custom_lines = []
    for planet, i in zip(planets, range(len(planets))):
        axs[0].plot(planet.age, planet.F_rey, c = colors[0])
        axs[0].plot(planet.age, planet.F_hori, c = colors[1], linestyle = '-.')
        axs[1].plot(planet.age, planet.Bdyn_rey, c = colors[0])
        axs[1].plot(planet.age, planet.Bdyn_hori, c = colors[1], linestyle = '-.')
        axs[2].plot(planet.age, planet.Bdip_rey, c = colors[0])
        axs[2].plot(planet.age, planet.Bdip_hori, c = colors[1], linestyle = '-.')
        
    # Make nice
    axs[0].set_ylabel('Efficiency Factor', fontsize = 14)
    axs[1].set_ylabel('Dynamo [G]', fontsize = 14)
    axs[2].set_ylabel('Dipole [G]', fontsize = 14)
    
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
        
    axs[0].set_xlim(500, 11_000)
    axs[0].set_yscale('log')
    if 'jup' in names[0]:
        #axs[0].set_ylim(1_000, 6_000)
        axs[1].set_ylim(0, 350)
        axs[2].set_ylim(0, 350)
    elif 'nep' in names[0]:
        axs[0].set_ylim(1E-12,1E1)
        axs[1].set_ylim(0, 25)
        axs[2].set_ylim(0, 25)
    
    #axs[0].set_yscale('log')
    #axs[1].set_yscale('log')

    fig.suptitle(title, 
                  fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)
    
    # Legend
    custom_lines = [ Line2D([0], [0], color = colors[0], linewidth = 3),
                     Line2D([0], [0], color = colors[1], linewidth = 3, linestyle = '-.'),
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
            fontsize =  9, ncols = 2, alignment = 'left', # Lawful Neutral
            bbox_to_anchor=(0.67, -0.03), bbox_transform = fig.transFigure,)
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
