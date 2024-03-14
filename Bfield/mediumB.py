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
from src.Utilities.profile_sorter import profile_sorter
from src.Bfield.reynolds import rey_mag, dynamo_region
from src.Bfield.B_reynolds import ReyB_doer
#%%

def convective_flux(cp, T, rho, vel_conv, P, delta):
    return 2 * cp * T * rho**2 * vel_conv**3 / (P * -delta)

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
    
def medB_doer(names):
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
            r = np.power(10, p.logR)
            R_dynamo_active, rmn, age = dynamo_region(p)
            try:
                Rdyn_end = R_dynamo_active[0] # it's the wrong way round
                Rdyn_start = R_dynamo_active[-1]
                # R_dynamo_active *= c.Rsol / c.Rearth
            except IndexError:
                # Save
                hold(age, hold.Bdyn[-1], hold.Bdip[-1])
                continue  
            
            start = np.argmin(np.abs(r - Rdyn_start))
            end = np.argmin(np.abs(r - Rdyn_end))
            
            # Neccecary quantities for B | NOTE: Assumes F=1
            density = np.power(10, p.logRho)
            mean_density = np.mean(density[end:start]) # in the dynamo region
            mean_density *= 1000 # g/cm^3 to kg/m^3
            T = np.power(10, p.logT)
            P = np.power(10, p.logP)
            delta = p.dlnRho_dlnT_const_Pgas
            
            idx = np.argmin(np.abs(r - Rdyn_start))
            q = convective_flux(p.cp[idx], T[idx], density[idx], 
                                p.conv_vel[idx], P[idx], delta[idx])
            q0 = q * Rdyn_start**2 / Rdyn_end**2
            q0 *= 1e-3 # mW/m^2 to W/m^2
            #print(q0, q)

            # Calc B | B^2 = 2mu_0 c <rho>^1/3 q(r)^2/3  | mu_0 cgs is 1
            B_dyn_squared_SI = 2  * c.mu_0_SI * c.porp_c * mean_density**(1/3)\
                                * q0**(2/3) # T
            B_dyn_SI = np.sqrt(B_dyn_squared_SI)
            B_dyn = B_dyn_SI * 10_000 # Tesla to Gauss
            
            #print(B_dyn)
            # Get Dipole
            dynamo = (r[0] - Rdyn_end) / r[0]
            B_dip =  B_dyn * np.power( 1 - dynamo, 3) / np.sqrt(2)
            
            # Save
            hold(age, B_dyn, B_dip)

        apothikh.append(hold)
    return apothikh


def plotter(names, names2, cols, labels, bigfirst = True):
    
    # Specify Palettes
    if cols == 1:
        colors = [c.AEK, 'k', 'maroon']
    elif cols == 2:
        colors = [c.darkb, c.cyan, c.prasinaki, c.yellow, c.kroki, c.reddish]
    elif cols == 3:
        colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c99]

    # Makes the calculations
    planets1 = medB_doer(names) 
    planets2 = ReyB_doer(names2)
    planets = planets1 + planets2
    fig, axs = plt.subplots(1,2, tight_layout = True, sharex = True,
                           figsize = (6,4))
    
    custom_lines = []
    for planet, color, i in zip(planets, colors, range(len(planets))):
        if i == 0 and bigfirst:
            axs[0].plot(planet.age, planet.Bdyn, color = color, lw=5)
            axs[1].plot(planet.age, planet.Bdip, color = color, lw=5)
                
        axs[0].scatter(planet.age, planet.Bdyn, color = color, s =10)
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
    
    # axs[0].set_ylim(300, 600)
    # axs[1].set_ylim(300, 550)
    axs[0].set_yscale('log')
    axs[1].set_yscale('log')

    fig.suptitle('Scaling vs F=1', 
                  fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)

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
        
    if len(planets) == 8:
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
    
    # name = 'jupenvtest1'
    # name2 = 'jupenvtest2'
    # name3 = 'jup17'
    # name4 = 'jupenvtest3'
    # name5 = 'jupenvtest4'
    # name6 = 'jupenvtest5'
    # name7 = 'jupenvtest6'
    # name8 = 'jupenvtest7'
    # labels = ['30', '40', '50', '60', '70', '80', '90', '95']
    # plotter([name, name2, name3, name4, name5, name6, name7, name8], 3, labels, False)
    
    name = 'jup17profiletest'
    name2 = 'jup17profiletest'
    labels = ['F=1', 'Scaling']
    plotter([name], [name2], 1, labels, False)
