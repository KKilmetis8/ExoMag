#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 19:55:13 2024

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
        self.Bdyn = []
        self.Bdip = []
        self.F = []
        self.radius = []

    def __call__(self, age_h, Bdyn_h, Bdip_h, F_h, radius_h):
        self.age.append(age_h)
        self.Bdyn.append(Bdyn_h)
        self.Bdip.append(Bdip_h)
        self.F.append(F_h)
        self.radius.append(radius_h)
        
def hardB_doer_single(profile):
    p = mr.MesaData(profile)
    r = np.power(10, p.logR)
    R_dynamo_active, _, age = dynamo_region(p)
    try:
        Rdyn_end = R_dynamo_active[0] # it's the wrong way round
        # R_dynamo_active *= c.Rsol / c.Rearth
    except IndexError:
        return -100, -100
    B_dyn, F = hardB(p, R_dynamo_active)
    # Get Dipole
    dynamo = (r[0] - Rdyn_end) / r[0]
    B_dip =  B_dyn * np.power( 1 - dynamo, 3) / np.sqrt(2)
    return B_dyn, B_dip, age

def hardB_doer(names):
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
                # R_dynamo_active *= c.Rsol / c.Rearth
            except IndexError:
                # Save
                hold(age, hold.Bdyn[-1], hold.Bdip[-1], 0, 0)
                continue  
            B_dyn, F = hardB(p, R_dynamo_active)
            # Get Dipole
            dynamo = (r[0] - Rdyn_end) / r[0]
            B_dip =  B_dyn * np.power( 1 - dynamo, 3) / np.sqrt(2)
            # Save
            hold(age, B_dyn, B_dip, F, r[0])
        apothikh.append(hold)
    return apothikh
    
def plotter(names, cols, labels, title, bigfirst = False):
    
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
    elif cols == 'blues':
        colors = [c.n6, c.n5, c.n4, c.n3, c.n2, c.n1]

    # Makes the calculations
    planets = hardB_doer(names) 
    fig, axs = plt.subplots(1,2, tight_layout = True, sharex = True,
                           figsize = (5,4))
    custom_lines = []
    for planet, color, i in zip(planets, colors, range(len(planets))):
        if i == 0 and bigfirst:
            #axs[0].plot(planet.age, planet.F, color = color, lw=10, zorder = 2)
            axs[0].plot(planet.age, planet.Bdyn, color = color, lw=10, zorder = 2)
            axs[1].plot(planet.age, planet.Bdip, color = color, lw=10, zorder = 2)
            
        if i == 1 and bigfirst:
            #axs[0].plot(planet.age, planet.F, color = color, lw=5, zorder = 2)
            axs[0].plot(planet.age, planet.Bdyn, color = color, lw=5, zorder = 2)
            axs[1].plot(planet.age, planet.Bdip, color = color, lw=5, zorder = 2)
        
        #axs[0].plot(planet.age, planet.F, color = color)
        axs[0].plot(planet.age, planet.Bdyn, color = color)
        axs[1].plot(planet.age, planet.Bdip, color = color)
    
        # Legend
        custom_lines.append( Line2D([0], [0], color = colors[i], 
                                    label = labels[i]))
        
    # Make nice
    #axs[0].set_ylabel('Efficiency Factor', fontsize = 14)
    axs[0].set_ylabel('Dynamo [G]', fontsize = 14)
    axs[1].set_ylabel('Dipole [G]', fontsize = 14)
    
    #axs[0].grid()
    axs[0].grid()
    axs[1].grid()
        
    axs[1].set_xlim(500, 10_000)
    #axs[0].set_yscale('log')
    if 'm317' in names[0]:
        #axs[0].set_ylim(1_000, 6_000)
        axs[0].set_ylim(50, 350)
        axs[1].set_ylim(50, 350)
    elif 'm17' in names[0]:
        # axs[0].set_ylim(1E-12,1E1)
        axs[0].set_ylim(0, 40)
        axs[1].set_ylim(0, 40)
    
    #axs[0].set_yscale('log')
    #axs[1].set_yscale('log')

    # fig.suptitle(title, 
    #               fontsize = 18, y = 0.98)
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
        
    fig.legend(custom_lines, labels,
            fontsize =  9, ncols = len(planets), alignment = 'left', # Lawful Neutral
            bbox_to_anchor=(0.94, -0.03), bbox_transform = fig.transFigure,)
#%%    
kind = 'nepenv_zero'
if __name__ == '__main__':
    if kind == 'jup-kai-nep':
        fenvs_jup = ['88', '9', '92', '94', '96']
        fenvs_nep = ['02', '04','06', '08', '1']
        labels_jup = ['88%', '2%', '90%', '4%', '92%', '6%' '94%', '8%', '96%',
                      '98%']
        
        names = []
        for fj,fn in zip(fenvs_jup, fenvs_nep):
            name_j = f'm317_env{fj}_zero_a01_s8'
            name_n = f'm317_env{fn}_zero_a01_s8'
            names.append(name_j)
            names.append(name_n)
            
            
    if kind == 'jupenv_zero_orb':
        name3 = 'm317_env94_zero_a01_s8'
        name4 = 'm317_e88_zero_a01_s8'
        name5 = 'm317_e91_zero_a01_s8'
        name6 = 'm317_e94_zero_a01_s8'
        name7 = 'm317_e97_zero_a01_s8'
        #name8 = 'jup_e97_zero'
        #name9 = 'jup3_e98_zero'
        names = [name3, name4, name5, name6, name7]#, name4, name9]
        labels = ['85', '90', '92', '94', '96',]#, '95', '98']
        plotter(names, 4, labels, r'Jupiter $317 M_\oplus$, $\alpha$=0.1 AU')
    if kind == 'jupenv_zero':
        name3 = 'm317_e85_zero_a01_s8'
        name4 = 'm317_e88_zero_a01_s8'
        name5 = 'm317_e91_zero_a01_s8'
        name6 = 'm317_e94_zero_a01_s8'
        name7 = 'm317_e97_zero_a01_s8'
        #name8 = 'jup_e97_zero'
        #name9 = 'jup3_e98_zero'
        names = [name3, name4, name5, name6, name7]#, name4, name9]
        labels = ['85', '90', '92', '94', '96',]#, '95', '98']
        plotter(names, 4, labels, r'Jupiter $317 M_\oplus$, $\alpha$=0.1 AU')
    if kind == 'nepenv_zero':
        name1 = 'm17_env0.02_zero_a0.1_s8'
        name2 = 'm17_env0.04_zero_a0.1_s8'
        name3 = 'm17_env0.06_zero_a0.1_s8'
        name4 = 'm17_env0.08_zero_a0.1_s8'
        name5 = 'm17_env0.1_zero_a0.1_s8'
        name6 = 'm17_env0.12_zero_a0.1_s8'
        names = [name1, name2, name3, name4, name5, name6]
        labels = ['2', '4', '6', '8', '10', '12',]
        plotter(names, 4, labels, 'Neptunes with Diff. Envelopes')

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
