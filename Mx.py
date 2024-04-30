#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:20:13 2024

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:10:44 2024

@author: konstantinos
"""

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
    R_dynamo_active, _, _ = dynamo_region(p)
    try:
        Rdyn_end = R_dynamo_active[0] # it's the wrong way round
        # R_dynamo_active *= c.Rsol / c.Rearth
    except IndexError:
        return -100, -100
    B_dyn, F = hardB(p, R_dynamo_active)
    # Get Dipole
    dynamo = (r[0] - Rdyn_end) / r[0]
    B_dip =  B_dyn * np.power( 1 - dynamo, 3) / np.sqrt(2)
    return B_dyn, B_dip

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

#%%
fenvs_jup = ['88', '9', '92', '94', '96']
fenvs_nep = ['02', '04','06', '08', '1']
labels_j = ['88\%', '90\%', '92\%', '94\%', '96\%']
labels_n = ['2\%', '4\%', '6\%', '8\%', '10\%']
names_j = []
names_n = []
for fj,fn in zip(fenvs_jup, fenvs_nep):
    fj = '0.' + fj
    fn = '0.' + fn
    name_j = f'm317_env{fj}_zero_a0.1_s8'
    name_n = f'm17_env{fn}_zero_a0.1_s8'
    names_j.append(name_j)
    names_n.append(name_n)
title = ''
# plotter(names_j, names_n, labels_j, labels_n, title)

colors = [c.c91, c.c92, c.c93, c.c95, c.c97, c.c99,]

# Makes the calculations
jups = hardB_doer(names_j) 
neps = hardB_doer(names_n) 
#%%
fig, axs = plt.subplots(2,2, tight_layout = True, sharex = True,
                       figsize = (6,6))
custom_lines = []

for i in range(len(jups)):
    dynflux_jup = np.array(jups[i].Bdyn) * (np.array(jups[i].radius) * c.Rsol)**2
    dipflux_jup = np.array(jups[i].Bdip) * (np.array(jups[i].radius) * c.Rsol)**2
    dynflux_nep = np.array(neps[i].Bdyn) * (np.array(neps[i].radius) * c.Rsol)**2
    dipflux_nep = np.array(neps[i].Bdip) * (np.array(neps[i].radius) * c.Rsol)**2
    
    #axs[0].plot(planet.age, planet.F, color = color)
    axs[0,0].plot(jups[i].age, dynflux_jup, 
                color = colors[i], ls = '-', label = labels_j[i])
    axs[1,0].plot(neps[i].age, dynflux_nep, 
                color = colors[i], ls = '-.', label = labels_n[i])
    axs[0,1].plot(jups[i].age, dipflux_jup, 
                color = colors[i], ls = '-', label = labels_j[i])
    axs[1,1].plot(neps[i].age, dipflux_nep, 
                color = colors[i], ls = '-.', label = labels_n[i])
    
# Make nice
#axs[0].set_ylabel('Efficiency Factor', fontsize = 14)
axs[0,0].set_ylabel(r'$B_\mathrm{dyn}R_p^2$  [Mx]', fontsize = 14)
axs[1,0].set_ylabel(r'$B_\mathrm{dip}R_p^2$ [Mx]', fontsize = 14)

#axs[0].grid()
axs[0,0].grid()
axs[0,1].grid()
axs[1,0].grid()
axs[1,1].grid()
#axs[:,:].set_yscale('log')
axs[0,0].set_xlim(500, 8000)

axs[0,0].set_ylim(1e21, 3e22)
axs[0,1].set_ylim(1e21, 3e22)
axs[1,0].set_ylim(1e15, 3e20)
axs[1,1].set_ylim(1e15, 3e20)


#axs[0].set_yscale('log')
#axs[1].set_yscale('log')

# fig.suptitle(title, 
#               fontsize = 18, y = 0.98)
fig.text(0.47, -0.02, r' Age [Myr]', 
          fontsize = 15, transform = fig.transFigure)

lines = axs[0,0].get_lines()
labels = [line.get_label() for line in lines]
fig.legend(lines, labels,
        fontsize = 12 , ncols = 1, alignment = 'center',# Lawful Neutral
        bbox_to_anchor=(0.95, 0.95), bbox_transform = fig.transFigure)
lines = axs[1,1].get_lines()
labels = [line.get_label() for line in lines]
fig.legend(lines, labels,
        fontsize = 12 , ncols = 1, alignment = 'center',# Lawful Neutral
        bbox_to_anchor=(0.95, 0.47), bbox_transform = fig.transFigure)