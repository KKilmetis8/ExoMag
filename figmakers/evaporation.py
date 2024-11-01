#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 13:25:59 2024

@author: konstantinos

Figure 7
"""

import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
# 
import src.prelude as c
from src.Bfield.reynolds import rey_mag, profile_sorter
from src.Bfield.hardB import hardB_doer_single


class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.start = []
        self.end = []
        self.radius = []
        self.relative = []
        self.Bdip = []
    def __call__(self, age, start_h, end_h, radius_h, bdip):
        self.age.append(age)
        self.start.append(start_h)
        self.end.append(end_h)
        self.radius.append(radius_h)
        self.Bdip.append(bdip)
        
def dynsize_doer(names):
    apothikh = []
    for name in names:
        # Instanciate Holder Object
        hold = apothicarios(name)
        
        # Get profiles
        p_path = name #'data/' + name
        profiles = os.popen('ls ' + p_path + '/profile*.data').read()
        profiles = list(profiles.split("\n"))
        profiles.pop() # Remove last
        profiles = profile_sorter(profiles) # Guess what that does
    
        # Do the thing
        for profile in profiles:
            p = mr.MesaData(profile)
            r, reynolds_mag_number, age = rey_mag(p)
            r *=  c.Rsol / c.Rearth
            R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
            Bdyn, B_dip, age = hardB_doer_single(profile)
            try:
                test = np.abs(R_dynamo_active[0] - R_dynamo_active[-1]) 
            except IndexError:
                hold(age,0,0,0,0)
                continue
            hold(age,
                 R_dynamo_active[-1], R_dynamo_active[0],
                 r[0], B_dip)
            
        # Keep em
        apothikh.append(hold)
        print('---')
    return apothikh

#%%
pre = 'data'
name1 = f'{pre}/m17_env0.06_zero_a0.1_s8'
name2 = f'{pre}/m17_env0.06_EL_a0.1_s8'
name3 = f'{pre}/m17_env0.06_HD_a0.1_s8'
names = [name1, name2, name3,]
labels = ['zero', 'EL', 'HD']
planets = dynsize_doer(names) 
#%%
plt.rcParams['text.usetex'] = True
colors = ['plum',c.yellow, c.cyan]

alphas = [0.7, 0.7, 0.6]
#alphas = [0.7, 0.5, 0.3,]
#alphas = [0.3, 0.3, 0.3]

hatches = ['|','-', '']
fig, axs = plt.subplots(1,2, tight_layout = True, sharex = False,
                       figsize = (7,4))
custom_lines = []
for planet, color, label, alpha, hatch in zip(planets, colors, labels, alphas, hatches):
    stop = np.argmin(planet.end)
    axs[0].plot(planet.age[:stop], planet.radius[:stop], ls = '-.',
              color = color)
    axs[1].plot(planet.age[:stop], planet.Bdip[:stop], ls = '-',
              color = color)
    axs[1].scatter(planet.age[stop-1], planet.Bdip[stop-1], ls = '', marker = 'X',
              color = color, s = 100, ec='k')
    axs[0].fill_between(planet.age[:stop], planet.start[:stop], planet.end[:stop], 
                        color = color, alpha = alpha, label = label,
                        hatch = hatch, edgecolor = 'k')

# # Make nice
axs[0].set_ylabel('Radial Distance[R$_\oplus$]', fontsize = 14)    
axs[1].set_ylabel('B$_\mathrm{dip}^\mathrm{(max)}$ [G]', fontsize = 14)    

#axs[0].grid()
axs[0].set_xlim(300, 3_800)
axs[1].set_xlim(300, 3_850)

axs[0].set_ylim(2.1,4)
axs[1].set_ylim(-1,12)
#fig.suptitle(title,  fontsize = 15, y = 0.98)

fig.text(0.47, -0.02, r' Age [Myr]', 
          fontsize = 15, transform = fig.transFigure)
# axs[0].legend(fontsize )






axs[0].text(2500, 3.7, r'$\mathrm{\mathbf{Planetary}}$ $\mathrm{\mathbf{Surface}}$',
              fontsize = 15, color = 'k', weight = 'bold',
              horizontalalignment='center', rotation = -5)
axs[0].text(2500, 2.3, r'$\mathrm{\mathbf{Dynamo}}$ $\mathrm{\mathbf{Region}}$',
              fontsize = 15, color = 'k', weight = 'bold',
              horizontalalignment='center', rotation = -25)

plt.rcParams['text.usetex'] = False
axs[1].text(3700, 8, 'No evaporation',
              fontsize = 15, color = colors[0], weight = 500,
              horizontalalignment='right')
axs[1].text(3700, 7, 'Energy Limit',
              fontsize = 15, color = colors[1], weight = 500,
              horizontalalignment='right')
axs[1].text(3700, 6, 'Hydrodynamic',
              fontsize = 15, color = colors[2], weight = 500,
              horizontalalignment='right')
plt.savefig('figs/evaporation_small.pdf', format = 'pdf', dpi = 300, bbox_inches='tight')

