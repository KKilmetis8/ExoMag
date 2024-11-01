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
from src.Bfield.hardB import hardB_doer

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

colors = [c.c91, c.c92, c.c93, c.c95, c.c97, c.c99,]

# Makes the calculations
jups = hardB_doer(names_j) 
neps = hardB_doer(names_n) 

#%% Pepper it

stop = [ np.argmin(nep.Bdyn)+1 for nep in neps]
#%%
import src.prelude as c

fig, axs = plt.subplots(2,2, tight_layout = False, sharex = True,
                       figsize = (6,6))
custom_lines = []

for i in range(len(jups)):
    #axs[0].plot(planet.age, planet.F, color = color)
    axs[0,0].plot(jups[i].age, jups[i].Bdyn, 
                color = colors[i], ls = '-', label = labels_j[i])
    axs[1,0].plot(neps[i].age[:stop[i]], neps[i].Bdyn[:stop[i]], 
                color = colors[i], ls = '-.', label = labels_n[i])
    axs[0,1].plot(jups[i].age, jups[i].Bdip, 
                color = colors[i], ls = '-', label = labels_j[i])
    axs[1,1].plot(neps[i].age[:stop[i]], neps[i].Bdip[:stop[i]], 
                color = colors[i], ls = '-.', label = labels_n[i])
    
    # X point
    axs[1,0].plot(neps[i].age[stop[i] -1], neps[i].Bdyn[stop[i]], 
                color = colors[i],
                marker = 'X', markeredgecolor = 'k', markersize = 10)
    axs[1,1].plot(neps[i].age[stop[i] -1], neps[i].Bdip[stop[i]], 
                color = colors[i],
                marker = 'X', markeredgecolor = 'k', markersize = 10)
# Make nice
#axs[0].set_ylabel('Efficiency Factor', fontsize = 14)
axs[0,0].set_ylabel(r'$B_\mathrm{dyn}^\mathrm{(max)}$ [G]', fontsize = 14)
axs[0,1].set_ylabel(r'$B_\mathrm{dip}^\mathrm{(max)}$ [G]', fontsize = 14)
axs[1,0].set_ylabel(r'$B_\mathrm{dyn}^\mathrm{(max)}$ [G]', fontsize = 14)
axs[1,1].set_ylabel(r'$B_\mathrm{dip}^\mathrm{(max)}$ [G]', fontsize = 14)


#axs[0].grid()
# axs[0,0].grid()
# axs[0,1].grid()
# axs[1,0].grid()
# axs[1,1].grid()
    
axs[0,0].set_xlim(300, 8000)

axs[0,0].set_ylim(120, 400)
axs[0,1].set_ylim(80, 250)
axs[1,0].set_ylim(-1, 37)
axs[1,1].set_ylim(-0.4, 12)

# axs[0,0].text(5000, 300, r'$\mathbf{1}$ $\mathbf{M_J}$',
#               fontsize = 30, color = 'k', weight = 'bold',
#               horizontalalignment='center')
# axs[1,0].text(5000, 25, r'$\mathbf{17}$ $\mathbf{M_\oplus}$',
#               fontsize = 30, color = 'k', weight = 'bold',
#               horizontalalignment='center',)

# fig.suptitle(title, 
#               fontsize = 18, y = 0.98)
fig.text(0.47, -0.02, r' Age [Myr]', 
          fontsize = 14, transform = fig.transFigure)

# lines = axs[0,0].get_lines()
# labels = [line.get_label() for line in lines]
# fig.legend(lines[::-1], labels[::-1],
#         fontsize = 12 , ncols = 1, alignment = 'center',# Lawful Neutral
#         bbox_to_anchor=(0.95, 0.97), bbox_transform = fig.transFigure)
# lines = axs[1,1].get_lines()
# labels = [line.get_label() for line in lines]
# fig.legend(lines[::-1], labels[::-1],
#         fontsize = 12 , ncols = 1, alignment = 'center',# Lawful Neutral
#         bbox_to_anchor=(0.95, 0.49), bbox_transform = fig.transFigure)

# plt.savefig('figs/fig2.pdf', format = 'pdf', dpi = 300, bbox_inches='tight')

# #%% 
# fig, axs = plt.subplots(1,2, tight_layout = False, sharex = True,
#                        figsize = (6,4))

# for i in range(0,5):
#     axs[0].plot(neps[i].age[:stop[i]], neps[i].Bdyn[:stop[i]], 
#                 color = colors[i], ls = '-', label = labels_n[i])
#     axs[0].plot(neps[i].age[stop[i] -1], neps[i].Bdyn[stop[i]], 
#                 color = colors[i],
#                 marker = 'X', markeredgecolor = 'k', markersize = 10)
#     axs[1].plot(neps[i].age[:stop[i]], neps[i].Bdip[:stop[i]], 
#                 color = colors[i], ls = '-', label = labels_n[i])
#     axs[1].plot(neps[i].age[stop[i] -1], neps[i].Bdip[stop[i]], 
#                 color = colors[i],
#                 marker = 'X', markeredgecolor = 'k', markersize = 10)

# axs[0].set_xlim(300,5200)
# axs[0].set_ylim(-2,30)
# axs[1].set_ylim(-2,30)

# axs[0].set_ylabel(r'$B_\mathrm{dyn}^\mathrm{(max)}$ [G]', fontsize = 14)
# axs[1].set_ylabel(r'$B_\mathrm{dip}^\mathrm{(max)}$ [G]', fontsize = 14)
# fig.text(0.47, -0.02, r' Age [Myr]', 
#           fontsize = 15, transform = fig.transFigure)
