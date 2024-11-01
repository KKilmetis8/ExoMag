#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 13:01:20 2024

@author: konstantinos

Figure 3
"""

# Vanilla
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mesa_reader as mr

# Choc
import src.prelude as c
import src.Utilities.planet_grids as grids
from src.Bfield.hardB import hardB_doer_single

class apothicarios:
    def __init__(self):
        self.masses = []
        self.orbsep = []
        self.Bdyn = []
        self.Bdip = []
        self.radii = []
    def __call__(self, m, orbsep, dyn, dip, r):
        self.masses.append(m)
        self.orbsep.append(orbsep)
        self.Bdyn.append(dyn)
        self.Bdip.append(dip)
        self.radii.append(r)
        
ages = [1500, 5500, 10000] # Myrs
orb_seps1 = list(np.arange(0.05, 0.3, 0.01))
orb_seps2 = list(np.arange(0.3, 2, 0.1))
orb_seps = orb_seps1 + orb_seps2

gas_giants_1 = apothicarios()
gas_giants_5 = apothicarios()
gas_giants_10 = apothicarios()
dunk = '/media/konstantinos/Dunkey/mesadata/'
#%% Calc
for age in ages:
    df = pd.read_csv(f'data/specific_ages/{age}.txt', header = None, 
                     delimiter = '\s', names = ['Planet', 'Profile'])
    for planet, profile_no in zip(df['Planet'], df['Profile']):
        modelname = planet
        underscores = [ i for i, x in enumerate(modelname) if x == '_']
        m = int( planet[ 1: underscores[0] ])
        env = float(planet[underscores[0]+4:underscores[1]])
        if m != 317 or env != 0.94:
            continue
        escape = planet[ underscores[1]+1:underscores[2] ]
        a = float( planet[ underscores[2]+2:underscores[3]])

        if int(profile_no) > 0:
            try:
                B_dyn, B_dip, age_this = hardB_doer_single(f'{dunk}data/{modelname}/profile{profile_no}.data')
                if np.abs(age-age_this) > 500:
                    continue
                p = mr.MesaData(f'{dunk}data/{modelname}/profile{profile_no}.data')
                r = np.power(10, p.logR)[2] *  c.Rsol / c.Rearth

            except (ValueError,FileNotFoundError) as e:
                print(f'{planet} not found')
                continue
            if  B_dyn>0:
                if age == ages[1]:
                    gas_giants_5(m, a, B_dyn, B_dip, r)
                elif age == ages[2]:
                    gas_giants_10(m, a, B_dyn, B_dip, r)
                elif age == ages[0]:
                    gas_giants_1(m, a, B_dyn, B_dip, r)

#%% Plot
# ---------------- Dynamo ----------------
fig, ax = plt.subplots(1,2, figsize = (6,4), sharex = True, sharey= True) 
a1, Bdyn1, Bdip1 = zip(*sorted(zip(gas_giants_1.orbsep, 
                                   gas_giants_1.Bdyn, 
                                   gas_giants_1.Bdip)))
a5, Bdyn5, Bdip5 = zip(*sorted(zip(gas_giants_5.orbsep, 
                                   gas_giants_5.Bdyn,
                                   gas_giants_5.Bdip)))
a10, Bdyn10, Bdip10 = zip(*sorted(zip(gas_giants_10.orbsep,
                                      gas_giants_10.Bdyn,
                                      gas_giants_10.Bdip)))
ax[0].plot(a10, Bdyn10, 
              c = 'maroon', ls = '-', lw = 1, markeredgecolor = 'k',
              marker = 'h', markersize = 7)

ax[0].plot(a5, Bdyn5, 
              c = 'k', ls = '-', lw = 1, 
              markeredgecolor = 'k', marker = 'h', markersize = 7)

ax[0].plot(a1, Bdyn1, 
              c = c.AEK, ls = '-', lw = 1, 
              markeredgecolor = 'k', marker = 'h', markersize = 7)
 

# ---------------- Dipole ----------------
ax[1].plot(a10, Bdip10, 
              c = 'maroon', ls = '-', lw = 1, markeredgecolor = 'k',
              marker = 'h', markersize = 7)

ax[1].plot(a5, Bdip5, 
              c = 'k', ls = '-', lw = 1, 
              markeredgecolor = 'k', marker = 'h', markersize = 7)

ax[1].plot(a1, Bdip1, 
              c = c.AEK, ls = '-', lw = 1, 
              markeredgecolor = 'k', marker = 'h', markersize = 7)
 

ax2 = ax[0].twiny()
ax2.plot(np.array(a10) * c.AU/c.Rsol, [-1] * len(a10))
ax2.set_xlabel(r'Orbital Separation  [$R_\odot$]', fontsize = 15, labelpad=10)
ax3 = ax[1].twiny()
ax3.plot(np.array(a10) * c.AU/c.Rsol, [-1] * len(a10))
ax3.set_xlabel(r'Orbital Separation  [$R_\odot$]', fontsize = 15, labelpad=10)


# Text
ax[0].text(1, 235, r'$\mathbf{1.5}$ $\mathbf{Gyr}$',
              fontsize = 15, color = c.AEK, weight = 'bold',
              horizontalalignment='center')
ax[0].text(1, 175, r'$\mathbf{5.5}$ $\mathbf{Gyr}$',
              fontsize = 15, color = 'k', weight = 'bold',
              horizontalalignment='center')
ax[0].text(1, 110, r'$\mathbf{10}$ $\mathbf{Gyr}$',
              fontsize = 15, color = 'maroon',
              weight = 'bold',
              horizontalalignment='center')

# Labels & Grid
ax[0].grid()
ax[0].set_xlabel(r'Orbital Separation  [au]', fontsize = 14)
ax[0].set_ylabel(r'$B_\mathrm{dyn}^\mathrm{(max)}$ [G]', fontsize = 14)
ax[1].grid()
ax[1].set_xlabel(r'Orbital Separation  [au]', fontsize = 14)
ax[1].set_ylabel(r'$B_\mathrm{dip}^\mathrm{(max)}$ [G]', fontsize = 14)
ax[0].set_xlim(0,2)
ax[0].set_ylim(40,280)
# Fig
# fig.suptitle(r'5.5 Gyrs', fontsize = 16, y = 0.95)

# Inset Dynamo ----------------------------------------------------------------
# axi = ax[0].inset_axes( [.35, .08, .5, .35], xlim=(0.03, 0.16), ylim=[80,275],
#                        xticks = [0.04, 0.09, 0.14], yticks = [125, 175, 225],)
# axi.set_facecolor('snow')
# ax[0].indicate_inset_zoom(axi, edgecolor="black")
# start = 0
# stop = 11
# axi.plot(a10[start:stop], Bdyn10[start:stop], 
#               c = 'maroon', ls = '-', lw = 1, markeredgecolor = 'k',
#               marker = 'h', markersize = 7)
# axi.plot(a5[start:stop], Bdyn5[start:stop], 
#               c = 'k', ls = '-', lw = 1, markeredgecolor = 'k',
#               marker = 'h', markersize = 7)
# axi.plot(a1[start:stop], Bdyn1[start:stop], 
#               c = c.AEK, ls = '-', lw = 1, markeredgecolor = 'k',
#               marker = 'h', markersize = 7)
# axi.set_xticklabels([r'$\mathrm{\mathbf{0.04}}$',
#                      r'$\mathrm{\mathbf{0.09}}$',
#                      r'$\mathrm{\mathbf{0.14}}$'], weight = 'bold', fontsize = 8)
# axi.set_yticklabels([r'$\mathrm{\mathbf{125}}$',
#                      r'$\mathrm{\mathbf{175}}$',
#                      r'$\mathrm{\mathbf{225}}$'], weight = 'bold', fontsize = 8)

# Inset Dipole ---------------------------------------------------------------
axi = ax[1].inset_axes( [.55, .65, .45, .35], xlim=(0.03, 0.16), ylim=[45,200],
                       xticks = [0.04, 0.09, 0.14], yticks = [75, 125, 175],)
axi.set_facecolor('snow')
ax[1].indicate_inset_zoom(axi, edgecolor="#3d3c3c")
start = 0
stop = 11
stop10 = 13
axi.plot(a10[start:stop10], Bdip10[start:stop10], 
              c = 'maroon', ls = '-', lw = 1, markeredgecolor = 'k',
              marker = 'h', markersize = 7)
axi.plot(a5[start:stop], Bdip5[start:stop], 
              c = 'k', ls = '-', lw = 1, markeredgecolor = 'k',
              marker = 'h', markersize = 7)
axi.plot(a1[start:stop], Bdip1[start:stop], 
              c = c.AEK, ls = '-', lw = 1, markeredgecolor = 'k',
              marker = 'h', markersize = 7)
axi.set_xticklabels([r'$\mathrm{\mathbf{0.04}}$',
                     r'$\mathrm{\mathbf{0.09}}$',
                     r'$\mathrm{\mathbf{0.14}}$'], weight = 'bold', fontsize = 8)
axi.set_yticklabels([r'$\mathrm{\mathbf{75}}$',
                     r'$\mathrm{\mathbf{125}}$',
                     r'$\mathrm{\mathbf{175}}$'], weight = 'bold', fontsize = 8)

plt.savefig('figs/orbsep_ages.pdf', format = 'pdf', dpi = 300, bbox_inches='tight')
