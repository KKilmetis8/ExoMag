#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 18:02:44 2024

@author: konstantinos

Figure 4
"""

# Vanilla
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mesa_reader as mr
from tqdm import tqdm

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
        
age = 1500 # Myrs
df = pd.read_csv(f'data/specific_ages/{age}.txt', header = None, 
                              delimiter = '\s', names = ['Planet', 'Profile'])
orb_seps1 = list(np.arange(0.05, 0.3, 0.01))
orb_seps2 = list(np.arange(0.3, 2, 0.1))
orb_seps = orb_seps1 + orb_seps2
ms = [317, 95, 17, 50]
fenvs =  [0.94, 0.9, 0.06, 0.7] # grids.j_envs
jup = apothicarios()
nep = apothicarios()
supnep = apothicarios()
sat = apothicarios()

#%% Calc
for planet, profile_no in tqdm(zip(df['Planet'], df['Profile'])):
    modelname = planet
    underscores = [ i for i, x in enumerate(modelname) if x == '_']
    m = int( planet[ 1: underscores[0] ])
    escape = planet[ underscores[1]+1:underscores[2] ]
    a = float( planet[ underscores[2]+2:underscores[3]])
    env = float(planet[underscores[0]+4:underscores[1]])

    try:
        profile_no = int(profile_no)
    except ValueError:
        continue
    
    if profile_no > 0 and escape == 'zero':
        try:
            B_dyn, B_dip, age_this = hardB_doer_single(f'data/{age}profile/{modelname}_profile{profile_no}.data')
            p = mr.MesaData(f'data/{age}profile/{modelname}_profile{profile_no}.data')
        except (ValueError, FileNotFoundError) as e:
            try:
                B_dyn, B_dip, age_this = hardB_doer_single(f'data/{age}profile/{modelname}_{profile_no}.data')
                p = mr.MesaData(f'data/{age}profile/{modelname}_{profile_no}.data')
            except (ValueError, FileNotFoundError) as e:
                continue
        if np.abs(age-age_this) > 500:
            print(modelname)
            continue
        r = np.power(10, p.logR)[2] *  c.Rsol / c.Rearth
        if B_dyn > 0:
            if m == ms[0] and env == fenvs[0]:
                jup(m, a, B_dyn, B_dip, r)
            if m == ms[1]:
                sat(m, a, B_dyn, B_dip, r)
            if m == ms[3]:
                supnep(m, a, B_dyn, B_dip, r)
            if m == ms[2] and env == fenvs[2]:
                nep(m, a, B_dyn, B_dip, r)

#%% Plot
jup_orbsep_plot , jup_Bdip_plot = zip(*sorted(zip(jup.orbsep, jup.Bdip)))
sat_orbsep_plot , sat_Bdip_plot = zip(*sorted(zip(sat.orbsep, sat.Bdip)))
nep_orbsep_plot , nep_Bdip_plot = zip(*sorted(zip(nep.orbsep, nep.Bdip)))
supnep_orbsep_plot , supnep_Bdip_plot = zip(*sorted(zip(supnep.orbsep, supnep.Bdip)))

#%% ---------------- Dynamo ----------------
fig, ax = plt.subplots(1,1, figsize = (4,4), sharex = True, sharey= True) 
ax.plot(jup_orbsep_plot, jup_Bdip_plot, 
              c = 'darkorange', ls = '-', lw = 1, 
              markeredgecolor = 'k', markeredgewidth = 0.5, 
              marker = 'h', markersize = 7)
ax.plot(sat_orbsep_plot, sat_Bdip_plot, 
              c = 'goldenrod', ls = '-', lw = 1,
              markeredgecolor = 'k', markeredgewidth = 0.5, 
              marker = 'h', markersize = 7, )
ax.plot(supnep_orbsep_plot[3:], supnep_Bdip_plot[3:], 
              c = 'forestgreen', ls = '-', lw = 1, 
              markeredgecolor = 'k', markeredgewidth = 0.5, 
              marker = 'h', markersize = 7)
ax.plot(nep_orbsep_plot[3:], nep_Bdip_plot[3:], 
              c = 'royalblue', ls = '-', lw = 1, 
              markeredgecolor = 'k', markeredgewidth = 0.5, 
              marker = 'h', markersize = 7)
ax2 = ax.twiny()
ax2.plot(np.array(jup_orbsep_plot) * c.AU/c.Rsol, [-1] * len(jup_Bdip_plot))
ax2.set_xlabel(r'Orbital Separation  [$R_\odot$]', fontsize = 14, labelpad=10)

# Labels & Grid
ax.grid()
ax.set_xlabel(r'Orbital Separation  [au]', fontsize = 14)
#ax[0].set_ylabel(r'$B_\mathrm{dyn}$ [G]', fontsize = 14)
#ax[1].grid()
#ax[1].set_xlabel(r'Orbital Separation  [AU]', fontsize = 14)
ax.set_ylabel(r'$B_\mathrm{dip}^\mathrm{(max)}$ [G]', fontsize = 14)
# ax[0].set_xlim(0.045,0.1)
# ax[0].set_ylim(1,21)
ax.set_yscale('log')
# Fig
#fig.suptitle(r'1.5 Gyrs', fontsize = 16, y = 0.95)

# Text
ax.text(1.8, 115, r'$\mathrm{\mathbf{1}}$ $\mathrm{\mathbf{M_J}}$',
              fontsize = 15, color = 'darkorange', weight = 'bold',
              horizontalalignment='center')
ax.text(1.75, 55, r'$\mathrm{\mathbf{0.3}}$ $\mathrm{\mathbf{M_J}}$',
              fontsize = 15, color = 'goldenrod', weight = 'bold',
              horizontalalignment='center')
ax.text(1.75, 36, r'$\mathrm{\mathbf{50}}$ $\mathrm{\mathbf{M_\oplus}}$',
              fontsize = 15, color = 'forestgreen', weight = 'bold',
              horizontalalignment='center')
ax.text(1.75, 6.75, r'$\mathrm{\mathbf{17}}$ $\mathrm{\mathbf{M_\oplus}}$',
              fontsize = 15, color = 'royalblue', weight = 'bold',
              horizontalalignment='center')
#%%
plt.savefig('figs/planets_orbsep.pdf', format = 'pdf', dpi = 300, bbox_inches='tight')
