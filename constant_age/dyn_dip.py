#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 16:15:30 2024

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mesa_reader as mr
import src.prelude as c
from src.Bfield.hardB import hardB_doer_single

class apothicarios:
    def __init__(self):
        self.masses = []
        self.envs = []
        self.entropies = []
        self.Bdyn = []
        self.Bdip = []
        self.radii = []
    def __call__(self, m, fenv, s, dyn, dip, r):
        self.masses.append(m)
        self.envs.append(fenv)
        self.entropies.append(s)
        self.Bdyn.append(dyn)
        self.Bdip.append(dip)
        self.radii.append(r)

        
age = 5500 # Myrs
df = pd.read_csv(f'data/specific_ages/{age}.txt', header = None, 
                              delimiter = '\s', names = ['Planet', 'Profile'])
ms = [8, 14, 17, 20, 25, 30, 40, 50, 60, 75, 150, 200, 250, 300, 350, 400]
fenvs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 60, 70, 80, 85, 90, 92, 94, 95, 96, 98]
entropies = [7, 8, 9]
gas_giants = apothicarios()
neptunes = apothicarios()
#%% Calc
for m in ms:
    for fenv in fenvs:
        for s in entropies:
            modelname = f'm{m}_e{fenv}_zero_a01_s{s}'
            if modelname in list(df['Planet']):
                profile = df.loc[df['Planet'] == modelname]
                profile_no = profile.iloc[0,1]
                try:
                    B_dyn, B_dip = hardB_doer_single(f'data/{modelname}/profile{profile_no}.data')
                    p = mr.MesaData(f'data/{modelname}/profile{profile_no}.data')
                    r = np.power(10, p.logR)[2] *  c.Rsol / c.Rearth
                except ValueError or FileNotFoundError:
                    continue
                if fenv > 29: 
                    gas_giants(m, fenv, s, B_dyn, B_dip, r)
                else:
                    neptunes(m, fenv, s, B_dyn, B_dip, r)
                
#%% Plot

fig, ax = plt.subplots(1,2, figsize = (6,4), sharex = True, sharey= True) 
ax[0].scatter(gas_giants.masses, gas_giants.Bdyn, ec = 'k', zorder = 3,
                 vmin = 60, vmax = 100, c = gas_giants.envs, cmap = 'Oranges')
ax[0].scatter(neptunes.masses, neptunes.Bdyn, ec = 'k', zorder = 3,
                 vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')
im1 = ax[1].scatter(gas_giants.masses, gas_giants.Bdip, ec = 'k', zorder = 3,
                 vmin = 60, vmax = 100, c = gas_giants.envs, cmap = 'Oranges')
im2 = ax[1].scatter(neptunes.masses, neptunes.Bdip, ec = 'k', zorder = 3,
                 vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')

# Limits and scale
ax[0].set_ylim(0,160)
ax[0].set_xscale('log')

# Labels & Grid
ax[0].grid()
ax[0].set_xlabel(r'Planet Mass [M$_\oplus$]', fontsize = 14)
ax[0].set_ylabel(r'Dynamo [G]', fontsize = 14)
ax[1].grid()
ax[1].set_xlabel(r'Planet Mass [M$_\oplus$]', fontsize = 14)
ax[1].set_ylabel(r'Dipole [G]', fontsize = 14)
# Fig
fig.colorbar(im1, cax = fig.add_axes([0.92, 0.12, 0.02, 0.36]),)
fig.colorbar(im2, cax = fig.add_axes([0.92, 0.52, 0.02, 0.36]),)
fig.suptitle(r'5.5 Gyrs, $\alpha$ 0.1 AU', fontsize = 16, y = 0.95)
fig.text(0.98, 0.35, 'Envelope Fraction', rotation = 90, fontsize = 14)
#%%
fig, ax = plt.subplots(1,2, figsize = (6,4), sharex = True, sharey= True) 
ax[0].scatter(gas_giants.radii, gas_giants.Bdyn, ec = 'k', zorder = 3,
                 vmin = 60, vmax = 100, c = gas_giants.envs, cmap = 'Oranges')
ax[0].scatter(neptunes.radii, neptunes.Bdyn, ec = 'k', zorder = 3,
                 vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')
im1 = ax[1].scatter(gas_giants.radii, gas_giants.Bdip, ec = 'k', zorder = 3,
                 vmin = 60, vmax = 100, c = gas_giants.envs, cmap = 'Oranges')
im2 = ax[1].scatter(neptunes.radii, neptunes.Bdip, ec = 'k', zorder = 3,
                 vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')

# Limits and scale
# ax[0].set_ylim(0,160)
ax[0].set_xlim(2,13)
#ax[0].set_xscale('log')

# Labels & Grid
ax[0].grid()
ax[0].set_xlabel(r'Planet Radius [R$_\oplus$]', fontsize = 14)
ax[0].set_ylabel(r'Dynamo [G]', fontsize = 14)
ax[1].grid()
ax[1].set_xlabel(r'Planet Radius [R$_\oplus$]', fontsize = 14)
ax[1].set_ylabel(r'Dipole [G]', fontsize = 14)
# Fig
fig.colorbar(im1, cax = fig.add_axes([0.92, 0.12, 0.02, 0.36]),)
fig.colorbar(im2, cax = fig.add_axes([0.92, 0.52, 0.02, 0.36]),)
fig.suptitle(r'5.5 Gyrs, $\alpha$ 0.1 AU', fontsize = 16, y = 0.95)
fig.text(0.98, 0.35, 'Envelope Fraction', rotation = 90, fontsize = 14)