#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 16:15:30 2024

@author: konstantinos
"""
# Vanilla
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mesa_reader as mr
from matplotlib.lines import Line2D

# Choc
import src.prelude as c
import src.Utilities.planet_grids as grids
from src.Bfield.hardB import hardB_doer_single
from src.radio.peaks import peak

class apothicarios:
    def __init__(self):
        self.masses = []
        self.envs = []
        self.entropies = []
        self.Bdyn = []
        self.Bdip = []
        self.radii = []
        self.peak = []
    def __call__(self, m, fenv, s, dyn, dip, r, p):
        self.masses.append(m)
        self.envs.append(fenv)
        self.entropies.append(s)
        self.Bdyn.append(dyn)
        self.Bdip.append(dip)
        self.radii.append(r)
        self.peak.append(p)

        
age = 5500 # Myrs
ext = True
if ext:
    pre = '/media/konstantinos/Dunkey/mesadata/'
else:
    pre = ''
df = pd.read_csv(f'{pre}data/specific_ages/{age}.txt', header = None, 
                              delimiter = '\s', names = ['Planet', 'Profile'])
orb_seps = [0.05, 0.1, 0.2]
ms = grids.n_masses + grids.j_masses
fenvs =  grids.n_envs + grids.j_envs
entropies = [8]
gas_giants_005 = apothicarios()
gas_giants_01 = apothicarios()
gas_giants_02 = apothicarios()

neptunes_01 = apothicarios()
neptunes_02 = apothicarios()
#%% Calc
for m in ms:
    for fenv in fenvs:
        for s in entropies:
            for a in orb_seps:
                modelname = f'm{m}_env{fenv}_zero_a{a}_s{s}'
                if modelname in list(df['Planet']):
                    profile = df.loc[df['Planet'] == modelname]
                    profile_no = profile.iloc[0,1]
                    try:
                        B_dyn, B_dip = hardB_doer_single(f'{pre}data/{modelname}/profile{profile_no}.data')
                        p = mr.MesaData(f'{pre}data/{modelname}/profile{profile_no}.data')
                        r = np.power(10, p.logR)[2] *  c.Rsol / c.Rearth
                        this_peak = peak(B_dip)
                    except ValueError or FileNotFoundError:
                        continue
                    if B_dyn>0 and m>35: 
                        if a == 0.05:
                            gas_giants_005(m, fenv, s, B_dyn, B_dip, r, this_peak)
                        if a == 0.1:
                            gas_giants_01(m, fenv, s, B_dyn, B_dip, r, this_peak)
                        if a == 0.2:
                            gas_giants_02(m, fenv, s, B_dyn, B_dip, r, this_peak)
                    if B_dyn > 0 and m<35:
                        if a == 0.1:
                            neptunes_01(m, fenv, s, B_dyn, B_dip, r, this_peak)
                        if a == 0.2:
                            neptunes_02(m, fenv, s, B_dyn, B_dip, r, this_peak)
                
#%% Plot
# Mass - Bdyn
fig, ax = plt.subplots(1,2, figsize = (6,4), sharex = True, sharey= True) 
ax[0].scatter(gas_giants_005.masses, gas_giants_005.Bdyn, ec = 'k', zorder = 3,
                  marker = 'o', c = gas_giants_005.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

ax[0].scatter(gas_giants_01.masses, gas_giants_01.Bdyn, ec = 'k', zorder = 3,
                  marker = 's', c = gas_giants_01.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

ax[0].scatter(gas_giants_02.masses, gas_giants_02.Bdyn, ec = 'k', zorder = 3,
                  marker = '^', c = gas_giants_02.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

ax[0].scatter(neptunes_01.masses, neptunes_01.Bdyn, ec = 'k', zorder = 3,
              marker = 's',
                  vmin = 0.02, vmax = 0.12, c = neptunes_01.envs, cmap = 'Blues')
ax[0].scatter(neptunes_02.masses, neptunes_02.Bdyn, ec = 'k', zorder = 3,
              marker = '^',
                  vmin = 0.02, vmax = 0.12, c = neptunes_02.envs, cmap = 'Blues')
# <--------------------------------------------------------------------------->
im1 = ax[1].scatter(gas_giants_005.masses, gas_giants_005.Bdip, ec = 'k', zorder = 3,
                  marker = 'o', c = gas_giants_005.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

ax[1].scatter(gas_giants_01.masses, gas_giants_01.Bdip, ec = 'k', zorder = 3,
                  marker = 's', c = gas_giants_01.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

ax[1].scatter(gas_giants_02.masses, gas_giants_02.Bdip, ec = 'k', zorder = 3,
                  marker = '^', c = gas_giants_02.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

im2 = ax[1].scatter(neptunes_01.masses, neptunes_01.Bdip, ec = 'k', zorder = 3,
                    marker = 's',
                  vmin = 0.02, vmax = 0.12, c = neptunes_01.envs, cmap = 'Blues')
ax[1].scatter(neptunes_02.masses, neptunes_02.Bdip, ec = 'k', zorder = 3,
              marker = '^',
                  vmin = 0.02, vmax = 0.12, c = neptunes_02.envs, cmap = 'Blues')

# Limits and scale
ax[0].set_ylim(0,250)
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
fig.suptitle(r'5.5 Gyrs', fontsize = 16, y = 0.95)
fig.text(1.02, 0.35, 'Envelope Fraction', rotation = 90, fontsize = 14)

custom_scatter = [ Line2D([0], [0], color = 'white', linestyle = '',
                          markeredgecolor = 'k', marker = 'o', markersize = 11),
                  Line2D([0], [0], color = 'white', linestyle = '',
                          markeredgecolor = 'k', marker = 's', markersize = 11),
                   Line2D([0], [0], color = 'white',  linestyle = '',
                          markeredgecolor = 'k', marker = '^', markersize = 11),
                 ]
labels = ['0.05 AU', '0.1 AU', '0.2 AU',]
ax[1].legend(custom_scatter, labels, fontsize = 10, ncols = 3,
          bbox_to_anchor=(0.80, 0.88), bbox_transform = fig.transFigure,)
#%% Peak
fig, ax = plt.subplots(1,1, figsize = (4,4), sharex = True, sharey= True) 
ax.scatter(gas_giants_005.masses, gas_giants_005.peak, ec = 'k', zorder = 3,
                  marker = 'o', c = gas_giants_005.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

im1 = ax.scatter(gas_giants_01.masses, gas_giants_01.peak, ec = 'k', zorder = 3,
                  marker = 's', c = gas_giants_01.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

ax.scatter(gas_giants_02.masses, gas_giants_02.peak, ec = 'k', zorder = 3,
                  marker = '^', c = gas_giants_02.envs, cmap = 'Oranges',
                  vmin = 0.85, vmax = 0.96, )

im2 = ax.scatter(neptunes_01.masses, neptunes_01.peak, ec = 'k', zorder = 3,
                    marker = 's',
                  vmin = 0.02, vmax = 0.12, c = neptunes_01.envs, cmap = 'Blues')
ax.scatter(neptunes_02.masses, neptunes_02.peak, ec = 'k', zorder = 3,
              marker = '^',
                  vmin = 0.02, vmax = 0.12, c = neptunes_02.envs, cmap = 'Blues')
ax.grid()
ax.set_xlabel(r'Planet Mass [M$_\oplus$]', fontsize = 14)
ax.set_ylabel(r'Peak Frequency [MHz]', fontsize = 14)
fig.colorbar(im1, cax = fig.add_axes([0.92, 0.12, 0.02, 0.36]),)
fig.colorbar(im2, cax = fig.add_axes([0.92, 0.52, 0.02, 0.36]),)
fig.text(1.05, 0.35, 'Envelope Fraction', rotation = 90, fontsize = 14)
ax.set_ylim(0,450)
ax.legend(custom_scatter, labels, fontsize = 10, ncols = 1,
          bbox_to_anchor=(0.98, 0.28), bbox_transform = ax.transAxes,)
#%%

fig, ax = plt.subplots(1,2, figsize = (6,4), sharex = True, sharey= True) 
ax[0].scatter(gas_giants_005.radii, gas_giants_005.Bdyn, ec = 'k', zorder = 3,
                  marker = 'o', c = gas_giants_005.envs, cmap = 'Blues',
                  vmin = 0.85, vmax = 1, )

ax[0].scatter(gas_giants_01.radii, gas_giants_01.Bdyn, ec = 'k', zorder = 3,
                  marker = '^', c = gas_giants_01.envs, cmap = 'Blues',
                  vmin = 0.85, vmax = 1, )

ax[0].scatter(gas_giants_02.radii, gas_giants_02.Bdyn, ec = 'k', zorder = 3,
                  marker = 's', c = gas_giants_02.envs, cmap = 'Blues',
                  vmin = 0.85, vmax = 1, )

# ax[0].scatter(neptunes.radii, neptunes.Bdyn, ec = 'k', zorder = 3,
#                   vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')
im1 = ax[1].scatter(gas_giants_005.radii, gas_giants_005.Bdip, ec = 'k', zorder = 3,
                  marker = 'o', c = gas_giants_005.envs, cmap = 'Blues',
                  vmin = 0.85, vmax = 1, )

ax[1].scatter(gas_giants_01.radii, gas_giants_01.Bdip, ec = 'k', zorder = 3,
                  marker = '^', c = gas_giants_01.envs, cmap = 'Blues',
                  vmin = 0.85, vmax = 1, )

ax[1].scatter(gas_giants_02.radii, gas_giants_02.Bdip, ec = 'k', zorder = 3,
                  marker = 's', c = gas_giants_02.envs, cmap = 'Blues',
                  vmin = 0.85, vmax = 1, )

# im2 = ax[1].scatter(neptunes.masses, neptunes.Bdip, ec = 'k', zorder = 3,
#                   vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')

# Limits and scale
# ax[0].set_ylim(40,220)
# ax[0].set_xscale('log')

# Labels & Grid
ax[0].grid()
ax[0].set_xlabel(r'Planet Radius [R$_\oplus$]', fontsize = 14)
ax[0].set_ylabel(r'Dynamo [G]', fontsize = 14)
ax[1].grid()
ax[1].set_xlabel(r'Planet Radius [R$_\oplus$]', fontsize = 14)
ax[1].set_ylabel(r'Dipole [G]', fontsize = 14)
# Fig
fig.colorbar(im1)#, cax = fig.add_axes([0.92, 0.12, 0.02, 0.36]),)
# fig.colorbar(im2, cax = fig.add_axes([0.92, 0.52, 0.02, 0.36]),)
fig.suptitle(r'5.5 Gyrs, $\alpha$ 0.1 AU', fontsize = 16, y = 0.95)
fig.text(0.98, 0.35, 'Envelope Fraction', rotation = 90, fontsize = 14)
#%% Plot
# Rad - Bdyn
fig, ax = plt.subplots(1,2, figsize = (6,4), sharex = True, sharey= True) 
ax[0].scatter(gas_giants.radii, gas_giants.Bdyn, ec = 'k', zorder = 3,
                  vmin = 60, vmax = 100, c = gas_giants.envs, cmap = 'Blues')
ax[0].scatter(neptunes.radii, neptunes.Bdyn, ec = 'k', zorder = 3,
                  vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')

im1 = ax[1].scatter(gas_giants.radii, gas_giants.Bdip, ec = 'k', zorder = 3,
                  vmin = 60, vmax = 100, c = gas_giants.envs, cmap = 'Blues')
im2 = ax[1].scatter(neptunes.radii, neptunes.Bdip, ec = 'k', zorder = 3,
                  vmin = 1, vmax = 15, c = neptunes.envs, cmap = 'Blues')

# Limits and scale
ax[0].set_xlim(2, 13)
ax[0].set_ylim(0, 160)
# ax[0].set_xscale('log')

# Labels & Grid
ax[0].grid()
ax[0].set_xlabel(r'Planet Radius [M$_\oplus$]', fontsize = 14)
ax[0].set_ylabel(r'Dynamo [G]', fontsize = 14)
ax[1].grid()
ax[1].set_xlabel(r'Planet Radius [M$_\oplus$]', fontsize = 14)
ax[1].set_ylabel(r'Dipole [G]', fontsize = 14)
# Fig
fig.colorbar(im1, cax = fig.add_axes([0.92, 0.12, 0.02, 0.36]),)
fig.colorbar(im2, cax = fig.add_axes([0.92, 0.52, 0.02, 0.36]),)
fig.suptitle(r'5.5 Gyrs, $\alpha$ 0.1 AU', fontsize = 16, y = 0.95)
fig.text(0.98, 0.35, 'Envelope Fraction', rotation = 90, fontsize = 14)
#%% Radii - Fenvs
fig, ax = plt.subplots(1,2, figsize = (5,4)) 
im1 = ax[0].scatter(gas_giants.radii, gas_giants.envs, ec = 'k', zorder = 3,
                 vmin = 50, vmax = 200, c = gas_giants.Bdyn, cmap = 'Blues')
cb = fig.colorbar(im1)
im2 = ax[1].scatter(neptunes.radii, neptunes.envs, ec = 'k', zorder = 3,
                  vmin = 1, vmax = 50, c = neptunes.Bdyn, cmap = 'Blues')
fig.colorbar(im2)
# im1 = ax[1].scatter(gas_giants.radii, gas_giants.envs, ec = 'k', zorder = 3,
#                  vmin = 60, vmax = 100, c = gas_giants.Bdip, cmap = 'Blues')
# im2 = ax[1].scatter(neptunes.radii, neptunes.envs, ec = 'k', zorder = 3,
#                  vmin = 1, vmax = 15, c = neptunes.Bdip, cmap = 'Blues')

# Limits and scale
# ax[0].set_ylim(0,160)
ax[0].set_xlim(7,13)
ax[1].set_xlim(1, 8)
ax[0].set_ylim(55, 102)
ax[1].set_ylim(0, 17)
#ax[0].set_xscale('log')

# Labels & Grid
ax[0].grid()
ax[1].grid()
ax[0].set_xlabel(r'Planet Radius [R$_\oplus$]', fontsize = 14)
ax[1].set_xlabel(r'Planet Radius [R$_\oplus$]', fontsize = 14)
ax[0].set_ylabel(r'Envelope Fraction', fontsize = 14)
# ax[1].grid()
# ax[1].set_xlabel(r'Planet Radius [R$_\oplus$]', fontsize = 14)
# ax[1].set_ylabel(r'Dipole [G]', fontsize = 14)
# Fig

# ax[1].colorbar(im2, cax = fig.add_axes([0.92, 0.52, 0.02, 0.36]),)
fig.suptitle(r'5.5 Gyrs, $\alpha$ 0.1 AU, Zero Evaporation', fontsize = 16, y = 0.95)
fig.text(0.92, 0.35, 'Dynamo [G]', rotation = 90, fontsize = 14,)