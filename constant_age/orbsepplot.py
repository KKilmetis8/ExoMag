#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:32:57 2024

@author: konstantinos
"""

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

# Choc
import src.prelude as c
import src.Utilities.planet_grids as grids
from src.Bfield.hardB import hardB_doer_single

class apothicarios:
    def __init__(self):
        self.masses = []
        self.orbsep = []
        self.entropies = []
        self.Bdyn = []
        self.Bdip = []
        self.radii = []
    def __call__(self, m, orbsep, s, dyn, dip, r):
        self.masses.append(m)
        self.orbsep.append(orbsep)
        self.entropies.append(s)
        self.Bdyn.append(dyn)
        self.Bdip.append(dip)
        self.radii.append(r)

        
age = 1500 # Myrs
df = pd.read_csv(f'data/specific_ages/{age}.txt', header = None, 
                              delimiter = '\s', names = ['Planet', 'Profile'])
orb_seps1 = list(np.arange(0.05, 0.3, 0.01))
orb_seps2 = list(np.arange(0.3, 2, 0.1))
orb_seps = orb_seps1 + orb_seps2
ms = [317, 17]
fenvs =  [0.94, 0.06] # grids.j_envs
entropies = [8]
escapes = ['zero', 'EL', 'HD']
gas_giants = apothicarios()
nep_zero = apothicarios()
nep_EL = apothicarios()
nep_HD = apothicarios()
#%% Calc
for m in ms:
    for fenv in fenvs:
        for s in entropies:
            for escape in escapes:
                for a in orb_seps:
                    modelname = f'm{m}_env{fenv}_{escape}_a{a}_s{s}'
                    if modelname in list(df['Planet']):
                        profile = df.loc[df['Planet'] == modelname]
                        profile_no = profile.iloc[0,1]
                        if int(profile_no) > 0:
                            try:
                                B_dyn, B_dip, age_this = hardB_doer_single(f'data/{modelname}/profile{profile_no}.data')
                                print(age_this)
                                if np.abs(age-age_this)>500:
                                    continue
                                p = mr.MesaData(f'data/{modelname}/profile{profile_no}.data')
                                r = np.power(10, p.logR)[2] *  c.Rsol / c.Rearth
                            except ValueError or FileNotFoundError:
                                continue
                            if B_dyn>0:
                                if m>300:
                                    gas_giants(m, a, s, B_dyn, B_dip, r)
                                else:
                                    if escape == 'zero':
                                        nep_zero(m, a, s, B_dyn, B_dip, r)
                                    elif escape == 'EL':
                                        nep_EL(m, a, s, B_dyn, B_dip, r)
                                    elif escape == 'HD':
                                        nep_HD(m, a, s, B_dyn, B_dip, r)
#%% Plot
# ---------------- Dynamo ----------------
fig, ax = plt.subplots(1,2, figsize = (6,4), sharex = True, sharey= True) 
ax[0].plot(gas_giants.orbsep, gas_giants.Bdyn, 
              c = 'darkorange', ls = '-', lw = 1, 
              markeredgecolor = 'k', marker = 'h', markersize = 7)
ax[0].plot(nep_zero.orbsep, nep_zero.Bdyn, 
              c = 'navy', ls = '-', lw = 1, 
              marker = 'h', markersize = 7)
ax[0].plot(nep_EL.orbsep, nep_EL.Bdyn, 
              c = 'b', ls = '-', lw = 1, 
              marker = 'h', markersize = 7)
ax[0].plot(nep_HD.orbsep, nep_HD.Bdyn, 
              c = 'skyblue', ls = '-', lw = 1, 
              markeredgecolor = 'k', marker = 'h', markersize = 7)
# ---------------- Dipole ----------------
ax[1].plot(gas_giants.orbsep, gas_giants.Bdip,  
                   c = 'darkorange', ls = '-', lw = 1, 
                   markeredgecolor = 'k', marker = 'h', markersize = 7)

ax[1].plot(nep_zero.orbsep, nep_zero.Bdip, 
              c = 'navy', ls = '-', lw = 1, 
              marker = 'h', markersize = 7)

ax[1].plot(nep_EL.orbsep, nep_EL.Bdip, 
              c = 'b', ls = '-', lw = 1, 
              marker = 'h', markersize = 7)

ax[1].plot(nep_HD.orbsep, nep_HD.Bdip, 
              c = 'skyblue', ls = '-', lw = 1, 
              markeredgecolor = 'k', marker = 'h', markersize = 7)

# Labels & Grid
ax[0].grid()
ax[0].set_xlabel(r'Orbital Separation  [AU]', fontsize = 14)
ax[0].set_ylabel(r'$B_\mathrm{dyn}$ [G]', fontsize = 14)
ax[1].grid()
ax[1].set_xlabel(r'Orbital Separation  [AU]', fontsize = 14)
ax[1].set_ylabel(r'$B_\mathrm{dip}$ [G]', fontsize = 14)
ax[0].set_xlim(0.045,0.1)
ax[0].set_ylim(1,21)
# ax[0].set_yscale('log')
# Fig
fig.suptitle(r'1.5 Gyrs', fontsize = 16, y = 0.95)


