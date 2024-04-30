#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 13:01:20 2024

@author: konstantinos
"""

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
ms = [317]#, 17]
fenvs =  [0.94]#, 0.06] # grids.j_envs
entropies = [8]
escapes = ['zero', 'EL', 'HD']
gas_giants_1 = apothicarios()
gas_giants_5 = apothicarios()
gas_giants_10 = apothicarios()

#%% Calc
for age in ages:
    df = pd.read_csv(f'data/specific_ages/{age}.txt', header = None, 
                     delimiter = '\s', names = ['Planet', 'Profile'])
    for planet, profile_no in zip(df['Planet'], df['Profile']):
        modelname = planet
        underscores = [ i for i, x in enumerate(modelname) if x == '_']
        m = int( planet[ 1: underscores[0] ])
        escape = planet[ underscores[1]+1:underscores[2] ]
        a = float( planet[ underscores[2]+2:underscores[3]])

        if int(profile_no) > 0 and m == 317:
            try:
                B_dyn, B_dip, age_this = hardB_doer_single(f'data/{modelname}/profile{profile_no}.data')
                #print(age_this)
                if np.abs(age-age_this) > 500:

                    continue
                p = mr.MesaData(f'data/{modelname}/profile{profile_no}.data')
                r = np.power(10, p.logR)[2] *  c.Rsol / c.Rearth
            except ValueError or FileNotFoundError:
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
 

# Labels & Grid
ax[0].grid()
ax[0].set_xlabel(r'Orbital Separation  [AU]', fontsize = 14)
ax[0].set_ylabel(r'$B_\mathrm{dyn}$ [G]', fontsize = 14)
ax[1].grid()
ax[1].set_xlabel(r'Orbital Separation  [AU]', fontsize = 14)
ax[1].set_ylabel(r'$B_\mathrm{dip}$ [G]', fontsize = 14)
ax[0].set_xlim(0,2)
#ax[0].set_ylim(75,250)
# Fig
# fig.suptitle(r'5.5 Gyrs', fontsize = 16, y = 0.95)


