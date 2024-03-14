#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 16:06:02 2024

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:45:21 2024

@author: konstantinos
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
# 
import src.prelude as c
from src.Bfield.reynolds import profile_sorter, rey_mag
from src.Bfield.hardB import convective_flux, temperature_scale_height

#%%
def q_calc(p, Rdyns):
    r = np.power(10, p.logR) * c.Rsol / c.Rearth
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

    return q

class apothicarios:
    def __init__(self, name, age):
        self.name = name
        self.age = age
        self.Rdyn = []
        self.q = []
        
    def __call__(self, r, qh):
        self.Rdyn.append(r)
        self.q.append(qh)
        
def doer(name, many = 9):
    apothikh = []
    # Get profiles
    p_path = 'data/' + name
    profiles = os.popen('ls ' + p_path + '/profile*.data').read()
    profiles = list(profiles.split("\n"))
    profiles.pop() # Remove last
    profiles = profile_sorter(profiles) # Guess what that does

    # Choose how many profiles you want
    step = len(profiles) // many
    indices = np.arange(0, len(profiles), step, dtype=int)
    # Do the thing
    for i in indices:
        p = mr.MesaData(profiles[i]) 
        r, reynolds_mag_number, age = rey_mag(p)
        r *= c.Rsol / c.Rearth
        # Get Rdyn
        R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
        try:
            Rdyn_end = R_dynamo_active[-1] # it's the wrong way round
        except IndexError:
            # Save
            hold = apothicarios(name, age)
            apothikh.append(hold)
            continue
        hold = apothicarios(name, age) # Instanciate Holder Object
        q = q_calc(p, R_dynamo_active)
        hold(R_dynamo_active, q) # Save

        # Keep em
        apothikh.append(hold)
    return apothikh
#%%
fig, axs = plt.subplots(2, 4, figsize = (12, 5), tight_layout = True, 
                        sharex = True, sharey = True)
planet = doer('nep_e1_zero')
fig.suptitle('Neptune 10 $M_\oplus$, 10 $\%$ $f_{env}$, 0.1 AU, No Escape', 
             fontsize = 20)
for ax, prof in zip(axs.reshape(-1), range(1,9)):
    try:
        ax.plot(planet[prof].Rdyn[0][1:], planet[prof].q[0], c ='k')
    except IndexError:
        continue
    ax.set_title(str(planet[prof].age) + ' Myrs')
    ax.grid()
    
#fig.text($R_\oplus$]')
fig.text(-0.008, 0.42, r'$F^{2/3}$', fontsize = 20, rotation = 90,
         transform = fig.transFigure)
fig.text(1.001, 0.35, r'Color Terms', fontsize = 20, rotation = 270,
         transform = fig.transFigure)
fig.text(0.45, -0.03, r'$R [R_\oplus]$', fontsize = 20, rotation = 0,
         transform = fig.transFigure)

labels = ['Convective Flux']
lines = [Line2D([0], [0], color = 'k',)]
fig.legend(lines, labels,
        fontsize = 15, ncols = 1,
        bbox_to_anchor=(0.8, 0.03), bbox_transform = fig.transFigure,)