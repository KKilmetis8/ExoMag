#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:32:06 2024

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.Bfield.reynolds import rey_mag, dynamo_region
from src.Utilities.profile_sorter import profile_sorter
from src.Bfield.hardB import convective_flux
from src.Bfield.hori import hori_doer

def q_explore(p, Rdyns):
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
    # Reference flux
    q0 = q[-1] # * Rdyns[-1]**2 / Rdyns[0]**2

    return q, q0, Rdyns[-1] / r[0], Rdyns[0] / r[0]
 
class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.rmax = []
        self.q = []
        self.q0 = []
        self.dyn_start = []
        self.dyn_end = []
        self.reynolds = []
    def __call__(self, age_h, q_h, q0_h, dyn_start_h, dyn_end_h, reynolds_h = 0):#), lum2_h):
        self.age.append(age_h)
        self.q.append(q_h)
        self.q0.append(q0_h)
        self.dyn_start.append(dyn_start_h)
        self.dyn_end.append(dyn_end_h)
        self.reynolds.append(reynolds_h)

def doer(names):
    # Count and generate profile lists
    apothikh = []
    for name in names:
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
            # r, rmn, age = rey_mag(p)
            # r *= c.Rsol / c.Rearth
            # Get Rdyn
            R_dynamo_active, rmn, age = dynamo_region(p)
            # R_dynamo_active = r[rmn > c.critical_rey_mag_num]
            # rmn = rmn[rmn > c.critical_rey_mag_num]
            try:
                Rdyn_end = R_dynamo_active[0] # it's the wrong way round
                R_dynamo_active *= c.Rsol / c.Rearth
            except IndexError:
                # Save
                hold(age, 0, 0, 0, 0 ,0)
                continue
            q, q0, Rdyns_start, Rdyn_end = q_explore(p, R_dynamo_active)
            hold(age, q, q0, Rdyns_start, Rdyn_end, rmn[0])#, L2)
        apothikh.append(hold)
    return apothikh


def plotter(names, cols, title):
    # Makes the calculations
    planets = doer(names) 
    hori = hori_doer(names, 100)[0]

    fig, ax = plt.subplots(1,1, tight_layout = True, sharex = True,
                           figsize = (4,4))
    
    custom_lines = []
    for planet in planets:
        #axs[0].plot(planet.age, planet.q0, color = 'k')
        ax.plot(planet.age, planet.dyn_start, 'k')
        ax.plot(planet.age, planet.dyn_end, 'k', linestyle = ':')
        
        #ax.plot(hori.age)
        
        ax2=ax.twinx()
        # ax2.plot(planet.age, planet.reynolds, color = c.c97)
        ax2.set_ylabel(r'Reynolds Mag', fontsize = 14 , rotation = 270, labelpad=16)
        
    # import Mors as mors 
    # axlum = ax.twinx()
    # star = mors.Star(Mstar=1.0, Omega=8.0)
    # axlum.plot(star.Tracks['Age'], np.log10(star.Tracks['Lbol']))
    # axlum.set_ylabel(r'$\log(L_{bol} [erg/s])$ ', fontsize = 14, rotation = 270, labelpad = 17)
    # # rotation = 270, labelpad = 17)

    # Make nice
    # ax2.set_ylabel(r'Convective Flux at Rdyn Start [W/m$^2$]', fontsize = 14, 
    # rotation = 270, labelpad = 17)
    ax.set_ylabel('Distance from core [r/R$_p$]', fontsize = 14,
                    labelpad = 17)
    # ax2.set_yscale('log')
    # ax.grid()
    # ax2.grid()
    
    # axs[0].set_yscale('log')
    # axs[2].set_yscale('log')
    ax.set_xlim(200, 11_000)
    # axs[1].set_ylim(7.2,14)
    # axs[2].set_ylim(0,1.5e-5)

    fig.suptitle(title, 
                  fontsize = 16, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)
    
    # Legend
    custom_lines.append( Line2D([0], [0], c = c.reddish))
    custom_lines.append( Line2D([0], [0], c = c.cyan))
    custom_lines.append( Line2D([0], [0], c = c.c97))
    labels = ['End', 'Start', 'Reynolds Mag']
    fig.legend(custom_lines, labels,
            fontsize =  9, ncols = 3, alignment = 'center', # Lawful Neutral
            bbox_to_anchor=(0.87, -0.04), bbox_transform = fig.transFigure,)  
#%%    
kind = 'saturn90'
both = True # hori or reynolds
if __name__ == '__main__':
    if kind == 'jup':
        name = 'jup_e94_zero'
        label = 'Jupiter'
        plotter([name], 1, r'Jupiter $M_\oplus$ 317, f$_0$ = 94 $\%$, a = 0.1 AU, Zero')
    if kind == 'jup96':
        name = 'jup_e96_zero'
        label = 'Jupiter'
        plotter([name], 1, r'Jupiter $M_\oplus$ 317, f$_0$ = 96 $\%$, a = 0.1 AU, Zero')
    if kind == 'jup92':
        name = 'jup_e92_zero'
        label = 'Jupiter'
        plotter([name], 1, r'Jupiter $M_\oplus$ 317, f$_0$ = 92 $\%$, a = 0.1 AU, Zero')
    if kind == 'jup_autoS':
        name = 'jup_e94_zero_autoS'
        label = 'Jupiter'
        plotter([name], 1, r'Jupiter $M_\oplus$ 317, f$_0$ = 94 $\%$, a = 0.1 AU, Zero AutoS',)
    if kind == 'nep10':
        name = 'nep_e10_zero_7s'
        label = 'Jupiter'
        plotter([name], 1, r'Neptune $M_\oplus$ 17, f$_0$ = 10 $\%$, a = 0.1 AU, Zero',)
    if kind == 'nep1':
        name = 'nep_e1_zero'
        label = 'Jupiter'
        plotter([name], 1, r'Neptune $M_\oplus$ 17, f$_0$ = 1 $\%$, a = 0.1 AU, Zero',)
    if kind == 'saturn95':
        name = 'm95_e95_zero_a01_s8'
        label = 'saturn'
        plotter([name], 1, r'Saturn $M_\oplus$ 95, f$_0$ = 95 $\%$, a = 0.1 AU, Zero',)
    if kind == 'saturn90':
        name = 'm95_e90_zero_a01_s8'
        label = 'saturn'
        plotter([name], 1, r'Saturn $M_\oplus$ 95, f$_0$ = 90 $\%$, a = 0.1 AU, Zero',)
    if kind == 'saturn85':
        name = 'm95_e85_zero_a01_s8'
        label = 'saturn'
        plotter([name], 1, r'Saturn $M_\oplus$ 85, f$_0$ = 90 $\%$, a = 0.1 AU, Zero',)