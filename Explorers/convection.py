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
from src.Bfield.reynolds import profile_sorter, rey_mag
from src.Bfield.hardB import convective_flux

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
    def __call__(self, age_h, rmax_h, q_h, q0_h, dyn_start_h, dyn_end_h):#), lum2_h):
        self.age.append(age_h)
        self.rmax.append(rmax_h)
        self.q.append(q_h)
        self.q0.append(q0_h)
        self.dyn_start.append(dyn_start_h)
        self.dyn_end.append(dyn_end_h)

      
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
            r, reynolds_mag_number, age = rey_mag(p)
            r *= c.Rsol / c.Rearth
            # Get Rdyn
            R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
            try:
                Rdyn_end = R_dynamo_active[0] # it's the wrong way round
            except IndexError:
                # Save
                hold(age, 0, 0, 0)
                continue
            q, q0, Rdyns_start, Rdyn_end = q_explore(p, R_dynamo_active)
            hold(age, r[0], q, q0, Rdyns_start, Rdyn_end)#, L2)
        apothikh.append(hold)
    return apothikh


def plotter(names, cols, title):
    # Makes the calculations
    planets = doer(names) 

    fig, ax = plt.subplots(1,1, tight_layout = True, sharex = True,
                           figsize = (4,4))
    
    custom_lines = []
    for planet in planets:
        ax.plot(planet.age, planet.q0, color = c.reddish)
        #axs[0].plot(planet.age, planet.q0, color = 'k')
        ax2=ax.twinx()
        ax2.plot(planet.age, planet.dyn_start, c.cyan)
        ax2.plot(planet.age, planet.dyn_end, c.kroki,)
        
    #Legend
    custom_lines.append( Line2D([0], [0], c = c.reddish, label = 'q'))
    custom_lines.append( Line2D([0], [0], c = c.cyan, label = 'Start'))
    custom_lines.append( Line2D([0], [0], c = c.kroki, label = 'End'))

    
    # Make nice
    ax.set_xlabel('Age [Myrs]')
    ax.set_ylabel(r'Convective Flux at Rdyn Start [W/m$^2$]', fontsize = 14)
    ax2.set_ylabel('Distance from core [r/R$_p$]', fontsize = 14, rotation = 270,
                   labelpad = 17)
    ax2.set_yscale('log')
    # ax.grid()
    # ax2.grid()
    
    # axs[0].set_yscale('log')
    # axs[2].set_yscale('log')
    # axs[0].set_xlim(100, 10_000)
    # axs[1].set_ylim(7.2,14)
    # axs[2].set_ylim(0,1.5e-5)

    fig.suptitle(title, 
                  fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r' Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)
    fig.legend(custom_lines, ['Start', 'End', 'q(Rdyn_start)'],
            fontsize =  9, ncols = 1, alignment = 'center', # Lawful Neutral
            bbox_to_anchor=(0.83, 0.65), bbox_transform = fig.transFigure,)  
#%%    
kind = 'jup'
if __name__ == '__main__':
    if kind == 'jup':
        name = 'jup_e95'
        label = 'Jupiter'
        plotter([name], 1, r'Jupiter $M_\odot$ 317, f$_0$ = 90 $\%$, a = 0.1, EL')
    if kind == 'jupenv':
        name1 = 'jup_e30'
        name2 = 'jup_e40'
        name3 = 'jup_e50'
        name4 = 'jup_e60'
        name5 = 'jup_e70'
        name6 = 'jup_e80'
        name7 = 'jup_e90'
        name8 = 'jup_e95'
        name9 = 'jup_e98'
        names = [name1, name2, name3, name4, name5, name6, name7, name8, name9]
        labels = ['30', '40', '50', '60', '70', '80', '90', '95', '98']
        plotter(names, 3, labels, 'Jupiter with Diff. Envelopes', 
                met_hyd_lines=True)
    if kind == 'nepenv':
        name1 = 'nep_1'
        name2 = 'nep_2'
        name3 = 'nep_3'
        name4 = 'nep_4'
        name5 = 'nep_5'
        name6 = 'nep_6'
        names = [name1, name2, name3, name4, name5, name6]
        labels = ['10', '20', '30', '40', '50', '60']
        plotter(names, 3, labels, 'Neptunes with Diff. Envelopes')
    if kind == 'se5env':
        name1 = 'se58'
        name2 = 'se59'
        name3 = 'se510'
        name4 = 'se511'
        name5 = 'se512'
        name6 = 'se513'
        name7 = 'se514'
        names = [name1, name2, name3, name4, name5, name6, name7]
        labels = ['1', '2', '6', '12', '18', '24', '30']
        plotter(names, 3, labels, '5xEarth with Diff. Envelopes')
    if kind == 'se5envEL':
        name1 = 'se5_e001_a03_EL'
        name2 = 'se5_e002_a03_EL'
        name3 = 'se5_e005_a03_EL'
        name4 = 'se5_e01_a03_EL'
        name5 = 'se5_e013_a03_EL'
        name6 = 'se5_e017_a03_EL'
        names = [name1, name2, name3, name4, name5, name6]
        labels = ['1', '2', '5', '10', '13', '17',]
        plotter(names, 3, labels, '5xEarth with Diff. Envelopes')
    if kind == 'jup-sep':
        name1 = 'jup_e95_a002'
        name3 = 'jup_e95_a008_2'
        name4 = 'jup_e95_a012'
        name5 = 'jup_e95_a016'
        name6 = 'jup_e95_a020'
        name7 = 'jup_e95_a024'
        name8 = 'jup_e95_a028'
        names = [name1, name3, name4, name5, name6, name7, name8]
        labels = ['0.02', '0.08', '0.12', '0.16', '0.20', '0.24', '0.28',]
        plotter(names, 3, labels, 'Jups far away')