#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 16:05:45 2024

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 20:49:45 2024

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mesa_reader as mr
import os
import src.prelude as c
from src.reynolds import profile_sorter, rey_mag

class apothicarios:
    def __init__(self, name):
        self.name = name
        self.age = []
        self.lum_mesa = []
        self.lum_stefan = [] # comparing MESA logT with photosphereL
        self.lum_dyn = []
        self.lum_rcb = []
        self.rcb = []
        self.dyn = []
        self.radius = []
        self.T_rcb = []
        self.T_dyn = []
        self.T_radius = []
      
    def __call__(self, age_h, lum_mesa_h, lum_stefan_h, lum_dyn_h, lum_rcb_h):
        self.age.append(age_h)
        self.lum_mesa.append(lum_mesa_h)
        self.lum_stefan.append(lum_stefan_h)
        self.lum_dyn.append(lum_dyn_h)
        self.lum_rcb.append(lum_rcb_h)

    def rs(self, radius_h, dyn_h, rcb_h):
        self.radius.append(radius_h)
        self.dyn.append(dyn_h)
        self.rcb.append(rcb_h)
            
    def Ts(self, T_surface_h, T_dyn_h, T_rcb_h):
        self.T_radius.append(T_surface_h)
        self.T_dyn.append(T_dyn_h)
        self.T_rcb.append(T_rcb_h)
      
def doer(names):
    # Count and generate profile lists
    apothikh = []
    for name in names:
        # History data wrangling
        # h_path = 'data/' + name + '/history_7.data'
        # h = mr.MesaData(h_path)
        # print(dir(h))
        # h_age = np.round(10**h.log_star_age /  1e9, 2) # Gyr
        
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
            r = np.power(10, p.logR)
            age = np.round(p.star_age /  1e6, 2)

            # The Phantom Luminosity. MESA
            L_mesa = p.photosphere_L
            #L_mesa /= c.Lsol
            
            # The Luminosity Wars. BB on surface
            r_cgs = r[0] * c.Rsol
            L_stefan = 4 * np.pi * r_cgs**2 * c.sigma* np.power(10, p.logT[0])**4
            L_stefan /= c.Lsol
            
            # Revenge of the Luminosity. BB on dynamo surface
            _, reynolds_mag_number, _= rey_mag(p)
            R_dynamo_active = r[reynolds_mag_number > c.critical_rey_mag_num]
            Rdyn = R_dynamo_active[0]
            idx = np.argmin(np.abs(r - Rdyn)) 
            Rdyn *= c.Rsol # [cgs]
            L_dyn = 4 * np.pi * Rdyn**2 * c.sigma* np.power(10, p.logT[idx])**4
            L_dyn /= c.Lsol
            # L_dyn = p.luminosity[idx]

            # A New Luminosity. BB on Radiative Convective Boundary
            rad = p.gradr
            conv = p.grada
            cmorethanb = r[conv < rad]
            RCB = cmorethanb[0]
            print(RCB)
            idx_rcb = np.argmin(np.abs(r - RCB))
            RCB *= c.Rsol
            L_rcb = 4 * np.pi * RCB**2 * c.sigma* np.power(10, p.logT[idx_rcb])**4
            L_rcb /= c.Lsol
            
            # Save
            hold(age, L_mesa, L_stefan, L_dyn, L_rcb)
            hold.rs(r_cgs, Rdyn, RCB)
            hold.Ts(np.power(10, p.logT[0]), np.power(10, p.logT[idx]),
                    np.power(10, p.logT[idx_rcb]))
        apothikh.append(hold)
    return apothikh


def plotter(names):
   
    colors = [c.c91, 'darkorange', 'k', 'r',]
    
    # Makes the calculations
    planets = doer(names) 
    
    fig, axs = plt.subplots(1,3, tight_layout = True, sharex = True,
                           figsize = (8,4))
    
    for planet in planets:
        axs[0].plot(planet.age, planet.radius, color = colors[1], 
                 linestyle = '-', label = 'Surface')
        axs[0].plot(planet.age, planet.dyn, color = colors[2],
                 linestyle = '-.', label = r'$R_{dyn}$')
        axs[0].plot(planet.age, planet.rcb, color = colors[3],
                 linestyle = ':', label = r'$R_{rcb}$')
    
        axs[1].plot(planet.age, planet.T_radius, color = colors[1], 
                 linestyle = '-', label = 'Surface')
        axs[1].plot(planet.age, planet.T_dyn, color = colors[2],
                 linestyle = '-.', label = r'$R_{dyn}$')
        axs[1].plot(planet.age, planet.T_rcb, color = colors[3],
                 linestyle = ':', label = r'$R_{rcb}$')
        
        axs[2].plot(planet.age, planet.lum_mesa, color = colors[0], 
                  label = 'MESA')
        axs[2].plot(planet.age, planet.lum_stefan, color = colors[1], 
                 linestyle = '-', label = 'BB Surface')
        axs[2].plot(planet.age, planet.lum_dyn, color = colors[2],
                 linestyle = '-.', label = r'BB $R_{dyn}$')
        axs[2].plot(planet.age, planet.lum_rcb, color = colors[3],
                 linestyle = ':', label = r'BB $R_{rcb}$')
        
    # Make nice
    axs[2].set_ylabel('Luminosity [$L_\odot$]', fontsize = 12)
    axs[1].set_ylabel('Temperature [K]', fontsize = 12)
    axs[0].set_ylabel('Radii [$R_\oplus$]', fontsize = 12)

    fig.text(0.5, -0.02, r'Age [Myr]', 
              fontsize = 15, transform = fig.transFigure)
    
    #axs.set_ylim(-1e-5,1e-5)
    axs[0].set_yscale('log')
    fig.suptitle('Different kinds of luminosities on hot jupiter')
    
    # Legend
    custom_lines = []
    custom_lines.append( Line2D([0], [0], linestyle = '-',
                         color = colors[1], label = 'Surface'))
    custom_lines.append( Line2D([0], [0], linestyle = '-.',
                         color = colors[2], label = r'$R_{dyn}$'))
    custom_lines.append( Line2D([0], [0], linestyle = ':',
                         color = colors[3], label = r'$R_{rcb}$'))
    labels = ['Surface', r'$R_{dyn}$', r'$R_{rcb}$' ]
    fig.legend(custom_lines, labels, 
               fontsize =  14, ncols = 3, alignment = 'left', 
               bbox_to_anchor=(0.77, -0.03), bbox_transform = fig.transFigure)

#%%    

name = 'jup17'
labels = ['Hot Jup']
plotter([name])
