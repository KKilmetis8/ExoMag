#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:03:36 2023

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 17:19:57 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [6 , 5]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
AEK = '#F1C410'
import mesa_reader as mr
# Constants
Rsol = 6.957e10 # [cm]
Rjup = 71.492e8 # [cm]
Rearth = 6.371e8 # [cm]
Msol = 1.98855e33 # [g]
Mjup = 1.8986e30 # [g]
Mjup_in_sol = Mjup / Msol
Mearth = 5.9722e27 # [g]
c = 2.99792458e10 #[cm/s]
h = 6.62607015e-27 #[gcm^2/s]
Kb = 1.380649e-16 #[gcm^2/s^2K]
Lsol = 3.826e33 # [erg/s]
sigma = 5.670e-5 # [cgs]
sqrt2 = np.sqrt(2)
#%%
path1 = 'xodros-puffy2'
title1 = r'320 M$_\oplus$, f$_{at}$ = 0.5, $\alpha$ = 0.05 AU, $M_*$ = 0.833 M$_\odot$, P$_*$ = 1 day'
path2 = 'xodros-noirr'
title2 = r'320 M$_\oplus$, f$_{at}$ = 0.5, $\alpha$ = 10 AU, $M_*$ = 0.1 M$_\odot$, P$_*$ = 60 days'
back = 'skyblue'
def cont_breaker(above, rs):
    ''' Given the mask where a>b, returns where this is in array c 
    as a list of lists'''
    true_ranges = []
    start = None

    for i, value in enumerate(above):
        if value == 1:
            if start is None:
                start = i
        elif start is not None:
            true_ranges.append(list(range(start, i)))
            start = None

    if start is not None:
        true_ranges.append(list(range(start, len(above))))

    result = []
    for ran in true_ranges:
        result.append(rs[ran])
    return result

def energy_transpport(path, title):
    back = '#cae8fa'
    p_path = 'data/' + path  + '/profile'
    fig, axs = plt.subplots(2,2, tight_layout = True, sharex=True,
                           figsize = (7,7))
    profiles = np.arange(3, 6+1, step = 1)
    profile_interval = 1
    for p, ax in zip(profiles, axs.reshape(-1)):
        
        # Load data
        p = mr.MesaData(p_path + str(p) + '.data')
        
        r = np.power(10, p.logR) * Rsol / Rearth
        
        grada = np.log10(p.grada)
        gradr = np.log10(p.gradr)
        
        # Plot
        ax.plot(r, grada, c = 'k', label = r'$\nabla_a$')
        ax.plot(r, gradr, c = 'b', label = r'$\nabla_r$')
        
        above = np.where(gradr > grada, 1, 0)
        mask = cont_breaker(above, r)
        for m in mask:
            ax.axvspan(m[0], m[-1], facecolor=back, alpha = 0.75)
        # Age text
        txt = str( int(np.round(p.star_age / 1e6, 0))) # Myr
        props = dict(boxstyle='round', facecolor='wheat', alpha=1)
        ax.text(0.05, 0.1, 'Age: ' + txt + ' Myr', 
                transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        ax.grid()
        # ax.set_yscale('log')
        # ax.set_xscale('log')
        
        
    fig.suptitle(title, fontsize = 18, y = 0.98)
    fig.text(0.47, -0.02, r'Radius $[R_\oplus]$', 
             fontsize = 20, transform = fig.transFigure)
    # fig.text(0.99, 0.4, 'Temperature [K]', rotation = 270, 
    #           fontsize = 15, transform = fig.transFigure)
    fig.text(-0.02, 0.4, r'log$\left(  \frac{dT}{dP} \right)$', rotation = 90, 
              fontsize = 20, transform = fig.transFigure)

    # Legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='k', lw = 3),
                    Line2D([0], [0], color = 'b', lw = 3),]
    
    props = dict(facecolor=back, alpha=1)
    fig.text(x = 0.4, y= -0.075, s ='Convective Regions', rotation='horizontal', 
          transform = fig.transFigure, horizontalalignment='center',
          bbox = props, fontsize=15)
 
    fig.legend(custom_lines, [r'$\nabla_a$', r'$\nabla_r$'],
                fontsize =  14, ncols = 3, alignment = 'left', # Lawful Neutral
                bbox_to_anchor=(0.87, -0.03), bbox_transform = fig.transFigure,)

    # Fill color

energy_transpport(path1, title1)
energy_transpport(path2, title2)
