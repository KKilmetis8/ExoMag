#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:43:40 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
from src.Extractors.time_extractor import linear_fit_days
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['axes.facecolor']='whitesmoke'
plt.rcParams['figure.figsize'] = [6.0, 3.0]
plt.rcParams['font.family'] = 'Times New Roman'
lei_blue = '#001158'
lei_red = '#be1908'
lei_green = '#aaad00'
AEK = '#F1C410' # Important color

change = True
if change:
    lei_blue = AEK
    lei_red = 'k'
    lei_green = 'maroon'
# Define plot
fig, ax = plt.subplots(1,2, tight_layout = True)
fig.suptitle('Orbital Energy Distribution', y=0.92)
def circ(e):
    return np.sqrt(1 - e**2)
def flat(e):
    return np.sqrt(1 - e**40)

test = np.linspace(-1, 1, num = 2000)
ax[0].plot(test, circ(test), 
           color = lei_blue, alpha = 0.9, label = 'Circular')
ax[0].plot(test, flat(test), 
           color = lei_red, alpha = 0.9, label = 'Flattish')
ax[0].hlines(0.5, 0, 1, color = lei_green)
ax[0].text(0.66, 0.35, r'$\Delta E$',
          color = lei_green,
          fontsize = 17,
          transform = ax[0].transAxes)

ax[0].vlines(0, -1, 1, color = 'k')
ax[0].hlines(0, -2, 2, color = 'k')
ax[0].set_xlim(-1.1, 1.1)
ax[0].set_ylim(-0.1, 1.1)
ax[0].grid()
ax[0].legend(loc = 'lower center')
ax[0].set_title('Theoretical')
ax[0].set_xlabel('Specific Energy')
ax[0].set_ylabel('dM/dE')

#%% Simulation Hist
# Conversion constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2
t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13

# Data Load
fix = '820'
fixdays = str(linear_fit_days(int(fix)))
mass = np.load(fix + '/Mass_' + fix + '.npy')
X = np.load(fix + '/CMx_' + fix + '.npy')
Y = np.load(fix + '/CMy_' + fix + '.npy')
Z = np.load(fix + '/CMz_' + fix + '.npy')
Vx = np.load(fix + '/Vx_' + fix + '.npy')
Vy = np.load(fix + '/Vy_' + fix + '.npy')
Vz = np.load(fix + '/Vz_' + fix + '.npy')
R = np.sqrt( np.power(X,2) + np.power(Y,2)+ np.power(Z,2))
V = np.sqrt( np.power(Vx,2) + np.power(Vy,2)+ np.power(Vz,2))
Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
Bound = np.where(Orbital < 0, 1, 0)
# Spec. energy for bound
spec_bound =  Orbital*Bound # *mass # The minus is me cheating

# Get unbound -1 turns 0 -> -1 and 1->0. 
# For a weird potential reason, 1008 gets negative energies
Unbound = (Bound - 1) * (-1)

# Spec. energy for unbound
spec_unbound = Orbital * Unbound # * mass

# Mask extreme values.
masking = True
if masking:
    lower_bound = -100 # 1000 for specific, 0.75e-5 for simunits
    upper_bound = 100
    bound_mask = np.logical_and(spec_bound>=lower_bound,
                                spec_bound<upper_bound)
    unbound_mask = np.logical_and(spec_unbound>=lower_bound,
                                  spec_unbound<upper_bound)
    mask = np.logical_and(bound_mask, unbound_mask) # Combine the 2 masks
    
    mass = mass[mask]
    Bound = Bound[mask]
    Unbound = Unbound[mask]
    spec_bound = spec_bound[mask]
    spec_unbound = spec_unbound[mask]

# Calculate bound and unbound mass
mass_bound = Bound*mass
mass_unbound = np.abs(Unbound*mass)

#%% Plot
dm, de, _ = ax[1].hist( spec_bound, bins=50, # BEEG PAPATZILIKI
          weights = mass_bound, 
          color = lei_blue,
          label = 'Bound',
          cumulative = False,
          log=True,
          alpha = 0.9)
ax[1].hist(  spec_unbound, bins=50,
          weights =  mass_unbound, 
          color = lei_red,
          label ='Unbound',
          cumulative = False,
          log=True,
          alpha = 0.9)

# Plotting boilerplate
ax[1].legend(loc = 'lower center')
ax[1].set_title('Simulated')
ax[1].set_xlabel(r'Specific Energy', fontsize=10)
ax[1].set_ylabel(r'Mass $\cdot$ Counts $\left[M_{\odot} \cdot \# \right]$', 
           fontsize=10)
ax[1].grid()