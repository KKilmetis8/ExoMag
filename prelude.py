#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 20:59:44 2024

@author: konstantinos
"""
import numpy as np

# Space Units
Rsol = 6.957e10 # [cm]
Rjup = 71.492e8 # [cm]
Rearth = 6.371e8 # [cm]
Msol = 1.98855e33 # [g]
Mjup = 1.8986e30 # [g]
Mjup_in_sol = Mjup / Msol
Mearth = 5.9722e27 # [g]
Lsol = 3.826e33 # [erg/s]

# Physics
c = 2.99792458e10 #[cm/s]
h = 6.62607015e-27 #[gcm^2/s]
Kb = 1.380649e-16 #[gcm^2/s^2K]
sigma = 5.670e-5 # [cgs]
sqrt2 = np.sqrt(2)

# Reynolds
critical_rey_mag_num = 50
inv_mu_0_cgs = c**2 / (4 * np.pi)
mu_0_SI = 4 * np.pi * 10**(-7)

# Hard B
porp_c = 0.63
# Colors
AEK = '#F1C410'
# 6 palette
darkb = '#264653'
cyan = '#2A9D8F'
prasinaki = '#6a994e'
yellow = '#E9C46A'
kroki = '#F2A261'
reddish = '#E76F51'

# 9 palette
c91 = '#e03524'
c92 = '#f07c12'
c93 = '#ffc200'
c94 = '#90bc1a'
c95 = '#21b534'
c96 = '#0095ac'
c97 = '#1f64ad'
c98 = '#4040a0'
c99 = '#903498'
# Plotting
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [6 , 5]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'