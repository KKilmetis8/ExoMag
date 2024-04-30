#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 18:08:32 2024

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 19:55:13 2024

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import src.prelude as c
from src.Bfield.hardB import hardB_doer
#%
def rescale_wind(orb_sep, T, M, vel_sw, vel_kep, num_den):
    ''' Give wind parameters at 1 AU, rescale them to the orbital separation
        of the planet. Assumes a parker wind. 
        orb_sep in AU.'''
        
    critical_vel = np.sqrt(c.ideal_gas * T / c.mu) # Isothermal wind, same T everywhere
    critical_point = c.G * 1.2 * c.Msol / (2 * critical_vel**2)
    critical_point *= c.AU
    if orb_sep > critical_point: # Grows as lnr
        scaler = np.sqrt(np.log(orb_sep)) 
        vel_sw_a = 2 * vel_sw * scaler
        num_a = (1/orb_sep)**2 * num_den/2 * 1/scaler
        rel_vel = vel_sw_a - vel_kep
    else: # Drops linearly
        slope = 2 * critical_vel**3 / (c.G * M)
        vel_sw_a = slope * orb_sep * (1/c.AU)
        num_a = (1/orb_sep)**2 * num_den/2 * vel_sw_a/vel_sw
        rel_vel = np.sqrt(vel_sw_a**2 + vel_kep**2)
        
    return num_a, rel_vel
        
def radio_power(number_den, vel, r_m):
    eta = 1e-5 # Efficiency Factor
    mp = 1.67e-24 # Proton mass
    rho = mp * number_den
    power = eta * rho * vel**3 * np.pi * r_m**2
    return power

def peak(B):
    Deltaf = 2.8 * B * 1e6 # Hz
    return Deltaf

def radio_flux(power, distance, peak, r, r_m):
    # Calculate cone
    delta_a = 17.5 * np.pi / 180
    a_0 = np.arcsin(np.sqrt(r/r_m))
    trig_stuff = np.cos(a_0 - delta_a/2) - np.cos(a_0 + delta_a/2)
    cone = 4 * np.pi * trig_stuff
    
    flux = power / (distance**2 * cone * peak)
    return flux

def radio_doer(names, num_den, vel, distance):
    planets = hardB_doer(names)[0]
    r_ms = np.array(planets.radius) * c.Rsol # r_m = Rp, cgs
    powers = np.array([radio_power(num_den, vel, r_m) for r_m in r_ms]) # 
    peaks = np.array([ peak(B) for B in planets.Bdip]) # Hz
    fluxes = np.array([ radio_flux(power, distance, peak, r_m, r_m) 
                for power, peak, r_m in zip(powers, peaks, r_ms)])
    fluxes *= c.jy * 1e6 # μJy
    return planets.age, powers, peaks, fluxes

def plotter(age, power, peaks, flux):
    start = np.argmin(np.abs(np.array(age) - 300)) # Myrs
    stop = len(age)
    
    age = age[start:stop]
    power = power[start:stop]
    peaks = peaks[start:stop]
    flux = flux[start:stop]

    fig, axs = plt.subplots(1,3, tight_layout = True, sharex = True,
                           figsize = (6,4))
    axs[0].plot(age, power , c = 'k') 
    axs[1].plot(age, peaks * 1e-6, c = 'k') # MHz
    axs[2].plot(age, flux, c = 'k')
    
    axs[0].set_ylabel('Radio Power [erg/s]', fontsize = 15)
    axs[1].set_ylabel('Maximum Frequency [MHz]', fontsize = 15)
    axs[2].set_ylabel(r'Radio Flux [$\mu$Jy]', fontsize = 15)

    axs[0].set_xlabel('Age [Myr]', fontsize = 15)
    axs[1].set_xlabel('Age [Myr]', fontsize = 15)
    axs[2].set_xlabel('Age [Myr]', fontsize = 15)

    #axs[0].set_xlim(300,)
    
    #axs[1].set_ylim(0, 1000)
    
    fig.suptitle('Saffar\'s Radio Emission', fontsize = 15, y=0.9)
    
if __name__ == '__main__':
    planet = 'Saffar Sun'
    if planet == 'Saffar':  
        names = ['Saffar']
        orb_sep = 0.06 # [AU]
        num_den = 23.3 # [cm-3]
        vel_sw = 400 * c.kms_to_cms # At 1 AU
        vel_kep = 55 * c.kms_to_cms # [cm]
        T = 100_000 # At 1 AU
        stellar_mass = 1.2 * c.Msol # [g]
        num_den, vel = rescale_wind(orb_sep, T, stellar_mass,
                                    vel_sw, vel_kep, num_den)
        distance = 44 * c.ly # [cm]
        age, power, peaks, flux = radio_doer(names, num_den, vel, distance)
        plotter(age, power, peaks, flux)
    if planet == 'Saffar Sun':  
        names = ['saffar_sun']
        orb_sep = 0.06 # [AU]
        num_den = 23.3 # [cm-3]
        vel_sw = 400 * c.kms_to_cms # At 1 AU
        vel_kep = 55 * c.kms_to_cms # [cm]
        T = 100_000 # At 1 AU
        stellar_mass = 1 * c.Msol # [g]
        num_den, vel = rescale_wind(orb_sep, T, stellar_mass,
                                    vel_sw, vel_kep, num_den)
        distance = 44 * c.ly  # [cm]
        age, power, peaks, flux = radio_doer(names, num_den, vel, distance)
        plotter(age, power, peaks, flux)