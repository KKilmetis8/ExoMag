#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 10:26:06 2024

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import src.prelude as c
from src.Bfield.hardB import hardB_doer
#%


def peak(B):
    Deltaf = 2.8 * B # MHz 1e6 # Hz
    return Deltaf

def radio_doer(names):
    planets = hardB_doer(names)
    all_peaks = []
    for planet in planets:
        this_peak = np.array([ peak(B) for B in planet.Bdip]) # Hz
        all_peaks.append(this_peak)
    return planets, all_peaks

def plotter(planets, peaks, labels):
    colors = [c.c91, c.c92, c.c93, c.c94, c.c95, c.c96, c.c97, c.c99, 'k']
    fig, ax = plt.subplots(1,1, tight_layout = True, sharex = True,
                            figsize = (4,4))
    custom_lines = []
    for planet, this_peak, color, label in zip(planets, peaks, colors, labels):
        start = np.argmin(np.abs(np.array(planet.age) - 300)) # Myrs
        stop = len(planet.age)
        
        age = planet.age[start:stop]
        this_peak = this_peak[start:stop]
        ax.plot(age, this_peak, c = color, label = label) # MHz
        
        custom_lines.append( Line2D([0], [0], color = color, 
                                    label = label))
    ax.set_ylabel('Maximum Frequency [MHz]', fontsize = 15)
    ax.set_xlabel('Age [Myr]', fontsize = 15)
    ax.legend(loc = 'upper right', fontsize = 15)
    ax.set_xlim(300,)
    # fig.legend(custom_lines, labels,
    #         fontsize =  9, ncols = len(planets), alignment = 'left', # Lawful Neutral
    #         bbox_to_anchor=(0.94, -0.03), bbox_transform = fig.transFigure,)
if __name__ == '__main__':
    kind = 'nep_env'
    if kind == 'nep_env':
        name1 = 'm17_env0.02_zero_a0.1_s8'
        name2 = 'm17_env0.04_zero_a0.1_s8'
        name3 = 'm17_env0.06_zero_a0.1_s8'
        name4 = 'm17_env0.08_zero_a0.1_s8'
        name5 = 'm17_env0.1_zero_a0.1_s8'
        name6 = 'm17_env0.12_zero_a0.1_s8'
        names = [name1, name2, name3, name4, name5, name6]
        planets, peaks = radio_doer(names)
        labels = ['2', '4', '6', '8', '10', '12',]
        plotter(planets, peaks, labels)
