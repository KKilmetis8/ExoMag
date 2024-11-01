#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 17:40:46 2024

@author: konstantinos

Figure 5
"""

# Vanilla
import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
from tqdm import tqdm
from scipy.interpolate import griddata, RBFInterpolator
from scipy.optimize import curve_fit
import colorcet
import matplotlib.colors as colors
# Choc
import src.prelude as c
from src.Bfield.hardB import hardB_doer_single

class apothicarios:
    def __init__(self):
        self.masses = []
        self.fenvs = []
        self.orbsep = []
        self.Bdyn = []
        self.Bdip = []
        self.radii = []
        self.peaks = []
    def __call__(self, m, env, orbsep, dyn, dip, r, peak):
        self.masses.append(m)
        self.fenvs.append(env)
        self.orbsep.append(orbsep)
        self.Bdyn.append(dyn)
        self.Bdip.append(dip)
        self.radii.append(r)
        self.peaks.append(peak)

#%% Calc
age = 1500
jups_close = apothicarios()
jups_mid = apothicarios()
jups_far = apothicarios()
neps_close = apothicarios()
neps_mid = apothicarios()
neps_far = apothicarios()

path = f'data/{age}profile'
models = os.listdir(path)
models_good = []
# pick good fenvs
good_fenvs_j = np.array([0.86, 0.88, 0.90, 0.92, 0.94])
good_fenvs_n = np.array([0.02, 0.04, 0.06, 0.08, 0.1])

# Filter by filesize
for model in models:
    size = os.path.getsize(f'{path}/{model}')
    if size < 500_000:
        models_good.append(model)
        
for modelname in tqdm(models_good):
    if 'env' in modelname:
        underscores = [ i for i, x in enumerate(modelname) if x == '_']
        m = int( modelname[ 1: underscores[0] ])
        escape = modelname[ underscores[1]+1:underscores[2] ]
        a = float( modelname[ underscores[2]+2:underscores[3]])
        fenv = modelname[ underscores[0]+1: underscores[1] ].replace('env', '')
        fenv = fenv.replace('e', '')
        fenv = fenv.replace('_', '')
        fenv = float(fenv)
        try:
            B_dyn, B_dip, age_this = hardB_doer_single(path+f'/{modelname}')
            peak = 2.8 * B_dip # MHz 1e6 # Hz
            if np.abs(age-age_this) > 500:
                # print(m)
                # print(age_this)
                continue
            p = mr.MesaData(path+f'/{modelname}')
            r = np.power(10, p.logR)[0] *  c.Rsol / c.Rearth
        except ValueError or FileNotFoundError:
            continue
        if  B_dyn>0 and escape == 'zero':
            if m > 150 and fenv in good_fenvs_j:
                if a == 0.05:
                    jups_close(m, fenv, a,  B_dyn, B_dip, r, peak)
                elif a == 0.1:
                    jups_mid(m, fenv, a,  B_dyn, B_dip, r, peak)
                elif a == 0.2: 
                    jups_far(m, fenv, a,  B_dyn, B_dip, r, peak)
            elif m < 30 and  B_dip<15 and fenv in good_fenvs_n:
                if a == 0.05:
                    neps_close(m, fenv, a,  B_dyn, B_dip, r, peak)
                elif a == 0.1:
                    neps_mid(m, fenv, a,  B_dyn, B_dip, r, peak)
                elif a == 0.2: 
                    neps_far(m, fenv, a,  B_dyn, B_dip, r, peak)
#%% Mass - Bdip
#fig, ax = plt.subplots(2,3, figsize = (6,6), tight_layout = True,) 
fig = plt.figure( figsize=(6,6))
ax1 = plt.subplot2grid((2,6), (0, 0), colspan = 2)
ax2 = plt.subplot2grid((2,6), (0, 2), colspan = 2)
ax3 = plt.subplot2grid((2,6), (0, 4), colspan = 2)
ax4 = plt.subplot2grid((2,6), (1, 0), colspan = 3)
ax5 = plt.subplot2grid((2,6), (1, 3), colspan = 3)

jupss = [jups_close, jups_mid, jups_far]
nepss = [neps_mid, neps_far]

def power_law(x,a,b):
    return x**a + b

def fitter(planets, fenv_target):
    b = []
    m = []
    for i, fenv in enumerate(planets.fenvs):
        if fenv == fenv_target:
            b.append(planets.Bdip[i])
            m.append(planets.masses[i])
    fit, _ = curve_fit(power_law, m, b)   
    xmass = np.linspace(7, 27, num = 1000)
    # plt.figure()
    # plt.plot(xmass, xmass**fit[0] + fit[1], c = 'r')
    # plt.scatter(m,b,c='k')
    return fit 

def plotter(planetss, kind, ax):
    if kind == 'jups':
        k = 0
        xmass_max = np.linspace(170, 500, num = 1000)
        xmass_min = np.linspace(170, 500, num = 1000)
        bounds = good_fenvs_j
        ymin = 70
        ymax = 250
        minfenv = np.min(good_fenvs_j)
        maxfenv = np.max(good_fenvs_j)
        cmin = minfenv * 100
        cmax = maxfenv * 100
        cmap = 'Oranges'
        cline = 'cyan'
        ax[1].set_xlabel(r'Mass [M$_\mathrm{\oplus}$]', fontsize = 14)
        x = 0.13
        cax = fig.add_axes([x, 0.99, 0.87 - x, 0.03])
        cb_label_pos = 'top'
        boundaries = good_fenvs_j * 100

    if kind == 'neps':
        k = 1 
        xmass_max = np.linspace(5, 28, num = 1000)
        xmass_min = np.linspace(17, 28, num = 1000)
        ymin = 0
        ymax = 13
        minfenv = 0.02 # np.min(good_fenvs_n)
        maxfenv = 0.1 # good_fenvs_n[-3]# np.max()
        cmin =  np.min(good_fenvs_n) * 100
        cmax =  np.max(good_fenvs_n) * 100
        cmap = 'cet_CET_L12'
        cline = 'r'
        x = 0.13
        cax = fig.add_axes([x, -0.07, 0.87 - x, 0.03])
        boundaries = good_fenvs_n  * 100
        cb_label_pos = 'bottom'

    for i, planets in enumerate(planetss):
        mass = np.array(planets.masses)
        Bdip = np.array(planets.Bdip)
        fenv = np.array(planets.fenvs) * 100
        
        import matplotlib as mpl        
        img = ax[i].scatter(mass, Bdip, c = fenv,
                        edgecolor = 'k', cmap = mpl.colormaps[cmap].resampled(5),
                        vmin = cmin-1, vmax = cmax+1)
        
        # Peak axis
        ax_t = ax[i].twinx()
        ticklabels = np.array(ax[i].get_yticks())
        ax_t.set_yticks(ticklabels * 2.8)

        # Limits
        ax[i].set_ylim(ymin, ymax)
        ax_t.set_ylim(ymin * 2.8, ymax * 2.8)

        # Line
        max_fit = fitter(planets, maxfenv)
        y_max_fit = xmass_max**max_fit[0] + max_fit[1]
        ax[i].plot(xmass_max, y_max_fit, c = cline, ls = '--', lw = 2)
        min_fit = fitter(planets, minfenv)
        y_min_fit = xmass_min**min_fit[0] + min_fit[1]
        ax[i].plot(xmass_min, y_min_fit, c = cline, ls = '--', lw = 2)
        
        # Tick labels
        if i>0:
            ax[i].set_yticklabels([])
        if i != len(planetss)-1:
            ax_t.set_yticklabels([])
        else:
            ax_t.set_ylabel('Peak Frequency [MHz]', fontsize = 14, labelpad = 10)

        
        # Text
        lawmax = f'{round(max_fit[0],2)}'
        lawmin = f'{round(min_fit[0],2)}'
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax[i].text(0.05, 0.9, f'$B_{{dip}}^{{max}} \propto M_p^{{{lawmax}}}$',
                      fontsize = 10, c = 'k', transform = ax[i].transAxes,
                      bbox = props)
        ax[i].text(0.95, 0.1, f'$B_{{dip}}^{{max}} \propto M_p^{{{lawmin}}}$',
                      horizontalalignment = 'right',
                      fontsize = 10, c = 'k', transform = ax[i].transAxes,
                      bbox = props)
        
    # Colorbar
    cb = fig.colorbar(img, cax = cax,
                      orientation = 'horizontal' )
    cb.set_ticks(boundaries)
    cax.xaxis.set_ticks_position(cb_label_pos)
    cax.xaxis.set_label_position(cb_label_pos)
    cb.set_label('Envelope Mass Fraction [$\%$]', fontsize = 14, labelpad = 10)
    ax[0].set_ylabel(r'B$_\mathrm{dip}^\mathrm{(max)}$ [G]', fontsize = 14)
    
# Labels
x_label_pos = ax2.xaxis.get_label()
x = x_label_pos.get_position()
fig.text(x[0], -1.38, r'Mass [M$_\mathrm{\oplus}$]', fontsize = 14,
         horizontalalignment = 'center',
         transform = ax2.transAxes)

plotter(jupss, 'jups', [ax1, ax2, ax3])
plotter(nepss, 'neps', [ax4, ax5] )
plt.tight_layout()
