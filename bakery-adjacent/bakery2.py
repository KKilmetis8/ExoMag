#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 16:35:51 2023

@author: konstantinos
"""

# Vanilla Imports
import math
import numpy as np
import os
wdir = os.getcwd()
if wdir[-5:] != 'daria':
    os.chdir('daria')
import shutil
import time
import random
import sys
# Custom Imports
import inlist_maker as inlist

#%% Constants

msun = 1.98855e33
mjup = 1.8986e30
mearth = 5.9722e27

rearth = 6.378e8
rjup = 7.14e9 
rsun = 6.9598e10

Lsun = 3.9e33
sigma=5.67e-5
au = 1.496e13
#%% Parameters 

# Planetary Parameters 
mp = 12.5*mearth	# Total planetary mass [g]
fat_0 = 0.3	# Init. atm. mass fraction in planetary masses 
maxage = 1.e10 # Final age [yr]

# Stellar. prameters
ms = 1. # Stellar Mass [Msol]
rs = 1. # Stellar radius in [Rsol]
Teff_star = 5780. # [K]
orb_sep = 0.05 #	orbital sepration [AU]
#BA= 0.0 # Bond albedo

# Atmospheric Escape
escape_model = 'HD'  # 'HD': HBA approximation (Kubyshkina 18b)
                     # 'EL': energy limited with eta = 0.15 
                     # (to change eta, change x_ctrl(2) in inlist_7)
                     # 'ELR': energy limited + radiation-recombination limited
                     #  (Chen&Rogers 16)
# Feuv = Cpow*age**(betapow), Ribas '05, solar values at 1 AU
escape_Cpow = 29.12  
escape_betapow = -1.24 
escape_sat = 8.85e3 # Feuv at saturation, 1 AU Kubyshkina 20

# Composition
z = 0.02 # Stellar & planetary metallicity 
y = 0.24 # (initial) Stellar & planetary He fraction

# Formation routine
mesh_delta_coeff = 0.5 # Grid points in MESA. 0.5 is a good init.guess, smaller for more points

rinitial = 3.0 # Inital planet radius in (formation routine) [rjup]
minitial= 25.6454 * mearth # For M_p < 25 M_earth, heavier are handled down
mcore = (1-fat_0)*mp	 # Mass [g]; not to change 
rhocore = inlist.ADEN_core( mp / mearth) # Density [g/cm^3]
# (uses rcore=mcore^0.27 approximation for numbers from Rogers et al., 2011)
mp_wo_core = mp - mcore	# Atmosphere mass [g]; not to change

maxEntropy = 9.	# Entopy of inflated atm. [kb/baryon] (init. thermal state; step iv in the paper)
coreLuminosity = 2.0e27	# backup value; ~works for ~4-30 ME; the true value is set in the setS procedure
coolingTime = 5.e6 # Cooling timescale [yr]
initage = 5e6 # Age when mass loss begins in years (the disc dispersal time)
irrad_col = 300.0 # Column depth for depositing stellar radiation as heat
flux_dayside = 0. # dayside flux set later in relax_irradiation (Lbol/4./pi/orb_sep[cm]**2)
# (only relevant for relaxing the initial irradiation)

#%% Flags

do_create_planet = True # Create a ball of gas
do_put_in_core = True 	# Put the core of the right size
do_relaxm = True # Remove the atmosphere to the desired value
do_set_entropy = True # Set artificial luminocity to inflate the planet
do_cool = True 	# Remove artificial luminocity and cool the planet
do_relax_irradiation = True 	# Stellar heating (Teq); 
do_evolve_planet = True # Evolve!

#%% MIST track for irradiation

Mssg_list = list(range(40,95,5))+list(range(92,132,2)) # [40:5:90 92:2:130];    #GRID OF STELLAR MASSES AVAILABLE in MIST ( 0.4..1.3 MSUN)
Mssg = np.array(Mssg_list)
n_tmp, = np.where(abs(ms*100-Mssg)==min(abs(ms*100-Mssg)))# take the closest available mass; these tracks change quite smooth, so it is accurate enough. However, one can also interpolate between two closest tracks (if you need assistance with it, please contact me)
n_tmp  = n_tmp[0]

if Mssg[n_tmp]>=100:
    nameEtrack = 'ev_inp/ei'+str(Mssg[n_tmp])+'.dat'
    nameEtrack1 = 'ev_inp/eho'+str(Mssg[n_tmp])+'.dat'
else:
    nameEtrack = 'ev_inp/ei0'+str(Mssg[n_tmp])+'.dat';
    nameEtrack1 = 'ev_inp/eho0'+str(Mssg[n_tmp])+'.dat';
		
Etrack = np.loadtxt(nameEtrack) # np.loadtxt('ev_inp/ei050.dat') 

#%% Naming

# Every case
createmodel = "p1_create_" + str(minitial/mearth)[0:5]+ "_ME"  + ".mod"
coremodel = "p2_core_" + str(mcore/mearth)[0:5] + "_ME_" + str(mp/mearth)[0:6] + "_ME" + ".mod"
relaxedmodel = "p3_relxdM_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]  + ".mod"
entropymodel = "p4_setS_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] + "_" +str(maxEntropy)  + ".mod"
removemodel = "p5_removeLc_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] + "_" +str(maxEntropy)  + ".mod"
relaxirradmodel = "p6_relaxirrad_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy)+"_d_"+ str(orb_sep)+ "_" +str(ms)  + "Msun.mod"
evolvedmodel = "p7_evolved_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy) +"_d_"+ str(orb_sep)+ "_" +str(ms)  + "Msun_Cpow_" +str(escape_Cpow)+ "_Fsat_" + str(escape_sat/1e4) +".mod"
# Big planets
coremodel1 = "p2_core_20.0_ME_.mod"
relaxedmodel1 = "p3_relxdM_" + str(mp/mearth)[0:6] + "_ME_core20.mod"

#%% Mass Split 

# Use normal seqeuence if <25 Mearth, else use custom order 
# 3 -> 2 -> 4 -> 5 -> 6 -> 7
# Makes things faster
shutil.copyfile("src/run_star_extras0.f","src/run_star_extras.f")
os.system('./clean')
if mp < 25 * mearth:
    # Do 1
    if do_create_planet:
        
    # Do 2
        inlist1 = "tmp_inlists/inlist_1_create_" + str(mp/mearth)[0:6] + "_ME"
        run_time = inlist.create_planet(minitial, y, z, inlist1, createmodel)
        success = True
        if not os.path.exists(createmodel):
            success=False	
        k=open('LOGS/history_1_create','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final model number in create=",last[0])
        #print "last[0]==1000",last[0]=="1000"
        if last[0]=="1000":
            success=False
        print( "step 1, Planet created: ", success)
    if do_put_in_core:
        inlist2 = "tmp_inlists/inlist_2_core_" + str(mcore/mearth)[0:5] + "_ME_" + str(mp/mearth)[0:6] + "_ME"
  
        run_time = inlist.put_core_in_planet(mesh_delta_coeff, mcore, rhocore, inlist2,   createmodel, coremodel)
        success = True
        if not os.path.exists(coremodel):
            success=False	
            mesh_delta_coeff = 0.0
            while success==False and mesh_delta_coeff<=2:
                mesh_delta_coeff = mesh_delta_coeff + 0.1
                print('trying mesh_delta_coeff =', mesh_delta_coeff)
                run_time = inlist.put_core_in_planet(mesh_delta_coeff, mcore, rhocore, inlist2, createmodel, coremodel)
                if os.path.exists(coremodel):
                     success = True
                     print('***mesh_delta_coeff ----->', mesh_delta_coeff)
        k=open('LOGS/history_2_core','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final age in core setting =", str(10.**float(last[0]))[0:6], "year")  
        print( "step 2, Core inserted: ", success)
    # Do 3
    if do_relaxm:
        inlist3 = "tmp_inlists/inlist_3_reducem_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]

        run_time = inlist.relaxm(mesh_delta_coeff, mp,inlist3,coremodel,relaxedmodel)
        success = True
        if not os.path.exists(relaxedmodel):
            success=False
            mesh_delta_coeff = 0.0
            while success==False and mesh_delta_coeff<=5: #ACHTUNG
                mesh_delta_coeff = mesh_delta_coeff + 0.1
                print('trying mesh_delta_coeff =', mesh_delta_coeff)
                run_time = inlist.relaxm(mesh_delta_coeff, mp,inlist3,coremodel,relaxedmodel)
                if os.path.exists(relaxedmodel):
                     success = True
                     print('***mesh_delta_coeff ----->', mesh_delta_coeff)	
        k=open('LOGS/history_3_reducemass','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final age in fat setting =", str(10.**float(last[0]))[0:6], "year")
        print( "step 3, Atmosphere reduced to fat0 = ", fat_0,": ", success)
        inlist3 = "tmp_inlists/inlist_3_reducem_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]
        run_time = inlist.relaxm(mesh_delta_coeff, mp,inlist3,
                                 coremodel,relaxedmodel)
        success = True
        if not os.path.exists(relaxedmodel):
            success=False	
        k=open('LOGS/history_3_reducemass','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final age in fat setting =", str(10.**float(last[0]))[0:6], "year")
        print( "step 3, Atmosphere reduced to fat0 = ", fat_0,": ", success)
else:   
    # Do 3
    # Start from pre-grown core, reduce the mass to desired
    if do_relaxm:
        inlist3 = "tmp_inlists/inlist_3a_reducem_" + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]
        
        run_time = inlist.relaxm(mesh_delta_coeff, mp,
                                 inlist3, coremodel1, relaxedmodel1)
        success = True
        if not os.path.exists(relaxedmodel1):
            success=False
            mesh_delta_coeff = 0.0
            while success==False and mesh_delta_coeff<=2:
                mesh_delta_coeff = mesh_delta_coeff + 0.1
                print('trying mesh_delta_coeff =', mesh_delta_coeff)
                run_time = inlist.relaxm(mesh_delta_coeff, mp, 
                                         inlist3,coremodel1,relaxedmodel1)
                if os.path.exists(relaxedmodel1):
                     success = True
                     print('***mesh_delta_coeff ----->', mesh_delta_coeff)
        k=open('LOGS/history_3_reducemass','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final age in fat setting =", str(10.**float(last[0]))[0:6], "year")
        #print "last[0]==1000",last[0]=="1000"
        print( "step 3, Atmosphere reduced to fat0 = ", fat_0,": ", success)
        
    # Do 2
    if do_put_in_core:
        inlist2 = "tmp_inlists/inlist_2a_core_" + str(mcore/mearth)[0:5] + "_ME_" + str(mp/mearth)[0:6] + "_ME"
  
        run_time = inlist.put_core_in_planet(mesh_delta_coeff, mcore, rhocore, 
                                             inlist2, relaxedmodel1, relaxedmodel)
        success = True
        if not os.path.exists(relaxedmodel):
            success=False	
            mesh_delta_coeff = 0.0
            while success==False and mesh_delta_coeff<=2:
                mesh_delta_coeff = mesh_delta_coeff + 0.1
                print('trying mesh_delta_coeff =', mesh_delta_coeff)
                run_time = inlist.put_core_in_planet(mesh_delta_coeff, mcore, rhocore, inlist2, relaxedmodel1, relaxedmodel)
                if os.path.exists(relaxedmodel):
                     success = True
                     print('***mesh_delta_coeff ----->', mesh_delta_coeff)
        k=open('LOGS/history_2_core','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final age in core setting =", str(10.**float(last[0]))[0:6], "year")
        #print "last[0]==1000",last[0]=="1000"
        print( "step 2, Core inserted: ", success)
#%%
# Do 4    
if do_set_entropy:
    with open('LOGS/history_3_reducemass', 'r') as f:
        for line in f:
            pass
        last = line.split()
		
    currentropy= float(last[7])

    if currentropy<float(maxEntropy):
        with open('LOGS/history_3_reducemass', 'r') as f:
            for line in f:
                pass
            last2 = line.split()
		
        if maxEntropy<13.6:
            coreLuminosity= 30*float(last2[1])*3.846e33
        print( coreLuminosity)
        inlist4 = "tmp_inlists/inlist_4_setS_"  + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] +"_"+ str(maxEntropy)
        run_time = inlist.set_initial_entropy(mesh_delta_coeff, maxEntropy,
                                              coreLuminosity,inlist4,
                                              relaxedmodel,entropymodel)
        success = True
        if not os.path.exists(entropymodel):
            success=False	
        k=open('LOGS/history_4_setS','r')
        for line in k.readlines():
            pass
        last_temp=line
        last=last_temp.split()
        print( "final entropy set =", last[7])
        print( "step 4, Initial entropy set to = ", maxEntropy,": ", success)
        
# Do 5
if do_cool:
    Lcore_low = mcore*5e-8 # Arbitrary, gets overwritten
    inlist5 = "tmp_inlists/inlist_5_removeLc_"  + str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6] +"_"+ str(maxEntropy)
    run_time = inlist.remove_Lcore(mesh_delta_coeff, coolingTime, 
                                   Lcore_low, inlist5,entropymodel,removemodel)
    success = True
    if not os.path.exists(removemodel):
        success=False	
    k=open('LOGS/history_5_removeLc','r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "final age in cooling =", str(10.**float(last[0]))[0:6], "year")
    print( "step 5, Artificial core dissipation removed: ", success)

# Do 6
if do_relax_irradiation:    
    Lbol = np.interp(coolingTime, Etrack[:,0], Etrack[:,1]);# [log10(Lsun)]
    Lbol = (10.**(Lbol))*Lsun; #[erg/s]
    
    Teff = np.interp(coolingTime, Etrack[:,0], Etrack[:,2]);# [log10(K)]
    Teff = 10.**(Teff); #[K]

    Rss  = np.interp(coolingTime, Etrack[:,0], Etrack[:,3]);# [log10(Rsun)]
    Rss  = 10.**(Rss); #[Rsun]

    Teq = Teff*((Rss*rsun/2./orb_sep/au)**0.5); #TEQ(Teff,Rss,d0,0);
    flux_dayside = Lbol/4.0/3.1416/orb_sep/orb_sep/au/au
    relaxage = coolingTime + 1e5
    initage = relaxage
    
    inlist6 = "tmp_inlists/inlist_6_relaxirrad_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy)+"_d_"+ str(orb_sep)+ "_" +str(ms)  + "Msun"
    run_time = inlist.relax_irradiation(mesh_delta_coeff, irrad_col,flux_dayside, relaxage, 
                                        inlist6, removemodel, relaxirradmodel)
    success = True
    if not os.path.exists(relaxirradmodel):
        success=False	
    k=open('LOGS/history_6_relaxsurfheat','r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "disk dissipation time set to", str(10.**float(last[0]))[0:6], "year")
    print( "step 6, Heating switched on: ", success)
  
# Do 7
if do_evolve_planet:
    inlist7 = "tmp_inlists/inlist_7_evolve_"  + str(mp/mearth)[0:6] + "_ME_fat_" + str(fat_0)[0:6] +"_S"+ str(maxEntropy)+"_d_"+ str(orb_sep)+ "_" +str(ms)  + "Msun"
    run_time = inlist.evolve_planet(mesh_delta_coeff, initage, orb_sep, maxage, 
                                    nameEtrack1, escape_model, escape_Cpow, 
                                    escape_betapow, escape_sat, inlist7, 
                                    relaxirradmodel, evolvedmodel)
    success = True
    if not os.path.exists(evolvedmodel):
        success=False	
    k=open('LOGS/history' + evolvedmodel[6:-4],'r')
    for line in k.readlines():
        pass
    last_temp=line
    last=last_temp.split()
    print( "the planet has evolved for ", str((10.**float(last[0]))/1e9)[0:6], "Gyr")
    print( "step 7, Evolution complete: ", success)
    # print( "atmospheric escape model: ", escape_model)
    # print( "Feuv = ", escape_Cpow, "* (age/1 Gyr)^(", escape_betapow, ")")
    # print( "Feuv_sat = ", str(escape_sat), "erg/cm**2/s")
    # print( "evolutionary profiles saved in ", '"LOGS/history' + evolvedmodel[6:-4] + '"')

# Finishing Sound
from src.Utilities.finished import finished
finished()


