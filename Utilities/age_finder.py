#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 14:14:50 2024

@author: konstantinos
"""
import numpy as np
import mesa_reader as mr
import os
import src.Utilities.planet_grids as grids


def age_finder(names):
    target = 10000 # Myrs
    for name in names:
        p_path = 'data/' + name
        profiles = os.popen('ls ' + p_path + '/profile*.data').read()
        profiles = list(profiles.split("\n"))
        profiles.pop() # Remove last
        profile_numbers = [] 
        # profiles = profile_sorter(profiles) # Guess what that does
        ages = []
        for profile, i in zip(profiles, range(len(profiles))):
            # Load data
            p = mr.MesaData(profile)
            age = np.round(p.star_age /  1e6, 0)
            # print(age)
            ages.append(age)
            prof_no = profile[-8:-4].replace('.','').replace('e','')
            profile_numbers.append(prof_no)
        ages = np.array(ages)
        target_idx = np.argmin( np.abs(ages - target))
        target_profile_no = profile_numbers[target_idx]
        print(ages[target_idx])
        print(target_profile_no)
    try:
        read = open(f'data/specific_ages/{target}.txt', 'r')
        if name not in read.read():
            read.close()
            append = open(f'data/specific_ages/{target}.txt', 'a')
            append.write(f'{name} {target_profile_no} \n')
            append.close()
    except FileNotFoundError:
        write = open(f'data/specific_ages/{target}.txt', 'w')
        write.write(f'{name} {target_profile_no} \n')
        write.close()

if __name__ == '__main__':
    #orb_seps1 = list(np.arange(0.07, 0.1, 0.001))
    orb_seps1 = list(np.arange(0.05, 0.3, 0.02))
    orb_seps2 = list(np.arange(0.3, 2, 0.2))
    orb_seps = orb_seps1 + orb_seps2
    # orb_seps = orb_seps2
    #orb_seps = list(np.arange(0.07, 0.1, 0.001))
    # orb_seps = [0.1, 0.2]
    # Nep grid
    masses = [317] # [95, 317, 17]
    fenvs = [0.94] # 0.9, 0.94, 0.06]
    ms = list(np.unique(masses))
    fenv = list(np.unique(fenvs))
    entropies = [8]
    for m in ms:
        for fenv in fenvs:
            #fenv = str(fenv).replace('.', '')
            for s in entropies:
                for orb_sep in orb_seps:
                    path = f'm{m}_env{fenv}_zero_a{orb_sep}_s{s}'
                    try:
                        age_finder([path])
                    except ValueError or FileNotFoundError:
                        continue
                    
    