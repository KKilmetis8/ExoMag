#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 14:14:50 2024

@author: konstantinos
"""
import numpy as np
import mesa_reader as mr
import os

def age_finder(names):
    target = 5500 # Myrs
    for name in names:
        p_path = 'data/' + name
        profiles = os.popen('ls ' + p_path + '/profile*.data').read()
        profiles = list(profiles.split("\n"))
        profiles.pop() # Remove last
        # profiles = profile_sorter(profiles) # Guess what that does
        ages = []
        for profile, i in zip(profiles, range(len(profiles))):
            # Load data
            p = mr.MesaData(profile)
            age = np.round(p.star_age /  1e6, 0)
            # print(age)
            ages.append(age)
        ages = np.array(ages)
        target_idx = np.argmin( np.abs(ages - target)) + 1
    try:
        read = open(f'data/specific_ages/{target}.txt', 'r')
        if name not in read.read():
            read.close()
            append = open(f'data/specific_ages/{target}.txt', 'a')
            append.write(f'{name} {target_idx} \n')
            append.close()
    except FileNotFoundError:
        write = open(f'data/specific_ages/{target}.txt', 'w')
        write.write(f'{name} {target_idx} \n')
        write.close()

if __name__ == '__main__':
    ms = [8, 14, 17, 20, 25, 30, 40, 50, 60, 75, 150, 200, 250, 300, 350, 400]
    fenvs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 60, 70, 80, 85, 90, 92, 94, 95, 96, 98]
    entropies = [7, 8, 9]
    for m in ms:
        for fenv in fenvs:
            for s in entropies:
                path = f'm{m}_e{fenv}_zero_a01_s{s}'
                try:
                    age_finder([path])
                except ValueError or FileNotFoundError:
                    continue
                    
    