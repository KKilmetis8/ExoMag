#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 14:14:50 2024

@author: konstantinos

This opens every profile and saves the desired age to a seperate folder. needed for the big plots
"""
import numpy as np
import mesa_reader as mr
import os
from tqdm import tqdm
from src.Utilities.profile_sorter import profile_sorter

def age_finder(names, target = 1500, path = ''):
    for name in tqdm(names):
        profiles = os.popen(f'ls {path}/{name}/profile*.data').read()
        profiles = list(profiles.split("\n"))
        profiles.pop() # Remove last
        profile_numbers = [] 
        profiles = profile_sorter(profiles) # Guess what that does
        ages = []
        for profile, i in zip(profiles, range(len(profiles))):
            p = mr.MesaData(profile)
            age = np.round(p.star_age /  1e6, 0)
            ages.append(age)
            prof_no = profile[-8:-4].replace('.','').replace('e','').replace('l','')
            profile_numbers.append(prof_no)
            # print(age)
            # print(prof_no)
            # print(i)
            # print('---')
            if age > 500+target:
                break
        try:
            ages = np.array(ages)
            target_idx = np.argmin( np.abs(ages - target))
            target_profile_no = profile_numbers[target_idx]
        except:
            continue
        try:
            read = open(f'{path}/specific_ages/{target}.txt', 'r')
            if name not in read.read():
                read.close()
                append = open(f'{path}/specific_ages/{target}.txt', 'a')
                append.write(f'{name} {target_profile_no} \n')
                append.close()
        except FileNotFoundError:
            write = open(f'{path}/specific_ages/{target}.txt', 'w')
            write.write(f'{name} {target_profile_no} \n')
            write.close()

if __name__ == '__main__':
    kind = 'many'
    if kind == 'single':
        path = '/media/konstantinos/Dunkey/mesadata/lotsaneps/'
        m = 5
        fenv = 0.06
        escape = 'zero'
        orb_sep = 0.2
        s = 8
        model = f'm{m}_env{fenv}_{escape}_a{orb_sep}_s{s}'
        age_finder([model], path = path)
    elif kind == 'many':
        # path = '/media/konstantinos/Dunkey/mesadata/jups01'
        # path = 'data/'
        path = 'bakery/LOGS/neps0.1'
        models = os.listdir(path)
        age_finder(models, path = path)

                        
    
