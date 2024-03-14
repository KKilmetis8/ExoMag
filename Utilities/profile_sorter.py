#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 15:09:51 2024

@author: konstantinos
"""
import numpy as np
def profile_sorter(profiles):
    ''' yes this is N^2, yes I could make this better, this is not a problem'''
    zs = np.arange(1, len(profiles)+1) # Zahlen!
    new_profiles = []
    for z in zs: 
        z = str(z)
        z = 'e' + z + '.' # don't mix 10 with 1
        profile = [s for s in profiles if z in s] 
        new_profiles.append(profile[0]) # only 1 match idealy
    return new_profiles