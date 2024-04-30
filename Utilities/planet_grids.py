#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:45:34 2024

@author: konstantinos
"""
import numpy as np
# Nep grid
masses = [4, 4, 4, 4, 4, 4,
          6, 6, 6, 6, 6, 6,
          8, 8, 8, 8, 8, 8, 
          10, 10, 10, 10, 10, 10,
          12, 12, 12, 12, 12, 12,
          14, 14, 14, 14, 14, 14,
          17, 17, 17, 17, 17, 17,
          20, 20, 20, 20, 20, 20,
          22, 22, 22, 22, 22, 22,
          24, 24, 24, 24, 24, 24,
          26, 26, 26, 26, 26, 26,
          28, 28, 28, 28, 28, 28,
          30, 30, 30, 30, 30, 30,
          ]
fenvs = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12] * (len(masses) // 6)
n_masses = list(np.unique(masses))
n_envs = list(np.unique(fenvs))

# Sat grid -> Pain
# masses = [60, 60, 60, 60, 60,
#           70, 70, 70, 70, 70,
#           80, 80, 80, 80, 80,
#           90, 90, 90, 90, 90,
#           100, 100, 100]

# Jup grid
masses = [  180, 180, 180, 180, 180, 180,
            190, 190, 190, 190, 190, 190,
            200, 200, 200, 200, 200, 200,
            210, 210, 210, 210, 210, 210,
            220, 220, 220, 220, 220, 220,
            230, 230, 230, 230, 230, 230,
            240, 240, 240, 240, 240, 240,
            250, 250, 250, 250, 250, 250,
            260, 260, 260, 260, 260, 260,
            270, 270, 270, 270, 270, 270,
            280, 280, 280, 280, 280, 280,
            290, 290, 290, 290, 290, 290,
            300, 300, 300, 300, 300, 300,
            310, 310, 310, 310, 310, 310,
            317, 317, 317, 317, 317, 317,
            325, 325, 325, 325, 325, 325,
            330, 330, 330, 330, 330, 330,
            340, 340, 340, 340, 340, 340,
            360, 360, 360, 360, 360, 360,
            380, 380, 380, 380, 380, 380,
            400, 400, 400, 400, 400, 400,
            425, 425, 425, 425, 425, 425,
            450, 450, 450, 450, 450, 450,
            475, 475, 475, 475, 475, 475,
            500, 500, 500, 500, 500, 500,
            ]
fenvs = [0.85, 0.88, 0.90, 0.92, 0.94, 0.96] * (len(masses) // 6)
j_masses = list(np.unique(masses))
j_envs = list(np.unique(fenvs))