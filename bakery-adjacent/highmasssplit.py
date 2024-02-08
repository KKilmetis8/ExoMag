#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 17:39:00 2023

@author: konstantinos
"""

# %%########################################################################################################
# split for light/heavy planets

   if mp_me <= 25.:
      #########################################################################################################
      if (not os.path.exists(coremodel)) and do_put_in_core:
          inlist2 = "tmp_inlists/inlist_2_core_" + \
              str(mcore/mearth)[0:5] + "_ME_" + str(mp/mearth)[0:6] + "_ME"

          run_time = inlist.put_core_in_planet(
              mesh_delta_coeff, mcore, rhocore, inlist2,   createmodel, coremodel)
          success = True
          if not os.path.exists(coremodel):
              success = False
              mesh_delta_coeff = 0.0
              while success == False and mesh_delta_coeff <= 2:
                  mesh_delta_coeff = mesh_delta_coeff + 0.1
                  print('trying mesh_delta_coeff =', mesh_delta_coeff)
                  run_time = inlist.put_core_in_planet(
                      mesh_delta_coeff, mcore, rhocore, inlist2, createmodel, coremodel)
                  if os.path.exists(coremodel):
                      success = True
                      print('***mesh_delta_coeff ----->', mesh_delta_coeff)
              if success == False:
                  return OTPfail
          k = open('LOGS/history_2_core', 'r')
          for line in k.readlines():
              pass
          last_temp = line
          last = last_temp.split()
          print("final age in core setting =", str(
              10.**float(last[0]))[0:6], "year")
          # print "last[0]==1000",last[0]=="1000"

          print("step 2, Core inserted: ", success)

  #########################################################################################################
      if (not os.path.exists(relaxedmodel)) and do_relaxm:
          inlist3 = "tmp_inlists/inlist_3_reducem_" + \
              str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]

          run_time = inlist.relaxm(
              mesh_delta_coeff, mp, inlist3, coremodel, relaxedmodel)
          success = True
          if not os.path.exists(relaxedmodel):
              success = False
              mesh_delta_coeff = 0.0
              while success == False and mesh_delta_coeff <= 5:  # ACHTUNG
                  mesh_delta_coeff = mesh_delta_coeff + 0.1
                  print('trying mesh_delta_coeff =', mesh_delta_coeff)
                  run_time = inlist.relaxm(
                      mesh_delta_coeff, mp, inlist3, coremodel, relaxedmodel)
                  if os.path.exists(relaxedmodel):
                      success = True
                      print('***mesh_delta_coeff ----->', mesh_delta_coeff)
              if success == False:
                  return OTPfail
          k = open('LOGS/history_3_reducemass', 'r')
          for line in k.readlines():
              pass
          last_temp = line
          last = last_temp.split()
          print("final age in fat setting =", str(
              10.**float(last[0]))[0:6], "year")
          # print "last[0]==1000",last[0]=="1000"

          print("step 3, Atmosphere reduced to fat0 = ", fat_0, ": ", success)

    else:  # mp_me>26. !!!mcore>=20.*mearth, fat>=0.3 for lowest masses!
      #########################################################################################################

      coremodel1 = "p2_core_20.0_ME_.mod"
      relaxedmodel1 = "p3_relxM_" + str(mp/mearth)[0:6] + "_ME_core20.mod"

      # start from pre-grown core, reduce the mass to desired
      if (not os.path.exists(relaxedmodel1)) and do_relaxm:
          inlist3 = "tmp_inlists/inlist_3a_reducem_" + \
              str(mp/mearth)[0:6] + "_ME_fat" + str(fat_0)[0:6]

          run_time = inlist.relaxm(mesh_delta_coeff, mp, inlist3, coremodel1, relaxedmodel1)
          success = True
          if not os.path.exists(relaxedmodel1):
              success = False
              mesh_delta_coeff = 0.0
              while success == False and mesh_delta_coeff <=2:
                  mesh_delta_coeff = mesh_delta_coeff + 0.1
                  print('trying mesh_delta_coeff =', mesh_delta_coeff)
                  run_time = inlist.relaxm(mesh_delta_coeff, mp, inlist3, coremodel1, relaxedmodel1)
                  if os.path.exists(relaxedmodel1):
                      success = True
                      print('***mesh_delta_coeff ----->', mesh_delta_coeff)
              if success == False:
                  return OTPfail
          k = open('LOGS/history_3_reducemass', 'r')
          for line in k.readlines():
              pass
          last_temp = line
          last = last_temp.split()
          print("final age in fat setting =", str(
              10.**float(last[0]))[0:6], "year")
          # print "last[0]==1000",last[0]=="1000"

          print("step 3, Atmosphere reduced to fat0 = ", fat_0, ": ", success)
  #########################################################################################################
      if (not os.path.exists(relaxedmodel)) and do_put_in_core:
          inlist2 = "tmp_inlists/inlist_2a_core_" + \
              str(mcore/mearth)[0:5] + "_ME_" + str(mp/mearth)[0:6] + "_ME"
          # relaxedmodel1 = "p3_relxM_" + str(mp/mearth)[0:6] + "_ME_core20.mod"

          run_time = inlist.put_core_in_planet(
              mesh_delta_coeff, mcore, rhocore, inlist2, relaxedmodel1, relaxedmodel)
          success = True
          if not os.path.exists(relaxedmodel):
              success = False
              mesh_delta_coeff = 0.0
              while success == False and mesh_delta_coeff <=2:
                  mesh_delta_coeff = mesh_delta_coeff + 0.1
                  print('trying mesh_delta_coeff =', mesh_delta_coeff)
                  run_time = inlist.put_core_in_planet(
                      mesh_delta_coeff, mcore, rhocore, inlist2, relaxedmodel1, relaxedmodel)
                  if os.path.exists(relaxedmodel):
                      success = True
                      print('***mesh_delta_coeff ----->', mesh_delta_coeff)
              if success == False:
                  return OTPfail
          k = open('LOGS/history_2_core', 'r')
          for line in k.readlines():
              pass
          last_temp = line
          last = last_temp.split()
          print("final age in core setting =", str(
              10.**float(last[0]))[0:6], "year")
          # print "last[0]==1000",last[0]=="1000"

          print("step 2, Core inserted: ", success)

# end split
