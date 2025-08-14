#!/usr/bin/env python3

""" find_depths.py

Aim: Cleaner python 3 version of my depth codes, make this nice!

Created: Fri 29th Nov 2019
Edited: Fri 7th Feb 2025

"""
# debug sep import issue
import sys
print('sys.path', sys.path)

###################### Import useful modules #########################
import numpy as np
from new_depth_codes import get_depths
import os

######################### Set-up ####################################
#fields = ['XMM1', 'XMM2', 'XMM3']
#fields = ['CDFS1', 'CDFS2', 'CDFS3'] #'PSFHOMO_COSMOS']
fields = ['COSMOS']

# required aperture diameter
#reqAper = 'all' # or specific one
#reqFilters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z'] #'J', 'H', 'Ks'] #['CFHTW1-u', 'CFHT-iy', 'CFHT-u', 'CFHT-z', 'J'] #
reqFilters = ['HSC-NB0816_DR3']
#reqFilters = ['HSC-R_DR3', 'HSC-Z_DR3', 'HSC-NB0816_DR3', 'HSC-NB0921_DR3']
#reqFilters = ['J', 'H', 'Ks', 'NB118']
#reqFilters = ['J_DR4', 'H_DR4', 'Ks_DR4']
#reqFilters = ['NB118_DR4', 'NB118']
#reqFilters = ['J', 'J_DR4', 'J_DR6', 'H', 'H_DR4', 'H_DR6', 'K_DR4', 'K', 'K_DR6', 'Y', 'Y_DR4', 'Y_DR6']
reqFilters = ['HK']
#reqFilters = ['NB118_DR4', 'NB118', 'N_DR6']

# enter the prefered queue.
queue =  'none'
overwrite = False # False is the default anyway

####################################################################
# required aperture diameters to run through
# change at your peril!
# apDiametersAS = [1.8, 2.0, 3.0, 4.0, 5.0]
# for IRAC [2.8, 3.8, 5.8, 9.8, 11.6]
#apDiametersAS = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])
#apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0]) # Testing an issue with NB0816 DR3 in COSMOS

#####################
########## Loop ################################
## Loop through the different fields
for ff, fieldName in enumerate(fields):
    
    print('#############################################')
    print("Analysing field ", fieldName)
    outputDir = '../depths/{0}/'.format(fieldName)

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    # get the depths, nice clean code
    get_depths(fieldName, queue = queue, reqFilters = reqFilters, \
               overwrite = overwrite, outputDir = outputDir, apDiametersAS = apDiametersAS)
