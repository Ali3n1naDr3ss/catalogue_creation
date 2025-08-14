#!/usr/bin/env python3

""" find_psf.py

Aim: cleaner version of code run_psf.py from before.  Batch run this like the depths.

Created: Wed 22nd Jan 2020
Modified: 7th Feb 2025

"""
# debug sep import issue
import sys
print('sys.path', sys.path)
###################### Import useful modules #########################
import numpy as np
from new_psf_codes import get_psf
import os

######################### Set-up ####################################
#fields = ['XMM1', 'XMM2', 'XMM3']
#fields = ['CDFS1', 'CDFS2', 'CDFS3']
fields = ['COSMOS']

# required aperture diameter
#reqAper = 'all' # or specific one
#reqFilters = ['HSC-R_DR3', 'HSC-Z_DR3', 'HSC-NB0816_DR3', 'HSC-NB0921_DR3']
#reqFilters = ['Y', 'J', 'H', 'Ks']
reqFilters = ['HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'HSC-Y_DR3']
# enter the prefered queue.
queue = 'none'
overwrite = True # False is the default anyway

# to run the first stage for each filter, set to true
# then check the pdf, update the star_param file
# then run again set to false
stars = False

############################### Loop ################################
## Loop through the different fields
for ff, fieldName in enumerate(fields):
    
    print('#############################################')
    print("Analysing field ", fieldName)
    outputDir = '../psf/{0}/'.format(fieldName)

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    # get the depths, nice clean code
    get_psf(fieldName, queue = queue, reqFilters = reqFilters, \
               overwrite = overwrite, outputDir = outputDir, starsOnly = stars)
