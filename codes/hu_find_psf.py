#!/usr/bin/env python3

""" find_psf.py

Aim: cleaner version of code run_psf.py from before.  Batch run this like the depths.

Created: Wed 22nd Jan 2020
Modified: 7th Feb 2025

"""

###################### Import useful modules #########################
import os
import sys
import numpy as np
from hu_psf_codes import get_psf
######################### Set-up ####################################
baseDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/'

fields = ['COSMOS']

# required aperture diameter
reqFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']
queue = 'none'
overwrite = False # <<<<<<<<<<<<<<<<< OVERWRITE

# to run the first stage for each filter, set to True
# then check the pdf, update the star_param file
# then run again set to False
stars = False

############################### Loop ################################
## Loop through the different fields
for ff, fieldName in enumerate(fields):
    
    print('#############################################')
    print("Analysing field ", fieldName)

    outputDir = os.path.join(baseDir ,'psf/')
    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    outputDir = os.path.join(baseDir ,'psf/{0}/'.format(fieldName))
    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    # get the depths, nice clean code
    get_psf(fieldName, queue = queue, reqFilters = reqFilters, \
               overwrite = overwrite, outputDir = outputDir, starsOnly = stars)



