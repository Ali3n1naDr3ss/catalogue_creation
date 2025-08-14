

""" hu_find_depths.py
17/07/25
Testing find_depths files.
"""

import sys
print('sys.path: ', sys.path)

########## Import useful modules###############################################################################
import numpy as np
from hu_depth_codes_1_4 import get_depths, detections_fig
from hu_depth_codes_1_3 import get_depths
import os

########## Set-up #############################################################################################
fields = ['COSMOS_test']
reqFilters = ['J','JH'] # NOTE: full JH_DR6 and weight images have been renamed to _backup.fits. The current JH files are cut-outs 30/05/25

queue =  'none'
overwrite = False # rename files with parameter values in title instead

########## aperture diameters ###############################################################################
# required aperture diameters to run through
# change at your peril!
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])

########## Run ###############################################################################################
for ff, fieldName in enumerate(fields):
    
    print('##########################################################################################')
    #outputDir = '../depths/{0}/'.format(fieldName)
    outputDir = '../depths/{0}/'.format(fieldName)
    print("OutputDir:  ", outputDir)

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    # get the depths
    outputCatalogues = get_depths(fieldName, queue = queue, reqFilters = reqFilters, \
               overwrite = overwrite, outputDir = outputDir, apDiametersAS = apDiametersAS)

    #detections_fig(['J','JH'], outputCatalogues) # warning: testing on J, JH only




