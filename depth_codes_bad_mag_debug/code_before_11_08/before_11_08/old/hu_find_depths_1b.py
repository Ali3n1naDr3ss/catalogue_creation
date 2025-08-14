

""" hu_find_depths.py
17/07/25
Testing find_depths files.
"""

########## Import useful modules###############################################################################
import os
import sys
import numpy as np
from hu_depth_codes import get_depths, detections_fig, detections_fig_full, make_cutout

########## Set-up #############################################################################################
fields = ['COSMOS_test']# HOLLY 17/07/25 doing some tests
reqFilters = ['J','JH']
back_size = 32
back_filtersize = 9

queue =  'none'
overwrite = False # rename files with parameter values in title instead
verbose = False
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

    # make random, small cutout to test on
    large_images = []
    for filt in reqFilters:
        large_image = '/raid/scratch/data/COSMOS_test/UVISTA_' + filt + '_DR6.fits'
        large_images.append(large_image)
    output_img_name = 'none'
    #output_img_name, output_wht_name = make_cutout(large_images, verbose=verbose)

    # get the depths
    outputCatalogues = get_depths(fieldName, back_size, back_filtersize, queue=queue, reqFilters=reqFilters, \
               overwrite=overwrite, outputDir=outputDir, apDiametersAS=apDiametersAS)
    
    if output_img_name == 'none':
        print(">>>>>>>>>>>>>>>>>>>>>>>>", outputCatalogues)
        detections_fig_full('/raid/scratch/data/COSMOS_test/UVISTA_J_DR6.fits', ['J','JH'], outputCatalogues) # warning: testing on J, JH only        
    else:
        detections_fig(output_img_name, ['J','JH'], outputCatalogues) # warning: testing on J, JH only




