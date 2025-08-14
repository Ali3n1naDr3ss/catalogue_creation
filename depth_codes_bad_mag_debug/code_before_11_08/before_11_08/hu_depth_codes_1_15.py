

""" hu_find_depths.py
17/07/25
Testing find_depths files.
"""

########## Import useful modules###############################################################################

import os
import sys
import numpy as np
from hu_depth_codes import get_depths, detections_fig, detections_fig_full, make_cutout, bkg_plotter, get_random_coords, show_in_ds9

########## Set-up #############################################################################################

fields = ['COSMOS_test']# HOLLY 17/07/25 doing some tests
reqFilters = ['J','JH']
back_size = 32
back_filtersize = 9
testing = True
#coords = 'full-size'
queue =  'none'
overwrite = False # rename files with parameter values in title instead
verbose = False
output_img_name = 'none'
random_coords_file = '/raid/scratch/data/COSMOS_test/random_cutouts/random_coords.txt'

########## aperture diameters ###############################################################################

# required aperture diameters to run through
# change at your peril!
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])

########## Run ###############################################################################################

for ff, fieldName in enumerate(fields):

    outputDir = '../depths/{0}/'.format(fieldName) 
    print("OutputDir:  ", outputDir)

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)

    if testing:
        # make random, small cutout to test on
        large_images = []
        for filt in reqFilters:
            large_image = '/raid/scratch/data/COSMOS_test/UVISTA_' + filt + '_DR6_backup.fits'
            large_images.append(large_image) # data and weight imgs
        output_img_name, output_wht_name = make_cutout(large_images, size_arcsec=500, verbose=verbose)
        
        # get set of random coords from .txt to test on
        coords = get_random_coords(random_coords_file)

    # get the depths
    outputCatalogues = get_depths(fieldName, back_size, back_filtersize, coords, queue=queue, reqFilters=reqFilters, overwrite=overwrite, outputDir=outputDir, apDiametersAS=apDiametersAS)
        

    if output_img_name == 'none':
        detections_fig_full('/raid/scratch/data/COSMOS_test/UVISTA_J_DR6.fits', ['J','JH'], outputCatalogues) # warning: testing on J, JH only 
        bkg_plotter([back_size, back_filtersize], coords, verbose=verbose)

    else:
        matching_coords = detections_fig(cutoutPath=output_img_name, detPath='/raid/scratch/depths/COSMOS_test/catalogues/', filters=['J','JH'], detCats=outputCatalogues) # warning: testing on J, JH only

        show_in_ds9(coords=coords, param_combo=[back_size, back_filtersize])
    
        bkg_plotter([back_size, back_filtersize], coords=matching_coords, verbose=verbose)



