

""" hu_find_depths.py
17/07/25
Testing find_depths files.
"""

########## Import useful modules###############################################################################

import os
import sys
import numpy as np
from hu_depth_codes import get_depths, detections_fig, detections_fig_full, make_cutout, bkg_plotter, get_random_coords, open_cats, show_in_ds9
#from hu_depth_codes_1_10 import get_depths, detections_fig, detections_fig_full, make_cutout, bkg_plotter, get_random_coords

########## Set-up #############################################################################################
testing = False

reqFilters = ['K']
queue =  'none'

if testing:
    fields = ['COSMOS_test']
    size_arcsec = 2000
    coords = 'random'
    random_coords_file = '/raid/scratch/data/data_test/COSMOS_test/random_cutouts/random_coords.txt'
    catDir ='/raid/scratch/depths/COSMOS_test/catalogues/'
    back_size = 32
    back_filtersize = 9
    overwrite = False # rename files with parameter values in title instead
    verbose = False
else:
    # uses _backup.fits from data_test/COSMOS_test ~1178
    # using data/COSMOS _backup.fits (original files)
    fields = ['COSMOS_test']
    catDir ='/raid/scratch/depths/COSMOS_test/catalogues/'
    back_size = 32
    back_filtersize = 9
    coords = 'full-size'
    overwrite = False
    verbose = False

########## aperture diameters ###############################################################################

# required aperture diameters to run through
# change at your peril!
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])

########## Run ###############################################################################################

for ff, fieldName in enumerate(fields):

    outputDir = '/raid/scratch/depths/{0}/'.format(fieldName) 
    print("OutputDir:  ", outputDir)

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)

    if testing:
        # make random, small cutout to test on
        large_images = []
        for filt in reqFilters:
            large_image = '/raid/scratch/data/data_test/COSMOS_test/UVISTA_' + filt + '_DR6_backup.fits'
            large_images.append(large_image) # data and weight imgs
        output_img_name, output_wht_name = make_cutout(large_images, size_arcsec=size_arcsec, ra='random', dec='random', corner='lower-left', output_path=None, verbose=verbose)
        
        # get set of random coords from .txt to test on
        coords = get_random_coords(random_coords_file)

    # get the depths
    outputCatalogues = get_depths(fieldName, back_size, back_filtersize, coords, queue=queue, reqFilters=reqFilters, overwrite=overwrite, dataDir='/raid/scratch/data/', outputDir=outputDir, apDiametersAS=apDiametersAS)

    #outputCatalogues =  ['/raid/scratch/depths/COSMOS_test/catalogues/dJbksz32bkfilt9_full-size.fits', '/raid/scratch/depths/COSMOS_test/catalogues/dJHbksz32bkfilt9_full-size.fits']

    if testing:
        matching_coords = detections_fig(cutoutPath=output_img_name, detPath=catDir, filters=['J','JH'], detCats=outputCatalogues) # warning: testing on J, JH only
        ra = matching_coords.split('_')[0]
        dec = matching_coords.split('_')[1]

        open_cats(catDir, reqFilters, size_arcsec=size_arcsec, bksz=back_size, bkfilt=back_filtersize, ra=ra, dec=dec, show_all=False)

        #show_in_ds9(coords=coords, param_combo=[back_size, back_filtersize])
    
        bkg_plotter([back_size, back_filtersize], coords=matching_coords, verbose=verbose)
    else:
        output_img_name = 'none'
        #open_cats(catDir, reqFilters, size_arcsec='full-size', bksz=back_size, bkfilt=back_filtersize, show_all=False)

        #show_in_ds9(param_combo=[back_size, back_filtersize], filters=reqFilters)

        #detections_fig_full('/raid/scratch/data/data_test/COSMOS_test/UVISTA_J_DR6_backup.fits', ['J','JH'], outputCatalogues) # warning: testing on J, JH only #TODO

        #bkg_plotter([back_size, back_filtersize], verbose=verbose) #TODO



