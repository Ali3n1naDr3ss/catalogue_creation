# visually compare the J, H, JH, and HK fluxes of a small sample of "bad mag" objects and store a Table with these objects
# 31/01/2025

############################### IMPORTS #######################################################
from functions_compare_fluxes_across_images import make_cutout, get_depths, show_in_ds9, open_cats, detections_fig
import os
import numpy as np 

############################### USER PARAMETERS #######################################################
baseDir = '/raid/scratch/'
dataDir = baseDir + 'data/data_test/COSMOS_flux_comp/'
subDir = 'flux_comparison_test/'
depthDir = 'depths/COSMOS_flux_comp/'
filters = ['J', 'H', 'JH', 'HK']
filters = ['J', 'JH']
size_arcsec = 1200
ra = 150.015
dec = 2.261




fieldName = 'COSMOS_flux_comp'
coords = [ra, dec]
back_size = 32
back_filtersize = 9
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])

verbose = False
save = False


################# FUNCTION CALLS - DO NOT MAKE CHANGES BELOW THIS LINE ##################
outputDir = '../depths/{0}/'.format(fieldName)
if os.path.isdir(outputDir) == False:
   os.system('mkdir ' + outputDir)

output_img_name, output_wht_name = make_cutout(dataDir, filters, subDir, size_arcsec, ra, dec)

outputCatalogues = get_depths(back_size=back_size, back_filtersize=back_filtersize, coords=coords, dataDir=dataDir, reqFilters=filters, outputDir=outputDir, apDiametersAS=apDiametersAS, size_arcsec=size_arcsec)

# get aperture fluxes table - then we can narrow it down into just a few objects and look at the numerical data (neg fluxes?)
show_in_ds9(coords=coords, filters=filters, dataDir=dataDir, size_arcsec=size_arcsec, depthDir=depthDir+'images/')

open_cats(baseDir+depthDir+'catalogues/', filters, size_arcsec=size_arcsec, show_all=False, bksz=back_size, bkfilt=back_filtersize, ra=ra, dec=dec)

detections_fig(dataDir+subDir, baseDir+depthDir+'catalogues/', filters, outputCatalogues, size_to_match=size_arcsec, cutout_coords_file='flux_comparison_test.txt', subsetCriterion='none')





