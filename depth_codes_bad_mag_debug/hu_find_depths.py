""" hu_find_depths.py
11/08/2025

Why do the stacked images have such a high proportion of "bad" magnitude measurements?
bad: MAG_APER[1]>50
bad proportion for JH: 21%
bad proportion for J: 1%

Are the segmentation maps being produced correctly?

"""

###################### Import useful modules #########################
import os
import sys
import numpy as np
from datetime import datetime
from hu_depth_codes import make_cutout, get_depths

######################### Set-up ####################################
baseDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/'

fields = ['COSMOS']
reqFilters = ['J', 'JH'] # <<<<<<<<<< may want to change this when running
testing = False          # <<<<<<<<<< may want to change this when running
size_arcsec = 100        # <<<<<<<<<< may want to change this when running

if testing:
    centre_ra = 149.5
    centre_dec = 2.2
    verbose = False
    overwrite = True # overwrites: cutouts, deptha/catalogues, depths/images, depths/results
else:
    cutouts = os.path.join(baseDir, 'data/COSMOS/cutouts/cutoutNames.txt')
    verbose = False
    overwrite = False # bc full-size result files are so large and time-consuming to produce, delete the file you don't want instead of writing over inadvertantly


###################### Apertures ##############################################
# required aperture diameters to run through
# change at your peril!
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])

########## Function call #######################################################

print("hu_find_depths.py run at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

for ff, fieldName in enumerate(fields):
    outputDir = '../depths/{0}'.format(fieldName) # HOLLY 17/07/25 doing some tests

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    dataDir = baseDir + 'data/' + fieldName + '/'

    if testing:
        cutouts = make_cutout(reqFilters, size_arcsec, centre_ra, centre_dec, dataDir, verbose=verbose, overwrite=overwrite) # todo: write filenames of cutouts to file so i dont have to generate them each time

    # get the depths
    # read in cutouts file as a list
    with open(cutouts, 'r') as f:
        cutouts = [line.strip() for line in f if line.strip()]  # strip removes \n and skips empty lines
    
    get_depths(fieldName, cutouts, size=str(size_arcsec), queue='none', reqFilters=reqFilters, overwrite=overwrite, outputDir='none', apDiametersAS=np.array([1.0, 1.8, 2.0, 3.0, 4.0]))


