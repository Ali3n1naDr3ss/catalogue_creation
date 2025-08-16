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
import re
import sys
import numpy as np
from datetime import datetime
from hu_depth_codes import make_cutout, get_depths

######################### Set-up ####################################
baseDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/'

fields = ['COSMOS']
reqFilters = ['J', 'JH'] # <<<<<<<<<< may want to change this when running
testing = True           # <<<<<<<<<< may want to change this when running
size_arcsec = 32 #'full_size'      # <<<<<<<<<< may want to change this when running
saveFig = True

overwrite = False # overwrites: cutouts, depth/catalogues, depths/images, depths/results, cataloguing/plots
cutoutFile = os.path.join(baseDir, 'data/COSMOS/cutouts/cutoutNames.txt')

if testing:
    centre_ra = 149.5
    centre_dec = 2.2
    verbose = False

else:
    verbose = False
    overwrite = False # bc full-size result files are so large and time-consuming to produce, delete the file you don't want instead of writing over inadvertantly

###################### Apertures ##############################################
# required aperture diameters to run through
# change at your peril!
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])

########## Function call #######################################################
print("hu_find_depths.py run at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

imagePaths = []
cutouts = []

centre_ra_str = str(centre_ra).replace('.', '')
centre_dec_str = str(centre_dec).replace('.', '')

if size_arcsec == 'full_size':
    for filt in reqFilters:
        dataDir = os.path.join(baseDir, 'data/COSMOS')
        for file_ in os.listdir(dataDir):
                pattern = (
                    rf"{baseDir}data/COSMOS/UVISTA_{filt}_DR6\.fits")
                if re.fullmatch(pattern, file_):
                    imagePaths.append(file_)
                whtpattern = (
                    rf"{baseDir}data/COSMOS/cutouts/UVISTA_{filt}_DR6_wht\.fits")
                if re.fullmatch(whtpattern, file_):
                    imagePaths.append(file_)

if size_arcsec != 'full_size' and os.path.isfile(cutoutFile):
    with open(cutoutFile, 'r') as c:
        for line in c:
            line = line.strip()
            for filt in reqFilters:
                pattern = (
                    rf"{baseDir}data/COSMOS/cutouts/"
                    rf"UVISTA_{filt}_DR6_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                if re.fullmatch(pattern, line):
                    cutouts.append(line)
                whtpattern = (
                    rf"{baseDir}data/COSMOS/cutouts/"
                    rf"UVISTA_{filt}_DR6_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                if re.fullmatch(whtpattern, line):
                    cutouts.append(line)

for ff, fieldName in enumerate(fields):
    outputDir = '../depths/{0}'.format(fieldName) # 17/07/25

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    dataDir = baseDir + 'data/' + fieldName + '/'

    if size_arcsec != 'full_size' and len(cutouts) == 0: # if the correct cutoutsdo not exist, make them
        imagePaths = make_cutout(reqFilters, size_arcsec, centre_ra, centre_dec, dataDir, verbose=verbose, overwrite=overwrite) # cutout image Paths

        # get the depths
        get_depths(fieldName, imagePaths, size=str(size_arcsec), queue='none', reqFilters=reqFilters, overwrite=overwrite, outputDir='none', apDiametersAS=np.array([1.0, 1.8, 2.0, 3.0, 4.0]), ra_str=centre_ra_str, dec_str=centre_dec_str, verbose=verbose, saveFig=saveFig)

get_depths(fieldName, cutouts, size=str(size_arcsec), queue='none', reqFilters=reqFilters, overwrite=overwrite, outputDir='none', apDiametersAS=np.array([1.0, 1.8, 2.0, 3.0, 4.0]), ra_str=centre_ra_str, dec_str=centre_dec_str, verbose=verbose, saveFig=saveFig)




