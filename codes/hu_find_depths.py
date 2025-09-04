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
# 'Y', 'J', 'H', 'K', 'JH', 'HK'
# 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y'
reqFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']  # <<<<<<<<<< may want to change 
size_arcsec = 100 #'full_size'         # <<<<<<<<<< may want to change 
saveFig = True            # <<<<<<<<<< may want to change 

centre_ra = 149.5         # <<<<<<<<<< may want to change 
centre_dec = 2.2          # <<<<<<<<<< may want to change        

verbose = False
overwrite = False # bc full-size result files are so large and time-consuming to produce, delete the file you don't want instead of writing over inadvertantly
# overwrite = True overwrites: cutouts, depth/catalogues, depths/images, depths/results
# NOTE: cataloguing/plots are always re-written unless hu_depth_codes: detections_fig and bg_plotter are updated

###################### Apertures ##############################################
# required aperture diameters to run through
# change at your peril!
apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])
apDiametersAS = np.array([1.8])
########## Function call #######################################################
print("hu_find_depths.py run at:", datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

dataDir = os.path.join(baseDir, 'data/{0}/'.format('COSMOS'))
cutoutDir = os.path.join(dataDir, f'cutouts/')
cutoutFile = os.path.join(cutoutDir, 'cutoutNames.txt')

imagePaths = []
cutouts = []
cutoutsToMatch = []

centre_ra_str = str(centre_ra).replace('.', '')
centre_dec_str = str(centre_dec).replace('.', '')

if size_arcsec == 'full_size':
    for filt in reqFilters:
        for file_ in os.listdir(dataDir):
                pattern = f"UVISTA_{filt}_DR6.fits"
                HSCpattern = f"{filt}_DR3.fits"
                if file_ == pattern:
                    imagePath = os.path.join(dataDir, file_)
                    imagePaths.append(imagePath)
                elif file_ == HSCpattern:
                    imagePath = os.path.join(dataDir, file_)
                    imagePaths.append(imagePath)
                whtpattern = f"UVISTA_{filt}_DR6_wht.fits"
                #if file_ == whtpattern:
                    #imagePath = os.path.join(dataDir, file_)
                    #imagePaths.append(imagePath)
            
if size_arcsec != 'full_size' and os.path.isfile(cutoutFile):
    with open(cutoutFile, 'r') as c:
        for line in c:
            line = line.strip()
            for filt in reqFilters:
                if 'HSC_' in filt:
                    cutToMatch = f"{cutoutDir}{filt}_DR3_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
                    whtToMatch = f"{cutoutDir}{filt}_DR3_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
                    if cutToMatch not in cutoutsToMatch:
                        cutoutsToMatch.append(cutToMatch)
                    if whtToMatch not in cutoutsToMatch:
                        cutoutsToMatch.append(whtToMatch)

                    pattern = (
                        rf"{cutoutDir}"
                        rf"{filt}_DR3_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                    if re.fullmatch(pattern, line):
                        cutouts.append(line)
                    whtpattern = (
                        rf"{cutoutDir}"
                        rf"{filt}_DR3_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                    if re.fullmatch(whtpattern, line):
                        cutouts.append(line)

                else:
                    cutToMatch = f"{cutoutDir}UVISTA_{filt}_DR6_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
                    whtToMatch = f"{cutoutDir}UVISTA_{filt}_DR6_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
                    if cutToMatch not in cutoutsToMatch:
                        cutoutsToMatch.append(cutToMatch)
                    if whtToMatch not in cutoutsToMatch:
                        cutoutsToMatch.append(whtToMatch)

                    pattern = (
                        rf"{cutoutDir}"
                        rf"UVISTA_{filt}_DR6_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                    if re.fullmatch(pattern, line):
                        cutouts.append(line)
                    whtpattern = (
                        rf"{cutoutDir}"
                        rf"UVISTA_{filt}_DR6_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                    if re.fullmatch(whtpattern, line):
                        cutouts.append(line)

elif os.path.isfile(cutoutFile) != True:
    print(">>>>>>>>>>>>>>>>>> WARNING: Could not find {cutoutFile}...")

for ff, fieldName in enumerate(fields):
    outputDir = os.path.join(baseDir,'depths/{0}/'.format(fieldName)) # 17/07/25

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)

    if size_arcsec != 'full_size':
        if cutouts != cutoutsToMatch: # if the correct cutouts do not exist, make them
            unmatchedCutouts = make_cutout(reqFilters, size_arcsec, centre_ra, centre_dec, dataDir, verbose=verbose, overwrite=overwrite) # cutout Paths
            cutouts = cutouts + unmatchedCutouts

        # get the depths
        print("doing small size")

        get_depths(fieldName, cutouts=cutouts, size=str(size_arcsec), queue='none', reqFilters=reqFilters, overwrite=overwrite, outputDir=dataDir, dataDir=dataDir, apDiametersAS=np.array([1.0, 1.8, 2.0, 3.0, 4.0]), ra_str=centre_ra_str, dec_str=centre_dec_str, verbose=verbose, saveFig=saveFig)

    elif size_arcsec == 'full_size':
        print("doing full size")
        get_depths(fieldName, fullsizeimages=imagePaths, size=str(size_arcsec), queue='none', reqFilters=reqFilters, dataDir=dataDir, overwrite=overwrite, outputDir='none', apDiametersAS=np.array([1.0, 1.8, 2.0, 3.0, 4.0]), ra_str=centre_ra_str, dec_str=centre_dec_str, verbose=verbose, saveFig=saveFig)




