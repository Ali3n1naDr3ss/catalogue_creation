#!/usr/bin/env python3

"""
stack_images.py

Create an image stack (e.g. Y+J), to use for a selection.

Created:  Thurs 3rd March 2022
Modified: Thurs 21st August 2025

"""

# Import modules ###############################################

import numpy as np
import os
import astropy.io.fits as fits
import gc

# Setup ########################################################

# Filters to be stacked
filter1 = 'UVISTA_J_DR6'
filter2 = 'UVISTA_H_DR6'
filter3 = 'None'    # Set filter3 to 'None' if you only want to stack 2.

# Rename the filters when saving since the above UVISTA ones are a mess
name1 = 'H'
name2 = 'K'

filters = [filter1, filter2, filter3]
fields = ['COSMOS']

# Image dir for UltraVISTA
imageDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/COSMOS/'
#stackDir = os.path.join(imageDir, 'test')
stackDir =imageDir

print('Stacking {0} in {1}'.format(filters, fields))

def load_fits(path):
    with fits.open(path, memmap=True) as hdul:
        # Memmap loads parts of data in as necessary.
        hdu = hdul[0]
        return hdu.header, hdu.data

# Loop through fields ##########################################

for i, fieldName in enumerate(fields):
    print('Running loop. Stacking {0} images in {1}'.format(filters, fieldName))

    # Set up images
    imagePath1 = imageDir + '{0}.fits'.format(filter1)
    imagePath2 = imageDir + '{0}.fits'.format(filter2)
    if filter3 != 'None':
        imagePath3 = imageDir + '{0}.fits'.format(filter3)

    whtPath1   = imageDir + '{0}_wht.fits'.format(filter1)
    whtPath2   = imageDir + '{0}_wht.fits'.format(filter2)
    if filter3 != 'None':
        whtPath3   = imageDir + '{0}_wht.fits'.format(filter3)

    print('First image: ' + imagePath1)
    print('With weight: ' + whtPath1)
    print('Second image: ' + imagePath2)
    print('With weight: ' + whtPath2)
    if filter3 != 'None':
        print('Third image: ' + imagePath3)
        print('With weight: ' + whtPath3)

    # Load images
    header1, data1 = load_fits(imagePath1)
    whtHeader1, whtData1 = load_fits(whtPath1)
    header2, data2 = load_fits(imagePath2)
    whtHeader2, whtData2 = load_fits(whtPath2)

    # Extract third image if stacking
    if filter3 != 'None':
                header3, data3 = load_fits(imagePath3)
                whtHeader3, whtData3 = load_fits(whtPath3)

    # Stack = (Y / (sigma_Y)**2) + (J / (sigma_J)**2)
    # If VIDEO or UVISTA: the conf.fits files are already in inverse variance
    # So stack = w_Y*Y + w_J*J

    print('########## Computing inverse variance of the images ###########')
    whtData1[whtData1 == 0.] = np.nan # mask out invalied weights? or set to placeholder? #previously set to 1.
    whtData2[whtData2 == 0.] = np.nan #previously set to 1.

    #whtData1 = 1 / (whtData1 ** 2) # inverse square of weight
    #whtData2 = 1 / (whtData2 ** 2)

    print('########## Computing weighted sum...###########')

    w1_1 = whtData1 * data1 # data * inverse variance of weight
    w2_2 = whtData2 * data2 # i.e. weight by inverse variance 

    # Manually free up some memory with garbage collector
    del data1
    del data2
    gc.collect()

    print('##### Weighted sum finished. Now computing weight sum...######')
    # Next, need w1 + w2 to normalise by.
    weightSum = whtData1 + whtData2 # This is also the stacked weight map, in terms of inverse variance.

    # Manually free up some memory with garbage collector
    del whtData1
    del whtData2
    gc.collect()

    # Compute stack!

    print('############ Stacking!!##############')
    finalImage = (w1_1 + w2_2) / weightSum

    # Convert nans to zero

    print('#################### Convert nans to zeroness...################')
    finalImage[np.isnan(finalImage)] = 0

    if filter3 != 'None':

        print('Now adding in third image')
        w3_3 = whtData3 * data3

        # Manually free up some memory with garbage collector
        del data3
        gc.collect()

        # Next, need w1 + w2 to normalise by.
        weightSum = weightSum + whtData3 # This is also the stacked weight map, in terms of inverse variance.

        # Manually free up some memory with garbage collector
        del whtData3
        gc.collect()

        finalImage = (w1_1 + w2_2 + w3_3) / weightSum

        finalImage[np.isnan(finalImage)] = 0

    # Save image
    print('########Saving image.#########')
    
    if filter3 == 'None':
        fits.writeto(stackDir + '/{0}_{1}{2}_DR6.fits'.format(fieldName.upper(), name1, name2), finalImage, header1, overwrite=False)
    if filter3 != 'None':
        fits.writeto(stackDir + '/{0}_{1}{2}{3}_drcoadd.fits'.format(fieldName.upper(), filter1, filter2, filter3), finalImage, header1, overwrite=True)

    print('###########Saved image to ' + stackDir + fieldName.upper())

    # Save weight
    print('########Saving weight#########')

    if filter3 == 'None':
        fits.writeto(stackDir + '/UVISTA_{1}{2}_DR6_wht.fits'.format(fieldName.upper(), name1, name2), weightSum, whtHeader1, overwrite=False)
    if filter3 != 'None':
        fits.writeto(stackDir + fieldName.upper() + '/{0}_{1}{2}{3}_coadd_wht.fits'.format(fieldName.upper(), filter1, filter2, filter3), weightSum, whtHeader1, overwrite=True)

    print('###########Saved weight to ' + stackDir + fieldName.upper())

