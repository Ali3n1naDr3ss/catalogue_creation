#!/usr/bin/env python3

"""
stack_images.py

Create an image stack (e.g. Y+J), to use for a selection.

Created: Thurs 3rd March 2022

"""

# Import modules ###############################################

import numpy as np
import os
import astropy.io.fits as fits
import gc

# Setup ########################################################

filter1 = 'HSC-G_DR3' # Filters to be stacked
filter2 = 'HSC-R_DR3'
filter3 = 'HSC-I_DR3' # Set filter3 to 'None' if you only want to stack 2.

filter1 = 'UVISTA_H_DR6'
filter2 = 'UVISTA_K_DR6'
filter3 = 'None'

# Rename the filters when saving since the above UVISTA ones are a mess
name1 = 'H'
name2 = 'K'

filters = [filter1, filter2, filter3]

#fields = ['cdfs1', 'cdfs2', 'cdfs3'] # Lowercase for VIDEO file names.
#fields = ['xmm1', 'xmm2', 'xmm3']

fields = ['COSMOS']
#fields = ['XMM1', 'XMM2', 'XMM3']
#fields = ['CDFS1', 'CDFS2', 'CDFS3']

# Image dir for VIDEO
#imageDir = '/mnt/hoy/VIDEO_FINAL/'

# Image dir for UltraVISTA
imageDir = '/raid/scratch/data/COSMOS/'

stackDir = imageDir

print('Stacking {0} in {1}'.format(filters, fields))

# Loop through fields ##########################################

for i, fieldName in enumerate(fields):

    print('Running loop. Stacking {0} images in {1}'.format(filters, fieldName))

    # Set up images

    # VIDEO
    #image1 = imageDir + '{0}_{1}.fits'.format(fieldName, filter1)
    #image2 = imageDir + '{0}_{1}.fits'.format(fieldName, filter2)
    #if filter3 != 'None':
    #    image3 = imageDir + '{0}_{1}.fits'.format(fieldName, filter3)

    #wht1   = imageDir + '{0}_{1}_conf.fits'.format(fieldName, filter1)
    #wht2   = imageDir + '{0}_{1}_conf.fits'.format(fieldName, filter2)
    #if filter3 != 'None':
    #    wht3   = imageDir + '{0}_{1}_conf.fits'.format(fieldName, filter3)

    # HSC
    #imageDir = '/mnt/vardy/vardygroupshare/data/{0}/'.format(fieldName)

    image1 = imageDir + '{0}.fits'.format(filter1)
    image2 = imageDir + '{0}.fits'.format(filter2)
    if filter3 != 'None':
        image3 = imageDir + '{0}.fits'.format(filter3)

    wht1   = imageDir + '{0}_wht.fits'.format(filter1)
    wht2   = imageDir + '{0}_wht.fits'.format(filter2)
    if filter3 != 'None':
        wht3   = imageDir + '{0}_weight.fits'.format(filter3)

    # VOICE
    #imageDir = '/mnt/vardy/vardygroupshare/data/{0}/'.format(fieldName)

    #image1 = imageDir + '{0}_VOICE_{1}.fits'.format(fieldName, filter1)
    #image2 = imageDir + '{0}_VOICE_{1}.fits'.format(fieldName, filter2)
    #if filter3 != 'None':
    #    image3 = imageDir + '{0}_VOICE_{1}.fits'.format(fieldName, filter3)

    #wht1   = imageDir + '{0}_VOICE_{1}_wht.fits'.format(fieldName, filter1)
    #wht2   = imageDir + '{0}_VOICE_{1}_wht.fits'.format(fieldName, filter2)
    #if filter3 != 'None':
    #    wht3   = imageDir + '{0}_VOICE_{1}_wht.fits'.format(fieldName, filter3)

    # UltraVISTA
    #image1 = imageDir + 'UVISTA_{0}_dr5_rc1.fits'.format(filter1)
    #image2 = imageDir + 'UVISTA_{0}_dr5_rc1.fits'.format(filter2)

    #wht1   = imageDir + 'UVISTA_{0}_dr5_rc1_wht.fits'.format(filter1)
    #wht2   = imageDir + 'UVISTA_{0}_dr5_rc1_wht.fits'.format(filter2)

    # Existing stacks
    #imageDir = '/mnt/hoy/temporaryFilesROHAN/VIDEO/stacking/{0}/'.format(fieldName)

    #image1 = imageDir + '{0}_{1}_coadd.fits'.format(fieldName, filter1)
    #image2 = imageDir + '{0}_{1}_coadd.fits'.format(fieldName, filter2)

    #wht1 = imageDir + '{0}_{1}_coadd_wht.fits'.format(fieldName, filter1)
    #wht2 = imageDir + '{0}_{1}_coadd_wht.fits'.format(fieldName, filter2)


    print('First image: ' + image1)
    print('With weight: ' + wht1)
    print('Second image: ' + image2)
    print('With weight: ' + wht2)
    if filter3 != 'None':
        print('Third image: ' + image3)
        print('With weight: ' + wht3)
    # Load images
    imHDU1 = fits.open(image1, memmap=True) # Memmap loads parts of data in as necessary.
    imHDU2 = fits.open(image2, memmap=True)
    whtHDU1 = fits.open(wht1, memmap=True)
    whtHDU2 = fits.open(wht2, memmap=True)

    # Extract things
    image1 = imHDU1[0]
    header1 = image1.header
    image1 = image1.data
    imHDU1.close()

    print('##########Check header of image 1!#############')
    print(header1)

    wht1 = whtHDU1[0]
    wheader1 = wht1.header
    wht1 = wht1.data
    whtHDU1.close()

    image2 = imHDU2[0]
    header2 = image2.header
    image2 = image2.data
    imHDU2.close()

    wht2 = whtHDU2[0]
    wheader2 = wht2.header
    wht2 = wht2.data
    whtHDU2.close()

    # Extract third image if stacking
    if filter3 != 'None':
            imHDU3 = fits.open(image3, memmap=True)
            whtHDU3 = fits.open(wht3, memmap=True)

            image3 = imHDU3[0]
            header3 = image3.header
            image3 = image3.data
            imHDU3.close()

            wht3 = whtHDU3[0]
            wheader3 = wht3.header
            wht3 = wht3.data
            whtHDU3.close()

    # Stack = (Y / (sigma_Y)**2) + (J / (sigma_J)**2)
    # If VIDEO: the conf.fits files are already in inverse variance
    # So stack = w_Y*Y + w_J*J


    print('########## Computing inverse variance of the images ###########')
    wht1[wht1 == 0.] = 1.
    wht2[wht2 == 0.] = 1.

    wht1 = 1 / (wht1 ** 2)
    wht2 = 1 / (wht2 ** 2)


    print('########## Computing weighted sum...###########')

    w1_1 = wht1 * image1
    w2_2 = wht2 * image2

    # Manually free up some memory with garbage collector
    del image1
    del image2
    gc.collect()

    print('##### Weighted sum finished. Now computing weight sum...######')
    # Next, need w1 + w2 to normalise by.
    weightSum = wht1 + wht2 # This is also the stacked weight map, in terms of inverse variance.

    # Manually free up some memory with garbage collector
    del wht1
    del wht2
    gc.collect()

    # Compute stack!

    print('############ Stacking!!##############')
    finalImage = (w1_1 + w2_2) / weightSum

    # Convert nans to zero

    print('#################### Convert nans to zeroness...################')
    finalImage[np.isnan(finalImage)] = 0

    if filter3 != 'None':

        print('Now adding in third image')
        w3_3 = wht3 * image3

        # Manually free up some memory with garbage collector
        del image3
        gc.collect()

        # Next, need w1 + w2 to normalise by.
        weightSum = weightSum + wht3 # This is also the stacked weight map, in terms of inverse variance.

        # Manually free up some memory with garbage collector
        del wht3
        gc.collect()

        finalImage = (w1_1 + w2_2 + w3_3) / weightSum

        finalImage[np.isnan(finalImage)] = 0

    # Save image
    print('########Saving image.#########')
    
    outFile = stackDir + '/UVISTA_{1}{2}_DR6.fits'.format(fieldName.upper(), name1, name2)
    if filter3 == 'None':
        fits.writeto(stackDir + '/{0}_{1}{2}_DR6.fits'.format(fieldName.upper(), name1, name2), finalImage, header1, overwrite=False)
    if filter3 != 'None':
        fits.writeto(stackDir + '/{0}_{1}{2}{3}_drcoadd.fits'.format(fieldName.upper(), filter1, filter2, filter3), finalImage, header1, overwrite=True)

    print('###########Saved image to ' + stackDir + fieldName.upper())

    # Save weight
    print('########Saving weight#########')

    if filter3 == 'None':
        fits.writeto(stackDir + '/UVISTA_{1}{2}_DR6_wht.fits'.format(fieldName.upper(), name1, name2), weightSum, wheader1, overwrite=False)
    if filter3 != 'None':
        fits.writeto(stackDir + fieldName.upper() + '/{0}_{1}{2}{3}_coadd_wht.fits'.format(fieldName.upper(), filter1, filter2, filter3), weightSum, wheader1, overwrite=True)

    print('###########Saved weight to ' + stackDir + fieldName.upper())

