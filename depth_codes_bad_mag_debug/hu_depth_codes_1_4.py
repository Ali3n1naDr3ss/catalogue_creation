

"""
17/07/25
hu_depths_codes.py
Testing new_depths_codes.py - is SE running correctly? Why aren't all the objects recorded with a measured magnitude?

"""
import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table


from astropy.visualization import ZScaleInterval

baseDir = '/raid/scratch/'
dataDir = baseDir + 'data/COSMOS/'
outputCatalogues = []
##############################################################################################################

# create random 1000^2 pixel cutouts of original image and run pipeline again to aid de-bugging

#!/usr/bin/env python3

def make_cutout(image_path, filt, size_arcsec, ra=None, dec=None, corner='lower-left', output_path=None):
    """
    Create a square cut-out from a large image, given the RA/Dec of a specified corner.
    
    image_path(str):    Full path to the original large FITS image.
    ra, dec(float):     RA and Dec (in degrees) of one corner of the square.
    size_arcsec(float): Size of the square cutout (arcsec on a side).
    corner(str):        Which corner the (ra, dec) point is (options: 'lower-left', 'upper-left', 'lower-right', 'upper-right').
    output_path(str):   Path to save the output cutout. If None, will auto-generate.
    """
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    #

    with fits.open(image_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

        # Convert RA/Dec corner to pixel coordinates
        skycoord = SkyCoord(ra, dec, unit='deg')
        x, y = wcs.world_to_pixel(skycoord)

        # Get scale in arcsec/pixel
        try:
            scale = abs(header['CD1_1']) * 3600  # deg to arcsec
        except KeyError:
            scale = abs(header['CDELT1']) * 3600

        size_pix = int(size_arcsec / scale)

        # Shift the position to the center based on corner
        if corner == 'lower-left':
            x_center = x + size_pix / 2
            y_center = y + size_pix / 2
        elif corner == 'upper-left':
            x_center = x + size_pix / 2
            y_center = y - size_pix / 2
        elif corner == 'lower-right':
            x_center = x - size_pix / 2
            y_center = y + size_pix / 2
        elif corner == 'upper-right':
            x_center = x - size_pix / 2
            y_center = y - size_pix / 2
        else:
            raise ValueError("Invalid corner type. Choose from: 'lower-left', 'upper-left', 'lower-right', 'upper-right'.")

        # Perform the cutout
        position = (x_center, y_center)
        cutout = Cutout2D(data, position, (size_pix, size_pix), wcs=wcs)

        # Create new FITS HDU
        hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
        hdul_out = fits.HDUList([hdu])

        # Save
        if output_path is None:
            output_path = os.path.splitext(image_path)[0].replace('_backup', f'{ra}_{dec}') +'.fits'
        hdul_out.writeto(output_path, overwrite=True)
        print(f"\n Saved cutout to {output_path}")

"""filter_ = filt
image_path = '/raid/scratch/data/COSMOS/UVISTA_'+filter_
image_type = '_DR6_backup.fits'
#image_type = '_DR6_wht_backup.fits'
#image_path = '/raid/scratch/depths/COSMOS/phot/'+filter_+'_aperPhot.fits'
cutout_size = size_arcsec    # arcseconds
if ra or dec not None:
    ra_corner = ra   # in degrees
    dec_corner = dec    # in degrees"""


def image_depth(imageName, zeropoint, apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtName = 'NONE', whtType = 'NONE', IRACapDiametersAS = np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segName = 'NONE', outputDir = 'none', filterName = 'NONE', numApertures = 300, step = 200, overwrite = False, inputSex = baseDir + 'data/bertin_config/video_mine.sex', strips = False, bgSub = True, mask = 'none', gridSepAS = 3.0):

    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    
    from astropy.io import fits
    # Tinsley SE env
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/share/sextractor'

    # define the output files
    if outputDir == 'none':
        outputDir = 'depths/'
        if os.path.isdir(outputDir) == False:
            os.system('mkdir ' + outputDir)

    # seg map, bkg map, bgsub images
    imageDir = outputDir + 'images/'
    if os.path.isdir(imageDir) == False:
        os.system('mkdir ' + imageDir)

    plotDir = outputDir + 'plots/'
    if os.path.isdir(plotDir) == False:
        os.system('mkdir ' + plotDir)

    catDir = outputDir + 'catalogues/'
    if os.path.isdir(catDir) == False:
        os.system('mkdir ' + catDir)

    aperDir = outputDir + 'phot/'
    if os.path.isdir(aperDir) == False:
        os.system('mkdir ' + aperDir)

    resultsDir = outputDir + 'results/'
    if os.path.isdir(resultsDir) == False:
        os.system('mkdir ' + resultsDir)
        
    parts = imageName.split('/') # split image dir by /
    if len(parts) > 1:
        baseName = parts[-1] # image filename string

    # also remove the file extension
    pparts = baseName.split('.')
    baseName = pparts[0] # image name string without file extension
    print("The base name is ", baseName)

    # check all the necessary files exist
    imyes = os.path.isfile(imageName)
    if imyes == False:
        print("ERROR: the image file does not exist, check the path:\n" + imageName)
        exit()

    if whtType != "NONE":
        whtyes = os.path.isfile(whtName) # find correct wight file for image
        if whtyes == False:
            print("ERROR: the weight file does not exist, check the path:\n" + whtName)            
            exit()
            
    ######### seg map, SE, and Bkg Subtraction ###########################################

    if filterName == 'NONE':
        filterName = baseName

    ### get the pixel scale ###
    hdulist = fits.open(imageName)
    imageHeader = hdulist[0].header

    if 'CD1_1' in imageHeader:
        cdone_o = -3600.0*imageHeader['CD1_1']
    else:
        cdone_o = 3600.0*np.abs(imageHeader['CDELT1'])
    pixScale = round(cdone_o, 5)
    print("The pixel scale is {0:.4f}".format(pixScale))
    
    if segName == 'NONE':
        
        # run source extractor
        ## define aperture sizes
        apStringPix = str(apDiametersAS[0]/pixScale) 
        for ai in range(1, apDiametersAS.size):
            apStringPix = apStringPix + ',' + str(apDiametersAS[ai]/pixScale)
            
        # now run SE
        #HOLLY 29/05/25 - improve bkg sub for VISTA
        back_size = 32
        back_filtersize = 9
        inputSex = baseDir + 'data/bertin_config/video_mine.sex'

        segName = imageDir + filterName + 'bk_sz' + str(back_size) + 'bkfisz' + str(back_filtersize) + '_seg.fits'  
        #bgSubName = imageDir  + filterName + '_bgsub.fits'
        bgSubName = imageDir  + filterName + 'bk_sz' + str(back_size) + 'bkfisz' + str(back_filtersize) + '_bgsub.fits'
        bgMapName = imageDir  + filterName + 'bk_sz' + str(back_size) + 'bkfisz' + str(back_filtersize) + '_bgmap.fits'
        outputCatalogue = catDir + 'd' + filterName + 'bk_sz' + str(back_size) + 'bkfisz' + str(back_filtersize) + '.fits'

        keywordsbase = ' -CATALOG_TYPE FITS_1.0 -CATALOG_NAME ' + \
                       outputCatalogue + \
                       ' -MAG_ZEROPOINT '+ str(zeropoint) + \
                       ' -WEIGHT_TYPE ' + whtType + \
                       ' -WEIGHT_IMAGE ' + whtName
        # -CHECKIMAGE_TYPE "BACKGROUND" for background map -BACKGROUND for bkg subtracted
        keywords = keywordsbase + \
                   ' -CHECKIMAGE_TYPE "BACKGROUND,-BACKGROUND,APERTURES" ' + \
                   '-CHECKIMAGE_NAME ' + bgMapName + ',' + bgSubName + ',' + segName + ' ' + \
                   '-PHOT_APERTURES ' + apStringPix + \
                   ' -BACK_SIZE ' + str(back_size) + \
                   ' -BACK_FILTERSIZE ' +str(back_filtersize)
    
        command = '/usr/local/sextractor/bin/sex '+ imageName +' -c ' + inputSex + keywords
        
        if os.path.isfile(bgSubName) == False or os.path.isfile(segName) == False or overwrite:
            print("The SEG and BG subtracted map do not exist.  Running SE like so: \n")
            print(command)
            os.system(command)
        else:
            print(f"The SEG and/or BG subtracted map exist at: {segName, bgSubName} \n")

    # next step after plotting is to place apertures down (see RB version) 17/07/25 ###########################################################
    # plotting
    # Fig 1. science img, detections, bad mags for J, JH
    # Fig 2. science img, bg sub, bag map for J, JH
    #bkg_fig()

    return outputCatalogue

def read_image_lis(dirHere):
    """ Read .lis file"""

    # read filters
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
            
    return Table.read(inputFile, format = 'ascii.commented_header')


def get_depths(fieldName, queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir=baseDir+'data/', outputDir='none', overwrite=False):

    
    # set the grid sep
    gridSepAS = 3.0
        
    # Read in the images.lis file
    dirHere = dataDir + fieldName + '/'
    imagedata = read_image_lis(dirHere)
    availableFilters = np.array(imagedata['Name'])
    #print("The available filters are:\n", availableFilters)

    if reqFilters[0] != 'all':

        # reduce the list of available filters
        keep = np.zeros(availableFilters.size, dtype=bool)
        for bf, keepfilt in enumerate(reqFilters):
            ii = np.where((availableFilters == keepfilt))
            keep[ii] = True
            
        imagedata = imagedata[keep]
        availableFilters = imagedata['Name']
        
    #print("Only running for filters: ",'\n', availableFilters)

    # check if we are in UltraVISTA (
    stripFilt = np.array(['Y', 'J', 'H', 'Ks', 'JH', 'YJ', 'HKs', 'Y_DR4', 'J_DR4', 'H_DR4', 'Ks_DR4', 'NB118', 'NB118_DR4'])

    ### Run image_depth for each filter ###
    # Run each filter as a separate file...
    for fi, filterName in enumerate(availableFilters):
        # define the images etc to send through
        imageName = imagedata['Image'][fi]
        whtName = imagedata['Weight'][fi]
        whtType = imagedata['Wht_type'][fi]
        zeropoint = imagedata['zeropoint'][fi]
        imageDir = imagedata['directory'][fi]
        maskName = '/raid/scratch/data/masks/{0}/{1}'.format(fieldName, imagedata['mask'][fi])
        
        if imageDir == 'here':
            imageDir = dataDir + fieldName + '/'

        # set up for image_depth
        if queue == 'none':
            #print("No queue. \n")

            strips = True
            if 'COSMOS' in fieldName:
                jj = (filterName == stripFilt) # VISTA observation strategy produces images with 'strips' with alternating depth
                if np.any(jj):
                    strips = True
                    
            outputCatalogue =  image_depth(imageDir + imageName, zeropoint, whtName = imageDir + whtName, whtType = whtType, outputDir = outputDir, strips = strips, filterName = filterName, overwrite = overwrite, mask = maskName, gridSepAS = gridSepAS, apDiametersAS = apDiametersAS)
        
            outputCatalogues.append(outputCatalogue)
            print(outputCatalogues)
        else:
            print("queue =! 'none': ", queue)
            strips = "False"
            if fieldName == 'COSMOS':
                jj = (filterName == stripFilt)
                if np.any(jj):
                    strips = "True"

            strips = "True"
            
            # make an ap diameters string
            apDiametersAS = np.array(apDiametersAS)
            apDiametersASstring = '{0:.2f}'.format(apDiametersAS[0])
            for i in range(apDiametersAS.size-1):

                apDiametersASstring = apDiametersASstring + ',{0:.2f}'.format(apDiametersAS[i+1])
            
            #print(apDiametersASstring)
            #print("Spawning in the queue...", queue)
            # make shell script
            tmpName = "tmp_{1}_{0}.sh".format(filterName, fieldName)
            f = open(tmpName, 'w')
            f.write('#!/bin/bash\n')
            f.write('python3 stupid.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(imageDir + imageName, imageDir +  whtName, whtType, zeropoint, outputDir, strips, filterName, overwrite, maskName, gridSepAS, apDiametersASstring))
            f.close()
            
            # now execute this
            command = "addqueue -c 'tmp_{0}' -m 9 -q {0} -d ./{1}".format(queue, tmpName)
            #print(command)
            os.system('chmod +x {0}'.format(tmpName))
            os.system(command)
    
    return outputCatalogues

def plot_detections(depthCat, save=True, overwrite=True):
    """
    Produces a sky-plot of all the detections that SE makes.
    """
    # detections
    plot_name = os.path.splitext(os.path.basename(depthCat))[0]
    
    # Load table
    table = Table.read(depthCat)

    # Extract RA and Dec
    ra = table['ALPHA_J2000']
    dec = table['DELTA_J2000']

    # Create scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(ra, dec, s=1, color='blue', alpha=0.5)
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.title(f"All detections: {plot_name}")
    plt.gca().invert_xaxis()  # Optional: invert RA axis to match sky convention
    plt.grid(True)
    if save==True and overwrite==True:
        saveDir = '/raid/scratch/depths/COSMOS_test/catalogues/' + 'testFigs/'
        os.makedirs(saveDir, exist_ok=True)
        plt.savefig(saveDir+plot_name, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {saveDir}{plot_name}")
    #plt.show()
    # TODO: show J and JH side by side


def make_subset(big_cat, outputDir, subset_tab_name, subsetCriterion='none'):
    """
     Make a subset from a larger astropy table
     accoring to a single criterion e.g. ["MAG_APER"][:, 1] > 50 
     subset_name(str) substring to modify filename for table subset """

    ### set up
    filename_root = os.path.splitext(os.path.basename(big_cat))[0]
    subset_tab_name = filename_root + subset_tab_name + '.fits'
    os.makedirs(filename_root, exist_ok=True) # create dir if it doesn't exist

    table = Table.read(big_cat)
    print(f"Objects in unfiltered catalogue: {len(table)}")
    
    if subsetCriterion is not'none':
        subset = table[table[subsetCriterion]]  # criterion by which to create a subset
    else:
        print("No subset criterion supplied. Using default: ['MAG_APER'][:, 1] > 50 ")
        subset = table[table["MAG_APER"][:, 1] > 50 ]  # criterion by which to create a subset
    outputPath = os.path.join(outputDir, subset_tab_name)  
    subset.write(outputPath, overwrite=True)
    print(f"Objects in filtered sub-catalogue: {len(subset)}")
    print(f"Saved to {outputPath}")
    return outputPath

def load_image(filepath):
    """Load FITS image and return the data array."""
    with fits.open(filepath) as hdul:
        data = hdul[0].data
    return data

def detections_fig(filters, detCats, subsetCriterion='none'):
    """ 
    Makes a figure showing the science image, detections,
     and subset detections for two filters in a 3x2 formation """
    ### Top row - filter:X science image, detections, subset-detections
    ### Find data
    
    sciImgs = []
    detNames = [] # all detections in a given list of cats
    detTables = []
    subsetDirs = [] # detections with bad magnitudes i.e. MAG_APER[1] > 50
    subsetTables = []
    subsetNames = []

    for f, filt in enumerate(filters):
        filename = f'UVISTA_{filt}_DR6.fits'
        sciDir = dataDir + filename # testing on vista only
        
        ## science image
        if os.path.isfile(sciDir):
            print(f"Found science image: {sciDir}")
            sciImgs.append(sciDir)
        elif os.path.isfile(sciDir)==False:
            print(f"Missing file: {sciDir}")

        ## detections
        if os.path.isfile(detCats[f]):
            print(f"Found detection catalogue: {detCats[f]}")
            det_table = Table.read(detCats[f]) # all detections in a given catalogue
            detTables.append(det_table)
            det_name = os.path.splitext(os.path.basename(detCats[f]))[0]
            detNames.append(det_name)
        elif os.path.isfile(detCats[f])==False:
            print(f"Missing file: {detCats[f]}")
        
        ## get subset data
        subset_dir = make_subset(detCats[f], baseDir+'depths/COSMOS_test/catalogues/', '_bad_mags', subsetCriterion=subsetCriterion)
        if os.path.isfile(subset_dir):
            print(f"Found detection catalogue: {subset_dir}")
            subset_table = Table.read(subset_dir) # a subset of detections
            subsetTables.append(subset_table)
            subsetDirs.append(subset_dir)
            subset_name = os.path.splitext(os.path.basename(subset_dir))[0]
            subsetNames.append(subset_name)
        elif os.path.isfile(subset_dir)==False:
            print(f"Missing file: {subset_dir}")

    ### initialise fig
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval() # z-scale all images
    for i, (sci_img, det_table, det_name, subset_table, subset_name) in enumerate(zip(sciImgs, detTables, detNames, subsetTables, subsetNames)):
        sci_data = fits.getdata(sci_img) ## load sci data
        vmin, vmax = zscale.get_limits(sci_data)
        ax = axes[i, 0]
        ax.imshow(sci_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        plot_name = os.path.basename(sci_img)
        ax.set_title(f"{plot_name}", fontsize=10)
        ax.axis('off')

    ## plot detections
        # Extract RA and Dec
        ra = det_table['ALPHA_J2000']
        dec = det_table['DELTA_J2000']
        # plot
        ax_det = axes[i, 1]
        ax_det.scatter(ra, dec, s=1, color='blue', alpha=0.5)
        ax_det.set_xlabel("RA (deg)")
        ax_det.set_ylabel("Dec (deg)")
        ax_det.set_title(f"All detections: {det_name}")
        ax_det.invert_xaxis()  # invert RA axis to match sky convention
        ax_det.grid(True)

    ## plot subset
        ra = subset_table['ALPHA_J2000']
        dec = subset_table['DELTA_J2000']
        ax_sub = axes[i, 2]
        ax_sub.scatter(ra, dec, s=1, color='blue', alpha=0.5)
        ax_sub.set_xlabel("RA (deg)")
        ax_sub.set_ylabel("Dec (deg)")
        ax_sub.set_title(f"subset of detections: {subset_name}")
        ax_sub.invert_xaxis()  # invert RA axis to match sky convention
        ax_sub.grid(True)

    plt.tight_layout()
    plt.show()



