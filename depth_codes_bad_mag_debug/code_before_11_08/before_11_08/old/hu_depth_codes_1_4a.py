

"""
17/07/25
hu_depths_codes.py
Testing new_depths_codes.py - is SE running correctly? Why aren't all the objects recorded with a measured magnitude?

"""
import os
import re
import sys
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from collections import defaultdict
from astropy.visualization import ZScaleInterval

baseDir = '/raid/scratch/'
dataDir = baseDir + 'data/COSMOS_test/'
outputCatalogues = []

##############################################################################################################


def get_fits_data(fits_img):
    with fits.open(fits_img) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
    return data, header, wcs

def make_random_coords(large_image, size_arcsec=1000, verbose=True):
    """ 
    Make random coords for de-bugging cutouts, 
    ensuring the random coords are actually in the large image  """
    
    from random import uniform

    data, header, wcs = get_fits_data(large_image)
    n_dec, n_ra = data.shape # shape = (rows, cols)

    # get pix scale
    try:
        scale = abs(header['CD1_1']) * 3600  # deg to arcsec
    except KeyError:
        scale = abs(header['CDELT1']) * 3600
    size_pix = int(size_arcsec / scale)

    ## get limits from original image
    corners_x = [0, n_ra - 1] # list of pixel values from zero to the right-most pixel (n_ra - 1)
    corners_y = [0, n_dec - 1]
    # get the 4 corners of the iage
    xx, yy = np.meshgrid(corners_x, corners_y) # create all combinations of x and y values
    ra_corners, dec_corners = wcs.all_pix2world(xx, yy, 0) # convert to RA/Dec
    
    # create RA/Dec bounds
    min_ra = np.min(ra_corners)
    max_ra = np.max(ra_corners)
    min_dec = np.min(dec_corners)
    max_dec = np.max(dec_corners)
    if verbose:
        print("Bounds of image: ", large_image)
        print("RA bounds: ", min_ra, max_ra)
        print("Dec bounds:", min_dec, max_dec)
        print("Pixel scale: ", scale, "arcsec/pixel")
    ra = uniform(min_ra, max_ra)
    dec = uniform(min_dec, max_dec)
    if verbose: 
        print("Random Ra/Dec:", ra, dec)

    return ra, dec


def make_cutout(image_paths, size_arcsec=1000, ra='random', dec='random', corner='lower-left', output_path=None, verbose=True):
    """
    Create a square cut-out from a large image, given the RA/Dec of a specified corner.
    
    image_path(str):    Full path to the original large FITS image.
    ra, dec(float):     RA and Dec (in degrees) of one corner of the square.
    size_arcsec(float): Size of the square cutout (arcsec on a side).
    corner(str):        Which corner the (ra, dec) point is (options: 'lower-left', 'upper-left', 'lower-right', 'upper-right').
    output_path(str):   Path to save the output cutout. If None, will auto-generate.
    """
    import astropy.units as u
    from astropy.nddata import Cutout2D
    from astropy.coordinates import SkyCoord
    
    data, header, wcs = get_fits_data(image_paths[0])
    # get coords of cutout
    if (ra and dec) == 'random':  ## get random positions 
        ra, dec = make_random_coords(image_paths[0], verbose=verbose)

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
    for image in image_paths:
        if output_path is None:
            ra_filename = str(ra).replace('.', '')
            dec_filename = str(dec).replace('.', '')
            base = os.path.splitext(image)[0].replace('_backup', '')
            output_img_name = f'{base}_{ra_filename}_{dec_filename}.fits'
            hdul_out.writeto(output_img_name, overwrite=True)
            print(f"\n Saved image cutout to: {output_img_name}")
            output_wht_name = f'{base}_wht_{ra_filename}_{dec_filename}.fits'
            hdul_out.writeto(output_wht_name, overwrite=True)
            print(f"\n Saved weight cutout to: {output_wht_name}")
    return output_img_name, output_wht_name

def image_depth(imageName, zeropoint, back_size = back_size, back_filtersize = back_filtersize, apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtName = 'NONE', whtType = 'NONE', IRACapDiametersAS = np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segName = 'NONE', outputDir = 'none', filterName = 'NONE', numApertures = 300, step = 200, overwrite = False, inputSex = baseDir + 'data/bertin_config/video_mine.sex', strips = False, bgSub = True, mask = 'none', gridSepAS = 3.0):

    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    
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
    #print("The base name is ", baseName)

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
    #print("The pixel scale is {0:.4f}".format(pixScale))
    
    if segName == 'NONE':
        
        # run source extractor
        ## define aperture sizes
        apStringPix = str(apDiametersAS[0]/pixScale) 
        for ai in range(1, apDiametersAS.size):
            apStringPix = apStringPix + ',' + str(apDiametersAS[ai]/pixScale)
            
        # now run SE

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


def get_depths(fieldName, back_size=32, back_filtersize=9, queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir=baseDir+'data/', outputDir='none', overwrite=False):
    
    
    # set the grid sep
    gridSepAS = 3.0
        
    # Read in the images.lis file
    dirHere = dataDir + fieldName + '/'
    imagedata = read_image_lis(dirHere)
    availableFilters = np.array(imagedata['Name'])
    #print("The available filters are:\n", availableFilters)
    print("Running get_depths", dirHere)
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
                    
            outputCatalogue =  image_depth(imageDir + imageName, back_size=back_size, back_filtersize=back_filtersize, zeropoint, whtName = imageDir + whtName, whtType = whtType, outputDir = outputDir, strips = strips, filterName = filterName, overwrite = overwrite, mask = maskName, gridSepAS = gridSepAS, apDiametersAS = apDiametersAS)
        
            outputCatalogues.append(outputCatalogue)

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
    bad_ratio = len(subset)/len(table)

    return outputPath, bad_ratio

def load_image(filepath):
    """Load FITS image and return the data array."""
    with fits.open(filepath) as hdul:
        data = hdul[0].data
    return data

def extract_filename(filename):
        """Extract the filter substring (e.g. 'JH') from a UVISTA cutout filename."""
        import re

        base = os.path.basename(filename)
        cutout_match = re.match(r"(UVISTA_)\w+(_DR6_\d+_\d+\.fits)", base) #\d+ matches digits
        full_img_match = re.match(r"(UVISTA_)\w+(_DR6\.fits)", base) # \w+ matches word characters e.g. J, JH

        if cutout_match:
            prefix = cutout_match.group(1)
            suffix = cutout_match.group(2)
            return prefix, suffix

        elif full_img_match:
            prefix = full_img_match.group(1)
            suffix = full_img_match.group(2)
            return prefix, suffix

        else:
            raise ValueError(f"Filename doesn't match expected pattern: \n {filename} \n {base} \n")

def detections_fig(cutout_filename, filters, detCats, subsetCriterion='none'):
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
    bad_ratios = []

    for f, filt in enumerate(filters):
        ## Retrieve science image data
        cutout_filename = os.path.basename(cutout_filename)
        prefix, suffix = extract_filename(cutout_filename)
        filename = f'{prefix}{filt}{suffix}'
        sciDir = dataDir + filename # testing on vista only
        
        ## science image
        if os.path.isfile(sciDir):
            print(f"Found science image: {sciDir}")
            sciImgs.append(sciDir)
        elif os.path.isfile(sciDir)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {sciDir}")

        ## detections
        if os.path.isfile(detCats[f]):
            print(f"Found detection catalogue: {detCats[f]}")
            det_table = Table.read(detCats[f]) # all detections in a given catalogue
            detTables.append(det_table)
            det_name = os.path.splitext(os.path.basename(detCats[f]))[0]
            detNames.append(det_name)
        elif os.path.isfile(detCats[f])==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {detCats[f]}")
        
        ## get subset data
        subset_dir, bad_ratio = make_subset(detCats[f], baseDir+'depths/COSMOS_test/catalogues/', '_bad_mags', subsetCriterion=subsetCriterion)
        bad_ratios.append(bad_ratio)
        if os.path.isfile(subset_dir):
            print(f"Found detection catalogue: {subset_dir}")
            subset_table = Table.read(subset_dir) # a subset of detections
            subsetTables.append(subset_table)
            subsetDirs.append(subset_dir)
            subset_name = os.path.splitext(os.path.basename(subset_dir))[0]
            subsetNames.append(subset_name)
        elif os.path.isfile(subset_dir)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {subset_dir}")

    ### initialise fig
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval() # z-scale all images

    for i, (sci_img, det_table, det_name, subset_table, subset_name) in enumerate(zip(sciImgs, detTables, detNames, subsetTables, subsetNames)):
        sci_data = fits.getdata(sci_img) ## load sci data
        vmin, vmax = zscale.get_limits(sci_data)
        ax = axes[i, 0]
        ax.imshow(sci_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        plot_name = os.path.basename(sci_img)[:20]
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
        ax_sub.text(np.max(ra), np.min(dec)-0.002, f'{bad_ratios[i]}', fontsize=12, color='red')
        ax_sub.invert_xaxis()  # invert RA axis to match sky convention
        ax_sub.grid(True)

    plt.tight_layout()
    plt.show()


def detections_fig_full(img_filename, filters, detCats, subsetCriterion='none'):
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
    bad_ratios = []
    print(f'>>>>>>>>>>>>>>>>>>>>> bad rats {bad_ratios}')
    for f, filt in enumerate(filters):
        ## Retrieve science image data
        img_filename = os.path.basename(img_filename)
        prefix, suffix = extract_filename(img_filename)
        filename = f'{prefix}{filt}{suffix}'
        sciDir = dataDir + filename # testing on vista only
        
        ## science image
        if os.path.isfile(sciDir):
            print(f"Found science image: {sciDir}")
            sciImgs.append(sciDir)
        elif os.path.isfile(sciDir)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {sciDir}")

        ## detections
        if os.path.isfile(detCats[f]):
            print(f"Found detection catalogue: {detCats[f]}")
            det_table = Table.read(detCats[f]) # all detections in a given catalogue
            detTables.append(det_table)
            det_name = os.path.splitext(os.path.basename(detCats[f]))[0]
            detNames.append(det_name)
        elif os.path.isfile(detCats[f])==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {detCats[f]}")
        
        ## get subset data
        subset_dir, bad_ratio = make_subset(detCats[f], baseDir+'depths/COSMOS_test/catalogues/', '_bad_mags', subsetCriterion=subsetCriterion)
        bad_ratios.append(bad_ratio)
        if os.path.isfile(subset_dir):
            print(f"Found detection catalogue: {subset_dir}")
            subset_table = Table.read(subset_dir) # a subset of detections
            subsetTables.append(subset_table)
            subsetDirs.append(subset_dir)
            subset_name = os.path.splitext(os.path.basename(subset_dir))[0]
            subsetNames.append(subset_name)
        elif os.path.isfile(subset_dir)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {subset_dir}")

    ### initialise fig
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval() # z-scale all images

    for i, (sci_img, det_table, det_name, subset_table, subset_name) in enumerate(zip(sciImgs, detTables, detNames, subsetTables, subsetNames)):
        sci_data = fits.getdata(sci_img) ## load sci data
        vmin, vmax = zscale.get_limits(sci_data)
        ax = axes[i, 0]
        ax.imshow(sci_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        plot_name = os.path.basename(sci_img)[:20]
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
        ax_sub.text(np.max(ra), np.min(dec)-0.002, f'{bad_ratios[i]}', fontsize=12, color='red')
        ax_sub.invert_xaxis()  # invert RA axis to match sky convention
        ax_sub.grid(True)

    plt.tight_layout()
    plt.show()


def bkg_plotter(param_combos, depthDir='/raid/scratch/depths/COSMOS_test/', dataDir='/raid/scratch/data/COSMOS_test', verbose=True):
    """
    what affect do each set of SE bkg modelling parameters have?
    Do visual check with this side-by-side plotter

    param_combos(list): SE parameters used in bkg modelling test: [BACK_SIZE, BACK_FILTERSIZE]
    depthDir(str):      directory of SE'd find_depth catalogues
    dataDir(str):       directory of science images (full-size and cutouts)
    """

    ImgDir = depthDir + 'images/' # bkg map and bkg-sub'd img

    ## get all filenames in ImgDir and verify given test parameters are present
    if os.listdir(ImgDir)==False:
        print(f">>>>>>>>>>>>>>>>>>>>>>>> WARNING: ImageDir NOT found: {testImgDir}")
    else:
        all_files = [f for f in os.listdir(ImgDir) if os.path.isfile(os.path.join(ImgDir, f))]
        comparison_imgs = []    
        for f in all_files:
            for param_combo in param_combos:
                if str(param_combo) in f:
                    comparison_imgs.append(ImgDir+f)
                    if verbose:
                        print(f"Bkg model comparison image found: {f}")
                else:
                     print(f">>>>>>>>>>>>>>>>>>>>>>>> WARNING: Bkg model comparison image NOT found: {f}")

    ## Specify which bkg params you want to compare  
    pattern = re.compile(r"^(?P<band>.+?)bk_sz(?P<bksz>.+?)bkfisz(?P<bkfsz>.+?)_(?P<type>.+?)\.fits")

    # Group files by (bksz, bkfsz)
    grouped = {}
    bands = []
    for fname in comparison_imgs:
        match = pattern.match(os.path.basename(fname)) # extracts info from filename if name matches supplied pattern
        if not match:
            continue  # Skip files that don't match the pattern
            if verbose:
                print(f"Skipping unmatched filename: {fname}")

        band = match.group("band") # filter/band
        bands.append(band)
        bksz = match.group("bksz")
        bkfsz = match.group("bkfsz")
        ftype = match.group("type")

        key = (bksz, bkfsz) # combination of test params

        if key not in grouped:
            grouped[key] = {"J": {}, "JH": {}} # add keys for subdicts 
            # TODO: make flexible in terms of filter choice
            print(f"Adding: grouped[{key}][{band}][{ftype}]")
        grouped[key][band][ftype] = os.path.join(ImgDir, fname)

    # get all filenames in science image dir and verify test bands are present
    all_files = [f for f in os.listdir(dataDir) if os.path.isfile(os.path.join(dataDir, f))]
    sci_imgs = []   
    orig_pattern = re.compile(r"^UVISTA_(?P<band>.+?)\_DR6.fits") # TODO: make flexible

    for f in all_files:
        basename = os.path.basename(f)
        match = orig_pattern.match(basename)
        if match:
            matched_band = match.group("band")
            if verbose:
                print(f"Science image found: {f}")
            if matched_band in ['J', 'JH']:
                sci_imgs.append(os.path.join(dataDir, f))
                if verbose:
                    print(f"Visualising background model for: {sci_imgs}")

    # display bgmap, bgsub, original for both filters

    # Find group with matching files for J and JH
    selected_group = None
    for (bksz, bkfsz), bands in grouped.items():
        if all(ftype in bands["J"] for ftype in ["bgmap", "bgsub"]) and \
           all(ftype in bands["JH"] for ftype in ["bgmap", "bgsub"]):
            selected_group = (bksz, bkfsz)
            break

    if not selected_group:
        raise ValueError("No group found with complete J and JH bgmap, and bgsub images.")

        bksz, bkfsz = selected_group

    ### Load bgmap and bgsub images TODO: hard-coded to J and JH only
    j_images = []
    jh_images = []
    ftypes = ["bgmap", "bgsub"]

    for ftype in ftypes:
        j_path = grouped[(bksz, bkfsz)]["J"][ftype]   # TODO: not flexible
        jh_path = grouped[(bksz, bkfsz)]["JH"][ftype] # TODO: not flexible

        j_data = fits.getdata(j_path)
        jh_data = fits.getdata(jh_path)

        j_images.append((j_data, ftype))
        jh_images.append((jh_data, ftype))

    ### Load science images
    sciImg_pattern = re.compile(r"^UVISTA_(?P<band>.+?)_DR6\.fits")

    for img_path in sci_imgs:
        basename = os.path.basename(img_path)
        match = orig_pattern.match(basename)
        if not match:
            continue

        band = match.group("band")
        data = fits.getdata(img_path)

        if band == "J":
            j_images.append((data, "original"))
        elif band == "JH":
            jh_images.append((data, "original"))

    # TODO: add check that science image corresponds to the bkg e.g. cutout coords
    
    ### Plot
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval()

    ## Top row: J-band
    for i, (data, ftype) in enumerate(j_images):
        vmin, vmax = zscale.get_limits(data)
        ax = axes[0, i]
        ax.imshow(data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"J bk_sz{bksz} bkfisz{bkfsz} {ftype}", fontsize=12)
        ax.axis('off')

    # Bottom row: JH-band
    for i, (data, ftype) in enumerate(jh_images):
        vmin, vmax = zscale.get_limits(data)
        ax = axes[1, i]
        ax.imshow(data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"JH bk_sz{bksz} bkfisz{bkfsz} {ftype}", fontsize=12, color='red')
        ax.axis('off')

    plt.tight_layout()
    plt.show()


