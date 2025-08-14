

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
def save_coordinates(directory, ra, dec, filename='random_coords.txt', verbose=True):
       
        # Ensure the directory exists
        os.makedirs(directory, exist_ok=True)
        # Full path to the file
        filepath = os.path.join(directory, filename)
        # Create the file if it doesn't exist
        if not os.path.exists(filepath):
            with open(filepath, 'w') as f:
                f.write("RA, Dec\n")  # header
        # Append the new line [RA, Dec]
        with open(filepath, 'a') as f:
            f.write(f"{ra}, {dec}\n")
        if verbose:
            print(f"Saved coordinates: RA = {ra}, Dec = {dec} to {filepath} \n")

def get_random_coords(file_):
    """ Output:
        random_coords(str): random line from the random_coords.txt in the format RA, Dec"""

    import random
    with open(file_, 'r') as f:
        lines = f.readlines()
    # Skip header if present
    if lines[0].strip().lower().startswith("ra"):
        lines = lines[1:]
    if not lines:
        raise ValueError("No coordinate lines found in the file.")
    random_coords = random.choice(lines).strip()
    print(f"Randomly selected RA/Dec: {random_coords}")
    print(type(random_coords))
    return random_coords

def get_fits_data(fits_img):
    with fits.open(fits_img) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
    return data, header, wcs

def make_random_coords(large_image, size_arcsec=1000, verbose=True):
    """ 
    Make random coords for de-bugging cutouts, 
    ensuring the random coords are within the bounds of the large image  """
    
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
        print("\n #################### RANDOM COORDS FOR TEST CUTOUT ####################")
        print("Bounds of image: ", large_image)
        print("RA bounds: ", min_ra, max_ra)
        print("Dec bounds:", min_dec, max_dec)
        print("Pixel scale: ", scale, "arcsec/pixel")
    ra = uniform(min_ra, max_ra)
    dec = uniform(min_dec, max_dec)
    if verbose: 
        print("Random Ra/Dec:", ra, dec, "\n")

    return ra, dec


def make_cutout(image_paths, size_arcsec=1000, ra='random', dec='random', corner='lower-left', output_path=None, verbose=True):
    """
    Makes square science and weight image cutouts, given the RA/Dec of a specified corner, and saves to data/{fieldName}/randomcutouts by default.
    
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
    if (ra and dec) == 'random':
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
            # make output path
            ra_filename = str(ra).replace('.', '')
            dec_filename = str(dec).replace('.', '')
            base = os.path.splitext(image)[0].replace('_backup', '')
            base_dir, instrument, fname_base = base.partition("UVISTA")

            # make dir if necessary
            randomCutoutDir = base_dir + 'random_cutouts/'
            if os.path.isdir(randomCutoutDir) == False:
                os.system('mkdir ' + randomCutoutDir)
            
            output_img_name = f'{randomCutoutDir}{instrument}{fname_base}_{ra_filename}_{dec_filename}.fits'
            save_coordinates(randomCutoutDir, ra_filename, dec_filename)
            hdul_out.writeto(output_img_name, overwrite=True)
            if verbose:
                print("\n #################### MAKING TEST CUTOUT ####################")
                print(f"Saved image cutout to: {output_img_name}\n")

            output_wht_name = f'{randomCutoutDir}{instrument}{fname_base}_wht_{ra_filename}_{dec_filename}.fits'
            hdul_out.writeto(output_wht_name, overwrite=True)
            if verbose:
                print(f"Saved weight cutout to: {output_wht_name} \n")

    return output_img_name, output_wht_name

def image_depth(imageName, zeropoint, back_size, back_filtersize, coords, apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtName = 'NONE', whtType = 'NONE', IRACapDiametersAS = np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segName = 'NONE', outputDir = 'none', filterName = 'NONE', numApertures = 300, step = 200, overwrite = False, inputSex = baseDir + 'data/bertin_config/video_mine.sex', strips = False, bgSub = True, mask = 'none', gridSepAS = 3.0):

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
            
        ### now run SE 
        coords = coords.replace(", ", "_") # get coords in correct format for filenaming

        inputSex = baseDir + 'data/bertin_config/video_mine.sex'

        segName = imageDir + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords + '_seg.fits'  
        bgSubName = imageDir  + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords + '_bgsub.fits'
        bgMapName = imageDir  + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords + '_bgmap.fits'
        outputCatalogue = catDir + 'd' + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords + '.fits'

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
    ra = coords.split('_')[0]
    dec = coords.split('_')[1]
    save_coordinates(catDir, ra, dec, verbose=True)
    return outputCatalogue

def read_image_lis(dirHere):
    """ Read .lis file"""

    # read filters
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print(f"Reading in images from {inputFile}")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
            
    return Table.read(inputFile, format = 'ascii.commented_header')


def get_depths(fieldName, back_size, back_filtersize, coords, queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir=baseDir+'data/', outputDir='none', overwrite=False):
    
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
    
    def get_cutout_img(imageName, coords):
        coords = coords.replace(", ", "_") # get coords in correct format for filenaming
        imageName = imageName.replace(".fits", "_.fits")
        imageName = imageName.replace(".fits", coords+".fits")
        return imageName

    ### Run image_depth for each filter ###
    # Run each filter as a separate file...
    for fi, filterName in enumerate(availableFilters):
        # define the images etc to send through
        imageName = imagedata['Image'][fi]
        imageName = get_cutout_img(imageName, coords)
        whtName = imagedata['Weight'][fi]
        whtName = get_cutout_img(whtName, coords)
        whtType = imagedata['Wht_type'][fi]
        zeropoint = imagedata['zeropoint'][fi]
        imageDir = imagedata['directory'][fi]
        maskName = '/raid/scratch/data/masks/{0}/{1}'.format(fieldName, imagedata['mask'][fi])
        
        if imageDir == 'here':
            imageDir = dataDir + fieldName + '/random_cutouts/'

        # set up for image_depth
        if queue == 'none':
            #print("No queue. \n")

            strips = True
            if 'COSMOS' in fieldName:
                jj = (filterName == stripFilt) # VISTA observation strategy produces images with 'strips' with alternating depth
                if np.any(jj):
                    strips = True
                    
            outputCatalogue =  image_depth(imageDir + imageName, zeropoint, back_size, back_filtersize, coords, whtName = imageDir + whtName, whtType = whtType, outputDir = outputDir, strips = strips, filterName = filterName, overwrite = overwrite, mask = maskName, gridSepAS = gridSepAS, apDiametersAS = apDiametersAS)
        
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
    subset_tab_name = filename_root + subset_tab_name
    #os.makedirs(filename_root, exist_ok=True) # create dir if it doesn't exist

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
    print(f"Saved subset to {outputPath}")
    bad_ratio = len(subset)/len(table)

    return outputPath, bad_ratio

def load_image(filepath):
    """Load FITS image and return the data array."""
    with fits.open(filepath) as hdul:
        data = hdul[0].data
    return data

def extract_filename(filename, verbose=False):
        """Extract substrings (e.g. 'JH') from a UVISTA cutout filename."""

        import re
        base = os.path.basename(filename)

        cutout_match = re.match(r"(UVISTA_)(\w+)(_DR6_)(\d+)_(\d+)\.fits", base) #\d+ matches digits
        wht_cutout_match = re.match(r"(UVISTA_)(\w+)(_DR6_wht_)(\d+)_(\d+)\.fits", base) #\d+ matches digits
        full_img_match = re.match(r"(UVISTA_)(\w+)(_DR6\.fits)", base) # \w+ matches word characters e.g. J, JH
        det_cat_match = re.match(r"d(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)", base)
        bad_cat_match = re.match(r"d(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)_bad_mags\.fits", base)
        bgmap_match = re.match(r"(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)_bgmap.fits", base)
        bgsub_match = re.match(r"(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)_bgsub.fits", base)
        bgseg_match = re.match(r"(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)_seg.fits", base)

        if cutout_match:
            prefix  = cutout_match.group(1)      # 'UVISTA_'
            filt    = cutout_match.group(2)        # 'JH'
            suffix  = cutout_match.group(3)      # '_DR6_'
            ra      = cutout_match.group(4)
            dec     = cutout_match.group(5)
            return prefix, suffix, filt, ra, dec

        if wht_cutout_match:
            prefix  = wht_cutout_match.group(1)      # 'UVISTA_'
            filt    = wht_cutout_match.group(2)        # 'JH'
            suffix  = wht_cutout_match.group(3)      # '_DR6_'
            ra      = wht_cutout_match.group(4)
            dec     = wht_cutout_match.group(5)
            return prefix, suffix, filt, ra, dec


        elif full_img_match:
            prefix  = full_img_match.group(1)    # 'UVISTA_'
            filt    = full_img_match.group(2)      # 'JH'
            suffix  = full_img_match.group(3)    # '_DR6.fits'
            return prefix, suffix, filt
        
        elif det_cat_match:
            filt    = det_cat_match.group(1)     # 'JH'
            bksz    = det_cat_match.group(2)     # '32'
            bkfilt  = det_cat_match.group(3)   # '9'
            ra      = det_cat_match.group(4)
            dec     = det_cat_match.group(5)
            return filt, bksz, bkfilt, ra, dec
        
        elif bad_cat_match:
            filt    = bad_cat_match.group(1)     # 'JH'
            bksz    = bad_cat_match.group(2)     # '32'
            bkfilt  = bad_cat_match.group(3)   # '9'
            ra      = bad_cat_match.group(4)
            dec     = bad_cat_match.group(5)
            return filt, bksz, bkfilt, ra, dec

        elif bgmap_match:
            filt    = bgmap_match.group(1)     # 'JH'
            bksz    = bgmap_match.group(2)     # '32'
            bkfilt  = bgmap_match.group(3)   # '9'
            ra      = bgmap_match.group(4)
            dec     = bgmap_match.group(5)
            return filt, bksz, bkfilt, ra, dec
        
        elif bgsub_match:
            filt    = bgsub_match.group(1)     # 'JH'
            bksz    = bgsub_match.group(2)     # '32'
            bkfilt  = bgsub_match.group(3)   # '9'
            ra      = bgsub_match.group(4)
            dec     = bgsub_match.group(5)
            return filt, bksz, bkfilt, ra, dec
        
        elif bgseg_match:
            filt    = bgseg_match.group(1)     # 'JH'
            bksz    = bgseg_match.group(2)     # '32'
            bkfilt  = bgseg_match.group(3)   # '9'
            ra      = bgseg_match.group(4)
            dec     = bgseg_match.group(5)
            return filt, bksz, bkfilt, ra, dec

        else:
            if verbose==True:
                print(f">>>>>>>>>>>>>>>>>> WARNING: Filename doesn't match expected pattern: \n {filename} \n {base} \n")
            return None


def find_last_matching_line(file1, file2):
    """
    Finds the last matching line between two text files.

    Args:
        file1 (str): Path to the first text file.
        file2 (str): Path to the second text file.

    Returns:
        str or None: The last matching line (stripped of whitespace), or None if no match is found.
    """
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = set(line.strip() for line in f1 if line.strip())
        lines2 = [line.strip() for line in f2 if line.strip()]

    # Iterate backwards through file2 lines to find the last match
    for line in reversed(lines2):
        if line in lines1:
            return line

    return None  # No matching lines found



def detections_fig(cutoutPath, detPath, filters, detCats, subsetCriterion='none'):
    """ 
    Makes a figure showing the science image, detections,
     and subset detections for two filters in a 3x2 formation. Subset creation happens here. 

        cutoutPath(str): /raid/scratch/data/COSMOS_test/random_cutouts/UVISTA_JH_DR6_15007186924767774_26432262759966996.fits"""
    
    sciImgs = []
    catPaths = [] # list of paths to cats which include all detections
    catTables = [] # data - all dets
    catNames = [] # subplot titles
    subsetPaths = [] # list of paths to cats which include a subset of dets  i.e. MAG_APER[1] > 50
    subsetTables = [] # data - only MAG_APER[1] > 50 dets
    subsetNames = [] # name of subset cat
    bad_ratios = [] # all detections : bad dets
    
    ### Get matching coords
    cutoutDir = cutoutPath.split(os.path.basename(cutoutPath))[0]
    cutout_coords_file = os.path.join(cutoutDir, 'random_coords.txt')
    detCat_coords_file = os.path.join(detPath, 'random_coords.txt')
    if os.path.isfile(cutout_coords_file):
        print(f"Found: {cutout_coords_file}")
    elif os.path.isfile(cutout_coords_file)==False:
        print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {cutout_coords_file}")
    if os.path.isfile(detCat_coords_file):
        print(f"Found: {detCat_coords_file}")
    elif os.path.isfile(detCat_coords_file)==False:
        print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {detCat_coords_file}")
    
    matching_coords = find_last_matching_line(cutout_coords_file, detCat_coords_file)
    if matching_coords:
        matching_coords = matching_coords.replace(', ', '_')
        print("Matching coords:", matching_coords)
    else:
        print("No matching coords found.")

    for f, filt in enumerate(filters):
        ## get science data
        # file with matching_coords in filename
        prefix, suffix, _, _, _ = extract_filename(cutoutPath) # get prefix and suffix to rebuild filename
        filename =  prefix + filt + suffix + matching_coords +'.fits' # rebuild filename with matched coords
        filePath = os.path.join(cutoutDir, filename)

        # ensure file exists
        if os.path.isfile(filePath):
            print(f"Found science image: {filePath}")
            sciImgs.append(filePath)
        elif os.path.isfile(filePath)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {filePath}")

        ## detections
        # get det with matching_coords in filename
        detName = detCats[f]
        _, bksz, bkfilt, ra, dec = extract_filename(detName)
        detName = 'd'+filt+'bksz'+bksz+'bkfilt'+bkfilt+'_'+matching_coords+'.fits'
        catPath = os.path.join(detPath, detName)
        # ensure file exists
        if os.path.isfile(catPath):
            print(f"Found det cat: {catPath}")
            catPaths.append(catPath)
            cat_table = Table.read(catPath) # a subset of detections
            catTables.append(cat_table)
            catNames.append(detName)
        elif os.path.isfile(catPath)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {catPath}")

    ## make subset data
        # file with matching_coords in filename
        subset_file, bad_ratio = make_subset(detCats[f], detPath, '_bad_mags.fits', subsetCriterion=subsetCriterion)
        bad_ratios.append(bad_ratio)
        print("bad rat bug",f, bad_ratios) # is this actually a bug or a limit on the smallest useful cutout (arcsec)
        # ensure file exists
        subsetPath = os.path.join(detPath, subset_file)
        if os.path.isfile(subsetPath):
            print(f"Found subset catalogue: {subsetPath}")
            subset_table = Table.read(subsetPath) # a subset of detections
            subsetTables.append(subset_table)
            subsetPaths.append(subsetPath)
            subset_name = os.path.splitext(os.path.basename(subsetPath))[0]
            subsetNames.append(subset_name)
        elif os.path.isfile(subsetPath)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {subsetPath}")
    
    ### initialise fig
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval() # z-scale all images
    from matplotlib.ticker import ScalarFormatter
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    for i, (sci_img, cat_table, subset_table, cat_name, subset_name) in enumerate(zip(sciImgs, catTables, subsetTables, catNames, subsetNames)):    
        sci_data = fits.getdata(sci_img) ## load sci data
        vmin, vmax = zscale.get_limits(sci_data)
        ax = axes[i, 0]
        ax.imshow(sci_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        plot_name = os.path.basename(sci_img)[:20]
        ax.set_title(f"{plot_name}", fontsize=10)
        ax.xaxis.set_major_formatter(formatter)
        ax.axis('off')

    ## plot detections
        # Extract RA and Dec
        ra = cat_table['ALPHA_J2000']
        dec = cat_table['DELTA_J2000']
        # plot
        ax_det = axes[i, 1]
        ax_det.scatter(ra, dec, s=1, color='blue', alpha=0.5)
        ax_det.set_xlabel("RA (deg)")
        ax_det.set_ylabel("Dec (deg)")
        filt, bksz, bkfilt, ra, dec = extract_filename(cat_name)
        ax_det.set_title(f"All detections: {filt} bksz:{bksz} bkfilt:{bkfilt}")
        ax_det.invert_xaxis()  # invert RA axis to match sky convention
        ax_det.xaxis.set_major_formatter(formatter)
        ax_det.grid(True)

    ## plot subset
        ra = subset_table['ALPHA_J2000']
        dec = subset_table['DELTA_J2000']
        ax_sub = axes[i, 2]
        ax_sub.scatter(ra, dec, s=1, color='blue', alpha=0.5)
        ax_sub.set_xlabel("RA (deg)")
        ax_sub.set_ylabel("Dec (deg)")
        filt, bksz, bkfilt, ra, dec = extract_filename(subset_name)
        ax_sub.set_title(f"subset: {filt} bksz:{bksz} bkfilt:{bkfilt}")
        ax_sub.xaxis.set_major_formatter(formatter)
        ax_sub.text(int(ra[:5])/100,int(dec[:5])/1e4, f'bad:all = {bad_ratios[i]:.5f}', fontsize=12, color='red')
        ax_sub.invert_xaxis()  # invert RA axis to match sky convention
        ax_sub.grid(True)

    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    #plt.show()

    return matching_coords

def bkg_plotter(param_combo, filters=['J', 'JH'], coords='full-size', depthDir='/raid/scratch/depths/COSMOS_test/images/', dataDir='/raid/scratch/data/COSMOS_test/random_cutouts/', verbose=True):
    """
    what affect do each set of SE bkg modelling parameters have?
    Do visual check with this side-by-side plotter

    param_combos(list): SE parameters used in bkg modelling test: [BACK_SIZE, BACK_FILTERSIZE]
    coords(list):       coordinates of random science cutout to match. [RA/Dec]
                        If default, mathces to full image (feature coming soon heheh) 
    depthDir(str):      directory of SE'd find_depth catalogues
    dataDir(str):       directory of science images (full-size and cutouts)
    """
    if coords == 'full-size':
        print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: i haven't done this part yet")
        exit()
    else:
        #get all images from dataDir/depthDir that have matching ra/dec, and matching param_combo
        dataFiles = []
        bgFiles = []
        matchSciImgs = []
        matchBgImgs = []

        all_data_files = [f for f in os.listdir(dataDir) if os.path.isfile(os.path.join(dataDir, f))]
        all_bg_imgs = [i for i in os.listdir(depthDir) if os.path.isfile(os.path.join(depthDir, i))] 

        for file_ in all_data_files:
            result = extract_filename(file_)
            if result is None:
                continue  # Skip non-matching files
            else:
                dataFiles.append(file_)

        for file_ in all_bg_imgs:
            result = extract_filename(file_)
            if result is None:
                continue  # Skip non-matching files 
            if re.match(r"\w+(bksz)\d+(bkfilt)\d+_\d+_\d+_bg(sub|map)\.fits", file_) or \
               re.match(r"\w+(bksz)\d+(bkfilt)\d+_\d+_\d+_seg\.fits", file_):
                #print(f"Found bkg image: {depthDir+file_}")
                bgFiles.append(depthDir+file_)

        # retireve data file with matching coords
        for file_ in dataFiles:
            prefix, suffix, filt, ra, dec = extract_filename(file_) # get bits to match to coords
            matchSciFname = 'UVISTA_' + filt + '_DR6_' + coords + '.fits'
            if file_ == matchSciFname:
                matchSciImgs.append(matchSciFname)
            else:
                if verbose:
                    print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: No matching coords for:", matchSciFname, '\n', file_)
       
        # retireve background file with matching coords and back_size, back_filtersize
        for file_ in bgFiles:
            filt, bksz, bkfilt, ra, dec = extract_filename(file_) # get bits to match to coords/param_combo
            matchBg = filt + 'bksz' +str(param_combo[0]) + 'bkfilt'+ str(param_combo[1]) + '_'+ coords
            if matchBg in file_:
                matchBgImgs.append(file_)
            else:
                if verbose:
                    print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: No matching params for:", matchBg, '\n', file_)


    ### initialise fig
    fig, axes = plt.subplots(2, 4, figsize=(12, 8))
    zscale = ZScaleInterval() # z-scale all images
    from matplotlib.ticker import ScalarFormatter
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    for i, filt in enumerate(filters):
        # Find the sci image for this filter
        sci_img = next((f for f in matchSciImgs if f"UVISTA_{filt}_" in f), None)
        if not sci_img:
            print(f">>>> Missing sci image for filter {filt}")
            continue
        ### load sci data
        sci_data = fits.getdata(dataDir+sci_img) 
        vmin, vmax = zscale.get_limits(sci_data)
        ax = axes[i, 0]
        ax.imshow(sci_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        plot_name = os.path.basename(sci_img)
        prefix, suffix, filt, ra, dec = extract_filename(plot_name)
        ax.set_title(f"Non-subtracted Image {filt}-band", fontsize=10)
        ax.xaxis.set_major_formatter(formatter)
        ax.axis('off')

        ### plot bg 
    for i, filt in enumerate(filters):
            # Find the bg image for this filter
            matchBgImgsFilt = [f for f in matchBgImgs if f"{filt}bksz" in f]
            if not matchBgImgsFilt:
                print(f">>>> Missing bg image for filter {filt}")
                continue

            breakpoint()
            for bg_img in matchBgImgsFilt:
                if bg_img.endswith('_bgsub.fits'):
                    bgsub_data = fits.getdata(bg_img) 
                    vmin, vmax = zscale.get_limits(bgsub_data)
                    ax_bgsub = axes[i, 1]
                    ax_bgsub.imshow(bgsub_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                    plot_name = os.path.basename(bg_img)
                    filt, bksz, bkfilt, ra, dec = extract_filename(plot_name)
                    ax_bgsub.set_title(f"Bg-sub'ed {filt} BKSZ:{bksz} BKFILT:{bkfilt}", fontsize=10)
                    ax_bgsub.xaxis.set_major_formatter(formatter)
                    ax_bgsub.axis('off')

                elif bg_img.endswith('_bgmap.fits'):
                    bgmap_data = fits.getdata(bg_img) 
                    vmin, vmax = zscale.get_limits(bgmap_data)
                    ax_bgmap = axes[i, 2]
                    ax_bgmap.imshow(bgmap_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                    plot_name = os.path.basename(bg_img)
                    filt, bksz, bkfilt, ra, dec = extract_filename(plot_name)
                    ax_bgmap.set_title(f"Bg-Map {filt} BKSZ:{bksz} BKFILT:{bkfilt}", fontsize=10)
                    ax_bgmap.xaxis.set_major_formatter(formatter)
                    ax_bgmap.axis('off')

                elif bg_img.endswith('_seg.fits'):
                    seg_data = fits.getdata(bg_img) 
                    vmin, vmax = zscale.get_limits(seg_data)
                    ax_seg = axes[i, 3]
                    ax_seg.imshow(seg_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                    plot_name = os.path.basename(bg_img)
                    filt, bksz, bkfilt, ra, dec = extract_filename(plot_name)
                    ax_seg.set_title(f"Seg Map {filt} BKSZ:{bksz} BKFILT:{bkfilt}", fontsize=10)
                    ax_seg.xaxis.set_major_formatter(formatter)
                    ax_seg.axis('off')
    
    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    plt.show()













def detections_fig_full(img_filename, filters, detCats, subsetCriterion='none'):
    """ 
    Makes a figure showing the science image, detections,
     and subset detections for two filters in a 3x2 formation """
    ### Top row - filter:X science image, detections, subset-detections
    ### Find data
    
    sciImgs = []
    catNames = [] # all detections in a given list of cats
    catTables = []
    subsetDirs = [] # detections with bad magnitudes i.e. MAG_APER[1] > 50
    subsetTables = []
    subsetNames = []
    bad_ratios = []

    for f, filt in enumerate(filters):
        ## Retrieve science image data
        img_filename = os.path.basename(img_filename)
        prefix, suffix, _ = extract_filename(img_filename)
        filename = f'{prefix}{filt}{suffix}{ra}{dec}'
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
        ax_sub.text(np.max(ra), np.min(dec)-0.05, f'bad:all = {bad_ratios[i]:.5f}', fontsize=12, color='red')
        ax_sub.invert_xaxis()  # invert RA axis to match sky convention
        ax_sub.grid(True)

    plt.tight_layout()
    plt.show()



