

"""
17/07/25
hu_depths_codes.py
Testing new_depths_codes.py - is SE running correctly? Why aren't all the objects recorded with a measured magnitude?

"""
import os
import re
import sys
import numpy as np
import warnings
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from collections import defaultdict
from astropy.visualization import ZScaleInterval

baseDir = '/raid/scratch/'
dataDir = baseDir + 'data/COSMOS_test/'
outputCatalogues = []

############################ MAKE RANDOM CUTOUTS #######################################################

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
    ra = round(ra, 5)
    dec = round(dec, 5)
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
        ra, dec = make_random_coords(image_paths[0], size_arcsec=size_arcsec, verbose=verbose)

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

    # Copy the original header and update it with the cutout WCS
    new_header = header.copy()
    new_header.update(cutout.wcs.to_header())

    hdu = fits.PrimaryHDU(data=cutout.data, header=new_header)
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
            print(f"Writing random coordinates file...")
            hdul_out.writeto(output_img_name, overwrite=True)
            if verbose:
                print("\n #################### MAKING TEST CUTOUT ####################")
                print(f"Saved image cutout to: {output_img_name}\n")

            output_wht_name = f'{randomCutoutDir}{instrument}{fname_base}_wht_{ra_filename}_{dec_filename}.fits'
            hdul_out.writeto(output_wht_name, overwrite=True)
            if verbose:
                print(f"Saved weight cutout to: {output_wht_name} \n")

    return output_img_name, output_wht_name

############################ RETREIEVE RANDOM CUTOUTS #######################################################

def extract_filename(filename, verbose=False):
        """Extract substrings (e.g. 'JH') from a UVISTA cutout filename."""

        import re
        base = os.path.basename(filename)

        cutout_match = re.match(r"(UVISTA_)(\w+)(_DR6_)(\d+)_(\d+)\.fits", base) #\d+ matches digits
        wht_cutout_match = re.match(r"(UVISTA_)(\w+)(_DR6_wht_)(\d+)_(\d+)\.fits", base) 
        full_img_match = re.match(r"(UVISTA_)(\w+)(_DR6\.fits)", base) # \w+ matches word characters e.g. J, JH
        full_wht_match = re.match(r"(UVISTA_)(\w+)(_DR6_wht\.fits)", base) 
        orig_img_match = re.match(r"(UVISTA_)(\w+)(_DR6_backup\.fits)", base)
        fullsize_cat_match = re.match(r"d(\w+)bksz(\d+)bkfilt(\d+)_full\-size\.fits", base)

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

        elif wht_cutout_match:
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

        elif orig_img_match:
            prefix  = orig_img_match.group(1)    # 'UVISTA_'
            filt    = orig_img_match.group(2)      # 'JH'
            suffix  = orig_img_match.group(3)    # '_DR6.fits'
            return prefix, suffix, filt

        elif full_wht_match :
            prefix  = full_wht_match.group(1)    # 'UVISTA_'
            filt    = full_wht_match.group(2)      # 'JH'
            suffix  = full_wht_match.group(3)    # '_DR6_backup.fits'
            return prefix, suffix, filt
        
        elif fullsize_cat_match:
            filt    = fullsize_cat_match.group(1)     # 'JH'
            bksz    = fullsize_cat_match.group(2)     # '32'
            bkfilt  = fullsize_cat_match.group(3)     # '9'
            return filt, bksz, bkfilt

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


def find_matching_coord_files(coords, param_combo, dataDir, depthDir, verbose=True):

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
    if verbose:
        print("dataFiles: ", dataFiles)

    for file_ in all_bg_imgs:
        result = extract_filename(file_)
        if result is None:
            continue  # Skip non-matching files 
        if re.match(r"\w+(bksz)\d+(bkfilt)\d+_\d+_\d+_bg(sub|map)\.fits", file_) or \
           re.match(r"\w+(bksz)\d+(bkfilt)\d+_\d+_\d+_seg\.fits", file_):
            #print(f"Found bkg image: {depthDir+file_}")
            bgFiles.append(depthDir+file_)
        if verbose:
            print("bgFiles: ", bgFiles)

    # retireve data file with matching coords
    for file_ in dataFiles:
        result  = extract_filename(file_) 
        if len(result) >3:
            prefix, suffix, filt, ra, dec = extract_filename(file_) # get bits to match to coords
        else:
            prefix, filt, suffix = extract_filename(file_)
        matchSciFname = 'UVISTA_' + filt + '_DR6_' + coords + '.fits'
        if file_ == matchSciFname:
            matchSciImgs.append(matchSciFname)
        #else:
            #if verbose:
                #print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: Check coord format; No matching coords for:", matchSciFname, '\n', file_)
   
    # retireve background file with matching coords and back_size, back_filtersize
    for file_ in bgFiles:
        filt, bksz, bkfilt, ra, dec = extract_filename(file_) # get bits to match to coords/param_combo
        matchBg = filt + 'bksz' +str(param_combo[0]) + 'bkfilt'+ str(param_combo[1]) + '_'+ coords
        if matchBg in file_:
            matchBgImgs.append(file_)
        #else:
            #if verbose:
                #print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: No matching params for:", matchBg, '\n', file_)

    return matchSciImgs, matchBgImgs 


def find_fullsize_files(param_combo, dataDir, depthDir, verbose=True):

    #get all images from dataDir/depthDir that have matching param_combo
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
    if verbose:
        print("dataFiles: ", dataFiles)

    for file_ in all_bg_imgs:
        result = extract_filename(file_)
        if result is None:
            continue  # Skip non-matching files 
        if re.match(r"\w+(bksz)\d+(bkfilt)\d+_full\-size_bg(sub|map)\.fits", file_) or \
           re.match(r"\w+(bksz)\d+(bkfilt)\d+_full\-size_seg\.fits", file_):
            #print(f"Found bkg image: {depthDir+file_}")
            bgFiles.append(depthDir+file_)
        if verbose:
            print("bgFiles: ", bgFiles)

    # retireve data file with matching coords
    for file_ in dataFiles:
        result = extract_filename(file_)
        if len(result) > 3:
            continue

        prefix, filt, suffix = extract_filename(file_) # get bits to match to coords
        matchSciFname = 'UVISTA_' + filt + '_DR6_backup.fits'
        if file_ == matchSciFname:
            matchSciImgs.append(matchSciFname)
        #else:
            #if verbose:
                #print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: Check coord format; No matching coords for:", matchSciFname, '\n', file_)
   
    # retireve background file with matching coords and back_size, back_filtersize
    for file_ in bgFiles:
        filt, bksz, bkfilt = extract_filename(file_) # get bits to match to coords/param_combo
        matchBg = filt + 'bksz' +str(param_combo[0]) + 'bkfilt'+ str(param_combo[1]) + '_full-size_'
        if matchBg in file_:
            matchBgImgs.append(file_)
        #else:
            #if verbose:
                #print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: No matching params for:", matchBg, '\n', file_)

    return matchSciImgs, matchBgImgs 




############################ SOURCE EXTRACTOR #######################################################

def image_depth(imageName, zeropoint, back_size, back_filtersize, coords, apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtName = 'NONE', whtType = 'NONE', IRACapDiametersAS = np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segName = 'NONE', outputDir = 'none', filterName = 'NONE', numApertures = 300, step = 200, overwrite = False, inputSex = baseDir + 'data/bertin_config/video_mine.sex', strips = False, bgSub = True, mask = 'none', gridSepAS = 3.0):

    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    
    # Tinsley SE env
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/share/sextractor'

    # define the output files
    if outputDir == 'none':
        outputDir = '../depths/'
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

    # get the pixel scale
    hdulist = fits.open(imageName)
    imageHeader = hdulist[0].header
    imageHeaderKeys = list(imageHeader.keys())

    # ensure header of data file is as expected 
    # i.e. if using a copy/cutout, header should == header of original file

    imgBase = os.path.basename(imageName)
    if 'DR6_' in imgBase:
        originalImageName = imgBase.split('DR6_')[0] + 'DR6_backup.fits'
        originalHDUlist = fits.open('/raid/scratch/data/COSMOS/'+originalImageName)
        originalImageHeaderKeys = list(originalHDUlist[0].header.keys())

        if imageHeaderKeys == originalImageHeaderKeys:
            print("ImageHeader matches expected format.")
        elif imageHeaderKeys != originalImageHeaderKeys:

            imageHeaderKeys = set(imageHeaderKeys)
            originalImageHeaderKeys = set(originalImageHeaderKeys)
            uniqueImageHeaderKeys = sorted(imageHeaderKeys - originalImageHeaderKeys)
            uniqueOrigImageHeaderKeys = sorted(originalImageHeaderKeys - imageHeaderKeys)

            warnings.warn(">>>>>>> WARNING: imageHeader does not match expected Header. Known issues: incorrect calculation of pixScale. \n")
            print("Keys only in imageHeader: ",uniqueImageHeaderKeys)
            print("Keys only in originalImageHeader: ",uniqueOrigImageHeaderKeys)
    else:
        cd1Test = ['CD1_1' in imageHeaderKeys]
        print(f">>>> WARNING: Passing file to SE. Header includes CD1_1?: {cd1Test} Is this the right file?", '\n', imageName, '\n') 
        cdone_o = -3600.0*imageHeader['CD1_1'] # TODO - this might break - not tested yet

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
            
        ### now run SE 
        coords = coords.replace(", ", "_") # get coords in correct format for filenaming

        inputSex = baseDir + 'data/bertin_config/video_mine.sex'
        param_file = baseDir + 'data/bertin_config/default.param'

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
                   ' -CHECKIMAGE_TYPE "BACKGROUND,-BACKGROUND,SEGMENTATION" ' + \
                   '-CHECKIMAGE_NAME ' + bgMapName + ',' + bgSubName + ',' + segName + ' ' + \
                   '-PHOT_APERTURES ' + apStringPix + \
                   ' -BACK_SIZE ' + str(back_size) + \
                   ' -BACK_FILTERSIZE ' +str(back_filtersize)
    
        command = '/usr/local/sextractor/bin/sex '+ imageName +' -c ' + inputSex + ' -PARAMETERS_NAME ' + param_file + keywords
        
        if os.path.isfile(bgSubName) == False or os.path.isfile(segName) == False or overwrite:
            print("The SEG and BG subtracted map do not exist, or overwrite:True.  Running SE like so: \n")

            print(command)
            os.system(command)
        else:
            print(f"The SEG and/or BG subtracted map exist at: {segName, bgSubName} \n")
 
    #######################################################################
    # Next step is to place apertures down
    aperPhotFile = aperDir + filterName + '_aperPhot.fits'
    overwrite = False
    if os.path.isfile(aperPhotFile) == False or overwrite == True:
        # define the grid separation
        gridSepPixels = np.ceil(gridSepAS/pixScale) # 5'' separation
        # gridSepPixels = 10.0
        
        # make this tunable...
        if bgSub == False:
            bgSubName = imageName
            print("Not using bg subtracted image.")
        print("Measuring the aperture photometry.")

        ii = segName.rfind('NIRSPEC')
        if ii > 0:
            field = 'NIRSPEC'
        else:
            field = 'NONE'
        # if aperphotfile doesn't exist, call ap phot blank - the output fits is aperphotfile
        aperture_photometry_blank(bgSubName, segName, whtName, apDiametersAS, gridSeparation = gridSepPixels, clean = True, outputFitsName = aperPhotFile, imageDir = imageDir, field = field, overwrite = overwrite)
        print("aperture_photometry_blank is running")

    
    #######################################################################
    # Then calculate the local depths, and make a nice plot
    # if COSMOS, I need to run in strips too
    recalculate = False

    # mask
    #print("aperPhotFile", aperPhotFile, bgSubName)
    #bgSubName = 'JHbk_sz32bkfisz9_bgmap.fits'
    regions, globaldepths, meddepths, modedepths = extract_local_depths(aperPhotFile, apDiametersAS, zeropoint, recalculate = recalculate, numApertures = numApertures, step = step, plotDir = plotDir, strips = strips, maskreg = mask, refimage = bgSubName) #, plot = True)
    
    ######################################################################
    # make a nice file with the output
    depthFile = resultsDir + '{0}_{1}.txt'.format(filterName, numApertures)
    f = open(depthFile, 'w')
    
    apString = ''
    for i in range(apDiametersAS.size):
        apString = apString + '{0:.1f}as\t'.format(apDiametersAS[i])

    f.write('#ap\t{0}\t{1}\n'.format(apString, 'type'))

    depthtype = ['median', 'global', 'mode']
    for di, deptht in enumerate(depthtype):
        
        for r, reg in enumerate(regions):
            
            apResultString = ''
            for i in range(apDiametersAS.size):
                if deptht == 'median':
                    apResultString = apResultString + '{0:.2f}\t'.format(meddepths[r, i])
                elif deptht == 'global':
                    apResultString = apResultString + '{0:.2f}\t'.format(globaldepths[r, i])
                elif deptht == 'mode':
                    apResultString = apResultString + '{0:.2f}\t'.format(modedepths[r, i])
                    
            printString = '{0}\t{1}\t{2}\n'.format(reg, apResultString, deptht)
            f.write(printString)

    f.close()
    print("Output file saved to ", depthFile)


    if coords != 'full-size':
        ra = coords.split('_')[0]
        dec = coords.split('_')[1]
        save_coordinates(catDir, ra, dec, verbose=True)
    return outputCatalogue


def extract_local_depths(inputTableFile, apDiametersAS, zeropoint, step = 500, numApertures = 200, strips = False, plot = False, local = True, recalculate = False, globalplot = False, clean = True, plotDir = '', maskreg = 'none', refimage = 'none'):
    ''' extractDepths.py
    
    Code to extract the depths from the input table of aperture photometry
    
    Modified: Dec 2019 '''

    import matplotlib
    #matplotlib.use('pdf')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    from matplotlib import gridspec
    from scipy.stats import norm
    from astropy import units as u
    import matplotlib.backends.backend_pdf
    from new_catalogue_codes import return_instrips, mask_column
    
    ###########################################
    # important setup
    nirspec = False
    edgebuffer = 0.1 ## 10% of the edge #500
    edgebuffer_full = 0 #3000
    # remove things close to zero!
    deltaZero = 1E-13
    
    # extract name for local depth results
    basename = inputTableFile[:-13]

    colourArray = ['Blue', 'Green', 'Green', 'Green', 'Green', 'Red', 'Red', 'Red', 'Red', 'Red']
    
    ## Now loop through the different regions of the image
    if strips:
        # for ultravista!
        regions = ['fullimage', 'stripone', 'striptwo', 'stripthree', 'stripfour', 'gap1', 'gap2', 'gap3', 'gap4']
        regions = ['full', 'str1', 'str2', 'str3', 'str4', 'gap1', 'gap2', 'gap3', 'gap4']
    else:
        regions = ['full']
        
    global_depth = np.zeros([len(regions), apDiametersAS.size])
    median_local_depth = np.zeros([len(regions), apDiametersAS.size])
    mode_local_depth = np.zeros([len(regions), apDiametersAS.size])
    
    # save plot
    startI = inputTableFile.rfind('/') + 1
    endI = inputTableFile.rfind('_')
    filterName = inputTableFile[startI:endI]
    plotName = plotDir + filterName + '_' + str(numApertures) + '_{0}.pdf'.format(step)
    pdf = matplotlib.backends.backend_pdf.PdfPages(plotName)
    
    # Loop through the available apertures
    for ai, apDiAS in enumerate(apDiametersAS):
        
        print("Extracting local and global depths for aperture = ", apDiAS, " as.")
        
        # The apertures are strings _0, _1 etc
        apString = '_{0:.1f}as'.format(apDiAS)
        
        localDepthsFile = basename + str(apDiAS) + 'as_gridDepths_{0}_{1}.fits'.format(numApertures, step)
  
        print("File will be saved to ", localDepthsFile)
        
        ## Always read in aperture file
        ## So i can do global depths
        
        ##############################################
        ## Read in the table file
        ## This has the aperture photometry (not the depths)
        inputTable = Table.read(inputTableFile)
        print("Reading file ", inputTableFile)
        
        if clean == False:
            # do the cleaning here
            smallNum = 0.0000001
         ## cut the table depending on apertures
         ## only accept coordinates where the seg map
         ## aperture is blank, and the wht map is > 0
            seg_sum = np.array(inputTable['SEG_flux' + apString])
            wht_sum = np.array(inputTable['WHT_flux' + apString])
            ap_sum = np.array(inputTable['IMAGE_flux' + apString])
            good_indicies = (seg_sum < 0.5)  & (wht_sum > smallNum) & \
                ((ap_sum > deltaZero) | (ap_sum < -deltaZero))
            #print("There are ", sum(good_indicies), " good indicies, out of ", len(seg_sum))
            good_indicies = np.where(good_indicies)
            
        ## This table has all the good apertures
            reducedTable = inputTable[good_indicies]
            
        else:
            reducedTable = inputTable
            #print("inputTable",inputTableFile)
        ## Get the x and y coordianates of the apertures
        apX = reducedTable['IMAGE_xcenter'] #/u.pix
        apY = reducedTable['IMAGE_ycenter'] #/u.pix
        
        ## Check if I need to run local depths from the
        ## apertures
        
        if (os.path.isfile(localDepthsFile) == False) or recalculate:
            ## Files doesn't exist or I want to recalculate it  
            
            ################## RUNNING LOCAL DEPTHS ####################

            xmax = max(apX)# - min(apX)#*3.0
            ymax = max(apY)# - min(apY)#*3.0
            
            numX = np.ceil(xmax/step)
            numY = np.ceil(ymax/step)

            x = min(apX) + np.arange(numX)*step ## modifed 21/9/2018 to add min(apX).
            y = min(apY) + np.arange(numY)*step
            #print("Step = ", step, " numx, y = ", numX, numY)
            #print("Max = ", xmax, ymax, max(apX), min(apX))
            
            # create x, y arrays
            x = np.zeros(1) #.value
            y = np.zeros(1) #.value                
            
            for xi in np.arange(step/2.0, numX*step, step):
                for yi in np.arange(step/2.0, numY*step, step):
                    x = np.append(x, xi)
                    y = np.append(y, yi)
                    
            # I want a constant grid over the image.
            
            # remove the first elements
            x = x[1:]
            y = y[1:]
            
            ## Now run local depths at those points
            depthsLocalFull, maskArray = local_depths(reducedTable, apString, x, y, numApertures, zeropoint = zeropoint, mask = True, sigmaClip = 3.0)#, plot = plotDir + filterName + '_' + str(numApertures))
            
            ## remove points that lie off the image
            
            ## Now save these results for faster calculating/plotting in future
            ## Create a table
            localTable = Table([x, y, depthsLocalFull, maskArray], names = ['x', 'y', 'depths', 'mask'], dtype = ['f4', 'f4', 'f4', 'f4'])
            localTable.write(localDepthsFile, overwrite=True)
            print("Local depths saved to ", localDepthsFile)
            
        else:
            # simply restore the results
            localTable = Table.read(localDepthsFile)
                
    # for plotting and median depths, remove negative objects!
        gg = (localTable['mask'] > 0)
        #gg = (localTable['depths'] > 0.0)
        localTable = localTable[gg]
        x = localTable['x']
        y = localTable['y']
        keys = localTable.colnames
        depthsLocalFull = localTable['depths']
        
    ################################################
    ## Extract median depths etc
    ## And for different sub-regions
    
    ## Loop through the different regions
        for ri, region in enumerate(regions):
        
            print("Calculating depths in region ", region)
                    
            if region != 'full':
                print('Splitting by strip')
                good_indicies = return_instrips(x, y, region = region)
            else:

                if maskreg == 'none':
                    ## put a buffer here
                    maxNumx = max(x)
                    maxNumy = max(y)
                    minNumx = min(x)
                    minNumy = min(y)
                    
                    minNumx = minNumx + edgebuffer*maxNumx
                    minNumy = minNumy + edgebuffer*maxNumy
                    maxNumx = maxNumx - edgebuffer*maxNumx
                    maxNumy = maxNumy - edgebuffer*maxNumy
                    
                    good_indicies = (x > minNumx) & (x < maxNumx) & (y > minNumy) & (y < maxNumy)
                    
                else:
                    # read in a header
                    from astropy.wcs import WCS
                    w = WCS(refimage)
                    ra, dec= w.all_pix2world(x,y, 1)
                    print('Masking with ', maskreg)
                                        
                    hsci = refimage.find('HSC')
                    
                    if hsci > -1:
                        hsc = True
                        print('Masking HSC')
                    else:
                        hsc = False

                    if dec[0] < -15.0:
                        # get the directory

                        kk = maskreg.rfind('/')
                        regDir = maskreg[:kk+1]
                        #print(regDir)
                        
                        fff = refimage.find('HSC-R')
                        if fff > -1:
                            circlesFile = regDir + 'HSC_circle_cdfs_R_xy.reg'
                        else:
                            circlesFile = regDir + 'HSC_circle_cdfs_xy.reg'
                            
                        good_indicies = mask_column(x, y, maskreg, tokeep = True, hsc = hsc, xy = True, circlesFile = circlesFile)
                    else:
                        #print(ra, dec, maskreg, hsc)
                        good_indicies = mask_column(ra, dec, maskreg, tokeep = True, hsc = hsc)
                    
                
                #print "There are ", x.shape, y.shape, " local depth positions..."

            #print('There are {0} good indicies'.format(np.sum(good_indicies)))

            good_indicies = good_indicies & np.logical_not(np.isnan(localTable['depths']))
            
            finalTable = localTable[good_indicies]
            apXregion = finalTable['x'] #/u.pix
            apYregion = finalTable['y'] #/u.pix
            depthsLocal = finalTable['depths'] #depthsLocalFull[good_indicies]

            ii = np.logical_not(np.isnan(depthsLocal))
            
            #print('Check for NANS', np.isnan(np.sum(depthsLocal)), np.sum(ii), depthsLocal[ii])
            
            # For the median depth of the local depths
            # I want to exclude the edges
            medianLocalDepth = np.median(depthsLocal)
            #           print(depthsLocal, medianLocalDepth)
            #            exit()
            magmin = medianLocalDepth - 1.2 # 1.2 before 19/05/22
            magmax = medianLocalDepth + 1.2 #1.8 before 16/10/23

            if strips:
                if filterName[0] == 'N':
                    magmin = 23.0
                    magmax = 25.5
                else:
                    magmin = 24.5
                    magmax = 27.0
                    
            if nirspec:
                magmin = medianLocalDepth - 0.9
                magmax = medianLocalDepth + 1.3
            
            if region == 'full':

                                
                fig = plt.figure(figsize=(6,8))
                gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
                ax = plt.subplot(gs[0])
                plt.axis('equal')

                if filterName[0] == 'N':
                    ax.set_xlim([-900, 38000])
                    ax.set_ylim([-4000, 34000])
                # get colours
                cm = plt.cm.get_cmap('RdYlBu')
                sc = plt.scatter(x, y, s = 5, c = depthsLocalFull, cmap = cm, linewidth = 0.0, vmin = magmin, vmax = magmax)

                plt.colorbar(sc)
                plt.title('Local depths for filter {0}\n Aperture diameter is {1:.1f}as'.format(filterName, apDiAS))
                
                # now make a histogram to go underneath, do this differently if in strips!
                ax = plt.subplot(gs[1])
                binwidth = 0.01 #0.01
                low = magmin
                
                bins = np.arange(low+binwidth/2.0, magmax+binwidth, binwidth)
                #print('Bins defined here ', bins)
                
                if strips:
                    ## plot subsets of the results
                    ## and colour based on strip/gap 
                    #print('In strips', bins.shape)
                    
                    ## 1) get the strips
                    #print("allgood")
                    ss = return_instrips(apXregion, apYregion)
                    #print("allgood2")
                    strip_histy, strip_histx, _ = plt.hist(depthsLocal[ss], facecolor= 'blue', alpha = 0.8, density = False, range = [magmin, magmax], bins = bins, label = ' Ultra-deep')

                    gg = return_instrips(apXregion, apYregion, region = 'puregap')# notinstrips       
                    gap_histy, gap_histx, _ = plt.hist(depthsLocal[gg], facecolor= 'red', alpha = 0.8, density = False, range = [magmin, magmax], bins = bins, label = ' Deep')
                    
                    # Make a caption
                    plt.legend(loc = 'upper left', handletextpad = 0.0, fontsize = 8, frameon = False)
                    
                    N  = 10
                    strip_smoothed = np.convolve(strip_histy, np.ones((N,))/N, mode = 'same')
                    plt.plot(bins[:-1]+binwidth/2.0, strip_smoothed, color= 'k', linestyle = ':')
                    gap_smoothed = np.convolve(gap_histy, np.ones((N,))/N, mode = 'same')
                    plt.plot(bins[:-1]+binwidth/2.0, gap_smoothed, color= 'k', linestyle = ':')

                    smoothed = strip_smoothed
                    #print('Definiting smoothed here, bins ', smoothed.shape, bins.shape)
                    
                else:

                    #print('Strips = false?')
                    #print(depthsLocal)
                    #print(magmin, magmax)
                    bad = True
                    if bad:
                        
                        histy, hists, _ = plt.hist(depthsLocal, facecolor= colourArray[ri], alpha = 0.8, density= True, range = [magmin, magmax], bins =bins, histtype = 'step')
                        hx, hy, _ = plt.hist(localTable['depths'], facecolor= 'grey', alpha = 0.8, density = True, range = [magmin, magmax], bins = bins, zorder = 1, histtype = 'step')
                        
                        # smooth this
                        # and plot again
                        N = 10
                        smoothx = bins[:-1]+binwidth/2.0
                        smoothed = np.convolve(histy, np.ones((N,))/N, mode = 'same')
                        plt.plot(smoothx, smoothed, color= 'k', linestyle = ':')
                        print('In bad loop')
                                   
            #####################################################
            ######## GLOBAL depths #############
            apertureResults = np.array(reducedTable['IMAGE_flux'+apString])
                        
            ## First clip the data to remove outliers
            ## that will skew the SD
            medianFlux = np.median(apertureResults)
            mad = np.median(abs(apertureResults - medianFlux))
            sigma_mad = 1.4826*mad
            sigmaClip = 2.0
            #print "The mad is ", mad, medianFlux
            
            ## cut the results for the FIRST time
            good_indicies = (apertureResults > medianFlux - sigmaClip*sigma_mad) & \
                            (apertureResults < medianFlux + sigmaClip*sigma_mad)
            
            cutTable = reducedTable[good_indicies]
            apertureResultsTwo = np.array(cutTable['IMAGE_flux'+apString])
            
            medianFlux = np.median(apertureResultsTwo)
            mad = np.median(abs(apertureResultsTwo - medianFlux))
            sigma_mad = 1.4826*mad
            
            ## cut the results for the SECOND time
            good_indicies = (apertureResultsTwo > medianFlux - sigmaClip*sigma_mad) & \
                            (apertureResultsTwo < medianFlux + sigmaClip*sigma_mad)
            
            cutTableTwo = cutTable[good_indicies]
            reducedResults = np.array(cutTableTwo['IMAGE_flux'+apString])

            ## Plot a histogram of the results.. is it a gaussian?!
            n, binsss = np.histogram(reducedResults, 20, density= True)

            ## Fit with Gaussian
            (mu, sigma) = norm.fit(reducedResults)
            
            #print "The global depth is then: MAD = {0:.2f}, GAUSSIAN = {1:.2f} from {2} apertures.".format(return_mag(mad, zeropoint, sigma = 5.0), return_mag(sigma, zeropoint, sigma = 5.0), reducedResults.shape)
            
            medianFlux = np.median(reducedResults)
            mad = np.median(abs(reducedResults - medianFlux))
            sigma_mad = 1.4826*mad
            global_depth[ri, ai] = return_mag(sigma_mad, zeropoint, sigma = 5.0)
            
            #print "The global depth is ", sigma, zeropoint, global_depth[ri,ai]

            ####################################################
            # Median local depth
            #print("The median local depth is then: MAD = {0:.2f}, Gaussian = {1:.2f}".format(medianLocalDepth, global_depth[ri, ai]))
            median_local_depth[ri, ai] = medianLocalDepth

            if region == 'full':
                # Mode depth
                #            mod = (histy == np.max(histy))
                #mod = (smoothed == np.max(smoothed))
                mod = np.argmax(smoothed)
                smallbins = bins[0:-1]
                #print(smoothed.shape, smallbins.shape, bins.shape, mod.shape)
                #print(bins)
                #print('smallbinds ', smallbins)
                #print(smallbins[mod])
                
                mode = smallbins[mod]+binwidth/2.0
                #print(mode)
                #mode = mode[0]
                mode_local_depth[ri, ai] = mode
            
            ## add these depths to the plot
            if region == 'full':
                if strips == False:
                    
                ## plot this global depth over the top!
                   
                    ylim = ax.get_ylim()

                    depthString = 'Mode depth = {0:4.2f}'.format(mode)
                    plt.text(magmin,ylim[1]*0.9, depthString, fontsize = 10, color = 'k')

                    if nirspec == False:
                        plt.plot([mode, mode], ax.get_ylim())
                        plt.plot([global_depth[ri, ai],global_depth[ri, ai]], ax.get_ylim(), 'r') 
                        plt.plot([medianLocalDepth, medianLocalDepth], ax.get_ylim(), 'g-', linewidth = 2.0) 
                        depthString = 'Median local depth = {0:4.2f}'.format(medianLocalDepth)
                        plt.text(magmin, ylim[1]*0.8, depthString, fontsize = 10)
                        depthString = 'Global depth = {0:4.2f}'.format(global_depth[ri, ai])
                        plt.text(magmin,ylim[1]*0.7, depthString, fontsize = 10)
                        # print('Mode = ', mode)
                    else:

                        from scipy.signal import find_peaks
                        peaks,_ = find_peaks(smoothed, width = 10)
                        #print(peaks)
                        #print('Smoothed!', smoothx[peaks])
                        plt.text(magmin+0.8, ylim[1]*0.9, 'Peak depths = {0}'.format(smoothx[peaks]))
                        plt.scatter(smoothx[peaks], smoothed[peaks], color = 'red')
                        
                pdf.savefig(fig)
                plt.close()
                
            ## also plot the results from global depths
            if globalplot:
                fig2 = plt.figure(figsize=(6,8))
                n, bins, patches = plt.hist(reducedResults, 20, density = True)
                
                #yll = mlab.normpdf(bins, mu, sigma)
                from scipy.stats import norm
                yll = norm(loc = mu, scale = sigma)
                l = plt.plot(bins, yll.pdf(bins), 'r--', linewidth=1)
                plt.xlabel('Flux/counts')
                plt.ylabel('Number of apertures')
                plt.title('Ap diameter  =' + str(apDiAS) + ' region = '+ region)
                pdf.savefig(fig2)
                                
        plt.close()
        ## loop through the different regions.

    # close the plot
    pdf.close()
    print("Plot saved to ", plotName)
    
    ## Before I return the results
    ## Collect the results into a nice table
    ## with ok formatting!
    for x, xi in enumerate(regions): print('{0}, {1:.2f}'.format(xi, global_depth[x,0]))
    
    return regions, global_depth, median_local_depth, mode_local_depth


def return_mag(flux_counts, zeropoint, sigma = 1.0):
    import math
    
    if flux_counts < 1E-15:
        return -99.0
    else:
        return -2.5*math.log(sigma*flux_counts, 10.0) + zeropoint

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


def get_depths(fieldName, back_size, back_filtersize, coords, queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir=baseDir+'data/data_test/', outputDir='none', overwrite=False):
    
    # set the grid sep
    gridSepAS = 3.0
        
    # Read in the images.lis file
    dirHere = dataDir + fieldName + '/'
    dirHere = dataDir + 'COSMOS' + '/'    # holly - change this
    imagedata = read_image_lis(dirHere)
    availableFilters = np.array(imagedata['Name'])
    #print("The available filters are:\n", availableFilters)
    #print("Running get_depths", dirHere)
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
        whtName = imagedata['Weight'][fi]

        if coords == 'full-size':
            # use original data files so that header is as expected
            imageName = imageName.split('.')[0]+'_backup.fits' 
            whtName = whtName.split('.')[0]+'_backup.fits'
        else:
            imageName = get_cutout_img(imageName, coords)
            whtName = get_cutout_img(whtName, coords)


        whtType = imagedata['Wht_type'][fi]
        zeropoint = imagedata['zeropoint'][fi]
        imageDir = imagedata['directory'][fi]
        maskName = '/raid/scratch/data/masks/{0}/{1}'.format(fieldName, imagedata['mask'][fi])
        
        if imageDir == 'here':
            if coords != 'full-size':
                imageDir = dataDir + fieldName + '/random_cutouts/'
            else:
                imageDir = dataDir + fieldName + '/'
                imageDir = dataDir + 'COSMOS/' #holly -  might cause break
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
            #print("queue =! 'none': ", queue)
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
    print(f"Objects in unfiltered catalogue {subset_tab_name}: \n {len(table)}")
    
    if subsetCriterion != 'none':
        subset = table[table[subsetCriterion]]  # criterion by which to create a subset
    else:
        print("No subset criterion supplied. Using default: ['MAG_APER'][:, 1] > 50 ")
        subset = table[table["MAG_APER"][:, 1] > 50 ]  # criterion by which to create a subset
    
    outputPath = os.path.join(outputDir, subset_tab_name)  
    subset.write(outputPath, overwrite=True)
    print(f"Objects in filtered sub-catalogue {subset_tab_name}: \n {len(subset)}")
    print(f"Saved subset to {outputPath}")

    if len(subset) == 0:
        bad_ratio = 0.0000000000001
    else:
        bad_ratio = len(subset)/len(table)

    return outputPath, bad_ratio

def load_image(filepath):
    """Load FITS image and return the data array."""
    with fits.open(filepath) as hdul:
        data = hdul[0].data
    return data


############################ PLOTTING #######################################################

def detections_fig(cutoutPath, detPath, filters, detCats, subsetCriterion='none'):
    """ 
    Makes a figure showing the science image, detections,
     and subset detections for two filters in a 3x2 formation. Subset creation happens here. 

        cutoutPath(str): /raid/scratch/data/COSMOS_test/random_cutouts/UVISTA_JH_DR6_15007186924767774_26432262759966996.fits"""

    import matplotlib
    matplotlib.use('TkAgg')  # or 'Qt5Agg', depending on your system
    
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
        subset_file, bad_ratio = make_subset(catPath, detPath, '_bad_mags.fits', subsetCriterion=subsetCriterion)
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
    plt.show()

    return matching_coords


def open_cats(catDir, filters, show_all=True, size_arcsec=None, bksz=None, bkfilt=None, ra=None, dec=None):
    """
        Opens catalogues for filters of interest in TOPCAT.


    INPUT(s)
        catDir(str)         lol
        filter(list[str])   Filters you want to see Cat for
        show_all(T/F)       Show all cats for the specified filter(s)
        size_arcsec(int)    The size of the science image cutout
    OUTPUT(s)

    """
    import re
    if (ra or dec) != None:
        ra = str(ra).replace('.', '')
        dec = str(dec).replace('.', '')

    all_cats = os.listdir(catDir) # all cats in this dir

    #if size_arcsec is not None:
        #show_all = False

    if show_all:
        # open all cats for filt? i.e. bgmap, seg...
        cats = [] 
        for f in filters:
            search_str = f
            for cat in all_cats:
                cat = os.path.join(catDir, cat)
                if search_str in cat and cat not in cats: # don't add the same cat multiple times
                    cats.append(cat)


    elif show_all==False:
        if size_arcsec == 'full-size':
            full_cats = []
            for f in filters:
                for cat in all_cats:
                    full_path = os.path.join(catDir, cat)
                    pattern = rf"d{f}bksz{bksz}bkfilt{bkfilt}_full\-size\.fits$"
                    subset_pattern = rf"d{f}bksz{bksz}bkfilt{bkfilt}_full\-size_bad_mags\.fits$"
                    basedOnOrig = rf"d{f}\.fits"

                    catPat = re.search(pattern, cat)
                    subsetCatPat = re.search(subset_pattern, cat)
                    basedOnOrigPat = re.search(basedOnOrig, cat)

                    if catPat is not None:
                        if full_path not in full_cats:
                            full_cats.append(full_path)
                    elif subsetCatPat is not None:
                        if full_path not in full_cats:
                            full_cats.append(full_path)
                    elif basedOnOrigPat is not None:
                        if full_path not in full_cats:
                            full_cats.append(full_path)
                    else:
                        continue

        if size_arcsec != None or bksz != None or bkfilt != None or ra != None or dec != None:
            cats = []
            #criteria for opening specific cats include: size_arcsec 
            for f in filters:
                for cat in all_cats:
                    full_path = os.path.join(catDir, cat)
                    # match: starts with d{filter}bksz..., ends with {size}arcsec.fits
                    #pattern = rf"d{f}bksz{bksz}bkfilt{bkfilt}_{ra}_{dec}_{size_arcsec}arcsec\.fits$" # in an earlier version I included the size in filename
                    pattern = rf"d{f}bksz{bksz}bkfilt{bkfilt}_{ra}_{dec}\.fits$"
                    subset_pattern = rf"d{f}bksz{bksz}bkfilt{bkfilt}_{ra}_{dec}_bad_mags\.fits$"
                    hello = re.search(pattern, cat)
                    goodbye = re.search(subset_pattern, cat)
                    if hello is not None:
                        filt, bksz, bkfilt, filename_ra, filename_dec = extract_filename(cat)
                        if (filename_ra == ra and filename_dec == dec) and (full_path not in cats):
                            cats.append(full_path)
                    elif goodbye is not None:
                        filt, bksz, bkfilt, filename_ra, filename_dec = extract_filename(cat)
                        if (filename_ra == ra and filename_dec == dec) and (full_path not in cats):
                            cats.append(full_path)
                    else:
                        continue

        # join all filenames into one space-separated string for shell command
        files_str = " ".join(cats)
        files_str = " ".join(full_cats)
        command = f"topcat {files_str} &"
        print(command)
        os.system(command)


def show_in_ds9_save(coords='full-size', param_combo=[32,9], filters=['J', 'JH'],
                images=['science', 'bgsub', 'bgmap'], 
                dataDir='/raid/scratch/data/data_test/COSMOS_test/random_cutouts/', 
                depthDir='/raid/scratch/depths/COSMOS_test/images/', verbose=True):
    sciImgs = []
    ds9_path = '/nvme/scratch/ullyott/ds9'
    command = f"{ds9_path}"
    if coords == 'full-size':
        dataDir = '/raid/scratch/data/data_test/COSMOS_test/'
    else:
        coords = coords.replace(', ', '_')

    matchSciImgs, matchBgImgs = find_matching_coord_files(
        param_combo=param_combo, coords=coords, dataDir=dataDir, depthDir=depthDir)

    for filt in filters:
        if not any(x in images for x in ['science', 'bgsub', 'bgmap']):
            continue

        # Get science image for this filter
        sciImg = next((img for img in matchSciImgs if f"UVISTA_{filt}_" in img), None)
        if sciImg is None:
            print(f">>>> No science image found for filter {filt}")
            continue

        if not sciImg.startswith('/'):
            sciImg = os.path.join(dataDir, sciImg)
            sciImgs.append(sciImg)
        #command += f" {sciImg} -zscale"

        # Get both bgsub and bgmap images for this filter
        bgImgs = [
            img for img in matchBgImgs
            if filt in img and coords in img and (
                img.endswith('_bgsub.fits') or img.endswith('_bgmap.fits'))
        ]

        if len(bgImgs) < 2:
            print(f">>>> Warning: expected 2 bg images for {filt}, found {len(bgImgs)}")
            continue

        #for bg in sorted(bgImgs):
            #command += f" -frame new {bg} -zscale"

    # Finalize and run
    for sciImg in sciImgs:
        command += f" {sciImg} -zscale"
        for bg in sorted(bgImgs):
            command += f" -frame new {bg} -zscale"
        
    command += " -tile grid layout 3 2 -zoom to fit -frame lock wcs &"

    if verbose:
        print("Running DS9 command:")
        print(command)

    os.system(command)

def show_in_ds9(coords='full-size', param_combo=[32,9], filters=['J', 'JH'],
                images=['science', 'bgsub', 'bgmap', 'seg'], 
                dataDir='/raid/scratch/data/data_test/COSMOS_test/random_cutouts/', 
                depthDir='/raid/scratch/depths/COSMOS_test/images/', verbose=True):

    sciImgs = []

    ds9_path = '/nvme/scratch/ullyott/ds9'
    command = f"{ds9_path}"
    total_frames = 0

    if coords == 'full-size':
        bgImgs_full = []
        dataDir = '/raid/scratch/data/COSMOS/'
        for f in os.listdir(dataDir):
            if f.startswith("UVISTA_") and f.endswith('_DR6_backup.fits'):
                for filt in filters:
                    f = dataDir+"UVISTA_" + filt + '_DR6_backup.fits'
                    sciImgs.append(f)
                    bgImg = filt + '_bgsub.fits'
                    segImg = filt + '_seg.fits'
                    bgImgs_full.append(depthDir+bgImg) 
                    bgImgs_full.append(depthDir+segImg)
                    command += f" {f} -zscale"
                    total_frames += 1
                    command += f" -frame new {bgImg} -zscale"
                    total_frames += 1
                    command += f" -frame new {segImg} -zscale"
                    total_frames += 1

    else:
        coords = coords.replace(', ', '_')

        matchSciImgs, matchBgImgs = find_matching_coord_files(
        param_combo=param_combo, coords=coords, dataDir=dataDir, depthDir=depthDir)

        for filt in filters:
            # Get science image
            if 'science' in images:
                sciImg = next((img for img in matchSciImgs if f"UVISTA_{filt}_" in img), None) or next((img for img in sciImgs if f"UVISTA_{filt}_" in img), None)
                if sciImg:
                    if not sciImg.startswith('/'):
                        sciImg = os.path.join(dataDir, sciImg)
                    command += f" {sciImg} -zscale"
                    total_frames += 1
                else:
                    print(f">>>> No science image found for filter {filt}")
                    continue

            # Get background-related images
            bgImgs = [
                img for img in matchBgImgs
                if filt in img and coords in img and any(t in img for t in ['_bgsub.fits', '_bgmap.fits', '_seg.fits'])
            ]

            for tag in ['bgmap', 'bgsub', 'seg']:
                if tag in images:
                    match = next((img for img in bgImgs if tag in img), None)
                    if match:
                        print("match:", match)
                        if not match.startswith('/'):
                            match = os.path.join(depthDir, match)
                        command += f" -frame new {match} -zscale"
                        total_frames += 1
                    else:
                        print(f">>>> Missing {tag} image for filter {filt}")

    # Compute layout
    import math
    columns = len(images)
    rows = math.ceil(total_frames / columns)
    # {columns} {rows}
    command += f" -tile grid layout 3 1 -zoom to fit -frame lock wcs &"

    if verbose:
        print("Running DS9 command:")
        print(command)

    os.system(command)

def bkg_plotter(param_combo, filters=['J', 'JH'], coords='full-size', depthDir='/raid/scratch/depths/COSMOS_test/images/', dataDir='/raid/scratch/data/data_test/COSMOS_test/', verbose=True):
    """
    what affect do each set of SE bkg modelling parameters have?
    Do visual check with this side-by-side plotter

    param_combos(list): SE parameters used in bkg modelling test: [BACK_SIZE, BACK_FILTERSIZE]
    coords(list):       coordinates of random science cutout to match. [RA/Dec]
                        If default, mathces to full image (feature coming soon heheh) 
    depthDir(str):      directory of SE'd find_depth catalogues
    dataDir(str):       directory of science images (full-size and cutouts)
    """
    # get files
    if coords=='full-size':
        matchSciImgs, matchBgImgs = find_fullsize_files(param_combo=param_combo, dataDir=dataDir, depthDir=depthDir, verbose=verbose)
    else:
        matchSciImgs, matchBgImgs = find_matching_coord_files(param_combo=param_combo, dataDir=dataDir, depthDir=depthDir, coords=coords, verbose=verbose)

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
            detCats.append(det_table)
            det_name = os.path.splitext(os.path.basename(detCats[f]))[0]
            catNames.append(det_name)
        elif os.path.isfile(detCats[f])==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {detCats[f]}")
        
        ## get subset data
        subset_dir, bad_ratio = make_subset(detCats[f], baseDir+'depths/COSMOS_test/catalogues/', '_bad_mags.fits', subsetCriterion=subsetCriterion)
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

    for i, (sci_img, det_table, det_name, subset_table, subset_name) in enumerate(zip(sciImgs, catTables,  catNames, subsetTables, subsetNames)):
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



