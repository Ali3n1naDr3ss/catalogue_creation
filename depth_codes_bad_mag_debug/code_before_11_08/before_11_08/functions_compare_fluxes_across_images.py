# visually compare the J, H, JH, and HK fluxes of a small sample of "bad mag" objects and store a Table with these objects
# 31/01/2025

import os
import re
import sys
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table

outputCatalogues = []

def get_fits_data(fits_img):
    with fits.open(fits_img) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
    return data, header, wcs

def save_coordinates(directory, ra, dec, filename='random_coords.txt', verbose=True):
    # Ensure the directory exists
    os.makedirs(directory, exist_ok=True)
    
    # Full path to the file
    filepath = os.path.join(directory, filename)
    
    # Line to be added
    new_line = f"{ra}, {dec}"
    
    # Create the file if it doesn't exist
    if not os.path.exists(filepath):
        with open(filepath, 'w') as f:
            f.write("RA, Dec\n")
    
    # Read the file and check if the line already exists
    with open(filepath, 'r') as f:
        lines = f.read().splitlines()

    # if new line not already in file, write line
    if new_line not in lines:
        with open(filepath, 'a') as f:
            f.write(f"{new_line}\n")
        if verbose:
            print(f"Saved coordinates: RA = {ra}, Dec = {dec} to {filepath}\n")
    else:
        if verbose:
            print(f"Coordinates already exist: RA = {ra}, Dec = {dec}\n")



def make_cutout(dataDir, filters, subDir, size_arcsec=1000, ra='random', dec='random', corner='lower-left',verbose=True):
    """
    Makes square science and weight image cutouts, given the RA/Dec of a specified corner, and saves to data/{fieldName}/randomcutouts by default.

    ra, dec(float):     RA and Dec (in degrees) of one corner of the square.
    size_arcsec(float): Size of the square cutout (arcsec on a side).
    corner(str):        Which corner the (ra, dec) point is (options: 'lower-left', 'upper-left', 'lower-right', 'upper-right').
    subDir(str):        output path = dataDir + subDir
    """

    import astropy.units as u
    from astropy.nddata import Cutout2D
    from astropy.coordinates import SkyCoord

    imgPaths = []
    outputImgNames = []
    outputWhtNames = []
    writtenFiles = []

    if (ra and dec) == 'random':
        subDir = 'random_cutouts/'
    else:
        subDir == subDir

    for f, filt in enumerate(filters):
        imgName = 'UVISTA_' + filt + '_DR6_backup.fits'
        whtName = 'UVISTA_' + filt + '_DR6_wht_backup.fits'
        imgPath = os.path.join(dataDir, imgName)
        whtPath = os.path.join(dataDir, whtName)
        if os.path.isfile(imgPath):
            imgPaths.append(imgPath)
        if os.path.isfile(whtPath):
            imgPaths.append(whtPath) # add wht to imgs, don't need to differentiate
    print("weight:",imgPaths)
    for i, img in enumerate(imgPaths):
        data, header, wcs = get_fits_data(imgPath)
        # get coords of cutout
        if (ra and dec) == 'random':
            ra, dec = make_random_coords(img, verbose=verbose)
            
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

        # make output path
        ra_filename = str(ra).replace('.', '')
        dec_filename = str(dec).replace('.', '')
        coords_filename = ra_filename + '_' + dec_filename
        img = os.path.basename(img).replace('_backup.fits', '')
        CutoutDir = dataDir + subDir

        # check path exists, make it if necessaray
        if os.path.isdir(CutoutDir) == False:
            os.system('mkdir ' + CutoutDir)
        
        output_img_name = f'{CutoutDir}{img}_{coords_filename}_{size_arcsec}arcsec.fits'


        save_coordinates(CutoutDir, ra_filename, dec_filename, subDir[:-1]+'.txt')
        hdul_out.writeto(output_img_name, overwrite=True)
        writtenFiles.append(output_img_name)
        if verbose:
            print(f"Saved image cutout to: {output_img_name}\n")

        if 'wht' in img:
            output_wht_name = f'{CutoutDir}{img}_{coords_filename}_{size_arcsec}arcsec.fits'
            hdul_out.writeto(output_wht_name, overwrite=True)
            writtenFiles.append(output_wht_name)
            if verbose:
                print(f"Saved weight cutout to: {output_wht_name} \n")


    print("Saved files: \n", writtenFiles)

    return outputImgNames, outputWhtNames




def image_depth(imageName, zeropoint, back_size, back_filtersize, coords, size_arcsec, apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtName = 'NONE', whtType = 'NONE', IRACapDiametersAS = np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segName = 'NONE', outputDir = 'none', filterName = 'NONE', numApertures = 300, step = 200, overwrite = False, inputSex = 'raid/scratch/data/bertin_config/video_mine.sex', strips = False, bgSub = True, mask = 'none', gridSepAS = 3.0):
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

        coords = str(coords[0]).replace('.','') + '_' + str(coords[1]).replace('.','') # get coords in correct format for filenaming

        inputSex = '/raid/scratch/data/data_test/COSMOS_flux_comp/video_mine.sex'

        segName = imageDir + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords +'_'+  size_arcsec + 'arcsec_seg.fits'  
        bgSubName = imageDir  + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords +'_'+  size_arcsec + 'arcsec_bgsub.fits'
        bgMapName = imageDir  + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords +'_'+  size_arcsec + 'arcsec_bgmap.fits'
        outputCatalogue = catDir + 'd' + filterName + 'bksz' + str(back_size) + 'bkfilt' + str(back_filtersize) +'_'+  coords +'_'+  size_arcsec + 'arcsec.fits'

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


def get_depths(back_size, back_filtersize, coords, dataDir, size_arcsec, queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), outputDir='none', overwrite=False):
    
    size_arcsec = str(size_arcsec)
    # set the grid sep
    gridSepAS = 3.0
        
    # Read in the images.lis file
    dirHere = dataDir
    print("Running get_depths", dirHere)

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
    
    def get_cutout_img(imageName, coords, size_arcsec):
        coords = str(coords[0]).replace('.','') + '_' + str(coords[1]).replace('.','')
        imageName = imageName.replace(".fits", '_'+coords+'_'+size_arcsec+"arcsec.fits")
        imageName = 'flux_comparison_test/' + imageName
        return imageName

    ### Run image_depth for each filter ###
    # Run each filter as a separate file...
    for fi, filterName in enumerate(availableFilters):
        # define the images etc to send through
        imageName = imagedata['Image'][fi]
        imageName = get_cutout_img(imageName, coords, size_arcsec)
        whtName = imagedata['Weight'][fi]
        whtName = get_cutout_img(whtName, coords, size_arcsec)
        whtType = imagedata['Wht_type'][fi]
        zeropoint = imagedata['zeropoint'][fi]
        imageDir = imagedata['directory'][fi]
        maskName = '/raid/scratch/data_test/masks/{0}/{1}'.format('COSMOS_flux_comp', imagedata['mask'][fi])

        if imageDir == 'here':
            imageDir = dataDir #+ 'flux_comparison_test/'

        # set up for image_depth
        if queue == 'none':
            #print("No queue. \n")

            strips = True
            if 'COSMOS' in 'COSMOS_flux_comp':
                jj = (filterName == stripFilt) # VISTA observation strategy produces images with 'strips' with alternating depth
                if np.any(jj):
                    strips = True
                    
            outputCatalogue =  image_depth(imageDir + imageName, zeropoint, back_size, back_filtersize, coords, size_arcsec, whtName = imageDir + whtName, whtType = whtType, outputDir = outputDir, strips = strips, filterName = filterName, overwrite = overwrite, mask = maskName, gridSepAS = gridSepAS, apDiametersAS = apDiametersAS)
        
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
            tmpName = "tmp_{1}_{0}.sh".format(filterName, COSMOS_flux_comp)
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

# get aperture fluxes table - then we can narrow it down into just a few objects and look at the numerical data (neg fluxes?)


def extract_filename(filename, verbose=False):
    """Extract substrings from various UVISTA file naming patterns."""
    import re
    base = os.path.basename(filename)

    patterns = [
        # cutouts
        (r"(UVISTA_)(\w+)(_DR6_)(\d+)_(\d+)\.fits",
                                              (1, 2, 3, 4, 5, None)),
        (r"(UVISTA_)(\w+)(_DR6_wht_)(\d+)_(\d+)\.fits",
                                              (1, 2, 3, 4, 5, None)),

        # cutouts with size_arcsec
        (r"(UVISTA_)(\w+)(_DR6_)(\d+)_(\d+)_(\d+arcsec)\.fits",
                                              (1, 2, 3, 4, 5, 6)),
        (r"(UVISTA_)(\w+)(_DR6_wht_)(\d+)_(\d+)_(\d+arcsec)\.fits",
                                              (1, 2, 3, 4, 5, 6)),

        # full image
        (r"(UVISTA_)(\w+)(_DR6\.fits)",
                                              (1, 2, 3, None, None, None)),

        # detection catalogues with size_arcsec
        (r"d(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)_(\d+arcsec)\.fits",
                                              (1, 2, 3, 4, 5, 6)),

        # detection catalogues
        (r"d(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)",
                                              (1, 2, 3, 4, 5, None)),


        # background + seg (with or without arcsec size, in any order)
        (r"(\w+)bksz(\d+)bkfilt(\d+)_(\d+)_(\d+)_"
         r"(?:(?:\d+arcsec_)?(?:bg(?:sub|map)|seg)|"
         r"(?:bg(?:sub|map)|seg)(?:_(\d+)arcsec)?)\.fits",
                                              (1, 2, 3, 4, 5, 6)),
    ]

    for pattern, group_map in patterns:
        match = re.match(pattern, base)
        if match:
            # Always return a tuple of length 6, filling missing with None
            result = []
            for idx in group_map:
                result.append(match.group(idx) if idx else None)
            return tuple(result)

    if verbose:
        print(f">>>>>>>>>>>>>>>>>> WARNING: Filename doesn't match expected pattern: \n {filename} \n {base} \n")
    return None



def find_matching_coord_files(coords, param_combo, size_arcsec, dataDir, depthDir, verbose=True):

    size_arcsec = str(size_arcsec)

    # get all images from dataDir/depthDir that have matching ra/dec, and matching param_combo
    dataFiles = []
    bgFiles = []
    matchSciImgs = []
    matchBgImgs = []

    # get all the files in data and depth dirs
    all_data_files  = []
    for f in os.listdir(dataDir):
        f = os.path.join(dataDir, f)
        if os.path.isfile(f):
            all_data_files.append(f)

    all_bg_imgs  = []
    for f in os.listdir(depthDir):
        f = os.path.join(depthDir, f)
        if os.path.isfile(f):
            all_bg_imgs.append(f)

    # narrow down the list of files so you only have the ones of interest
    for file_ in all_data_files:
        file_ = os.path.basename(file_)
        result = extract_filename(file_)
        if result is None:
            continue  # Skip non-matching files
        else:
            dataFiles.append(file_)
    if verbose:
        print("dataFiles: ", dataFiles)

    for file_ in all_bg_imgs:
        file_ = os.path.basename(file_)
        result = extract_filename(file_)
        if result is None:
            continue  # Skip non-matching files 
        else:
            bgFiles.append(depthDir+file_)
        if verbose:
            print("bgFiles: ", bgFiles)

    # retireve data file with matching coords, size_arcsec
    for file_ in dataFiles:
        result = extract_filename(file_)
        if result is None or len(result) != 6:
            continue  # Skip non-matching files
        else:
            prefix, filt, suffix, ra, dec, arcsec = extract_filename(file_) # get bits to match to coords
            matchSciFname = 'UVISTA_' + filt + '_DR6_' + coords + '_'+size_arcsec + 'arcsec.fits'
            #matchSciFname_fullsize = 'UVISTA_' + filt + '_DR6_' + coords + '.fits'
            if file_ == matchSciFname:
                matchSciImgs.append(matchSciFname)
        #else:
            #if verbose:
                #print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: Check coord format; No matching coords for:", matchSciFname, '\n', file_)
   
    # retireve background file with matching coords and back_size, back_filtersize
    for file_ in bgFiles:
        result = extract_filename(file_)
        if result is None or len(result) != 6:
            continue  # Skip non-matching files
        else:
            filt, bksz, bkfilt, ra, dec, arcsec= extract_filename(file_) # get bits to match to coords/param_combo
            matchBg = filt + 'bksz' +str(param_combo[0]) + 'bkfilt'+ str(param_combo[1]) + '_'+ coords + '_'+size_arcsec 
            if matchBg in file_:
                matchBgImgs.append(file_)
        #else:
            #if verbose:
                #print(">>>>>>>>>>>>>>>>>>>>>>>> WARNING: No matching params for:", matchBg, '\n', file_)

    return matchSciImgs, matchBgImgs 

def show_in_ds9(coords='full-size', param_combo=[32,9], filters=['J', 'JH'], size_arcsec=None,
                images=['science', 'bgsub', 'bgmap'], 
                dataDir='/raid/scratch/data/COSMOS_test/random_cutouts/', 
                depthDir='/raid/scratch/depths/COSMOS_test/images/', verbose=True):

    ds9_path = '/nvme/scratch/ullyott/ds9'
    command = f"{ds9_path}"
    if coords == 'full-size':
        dataDir = '/raid/scratch/data/COSMOS_test/'
    else:

        ra_filename = str(coords[0]).replace('.', '')
        dec_filename = str(coords[1]).replace('.', '')
        coords = ra_filename + '_' + dec_filename


    matchSciImgs, matchBgImgs = find_matching_coord_files(
        param_combo=param_combo, coords=coords, size_arcsec=size_arcsec, dataDir=dataDir+'flux_comparison_test/', depthDir='/raid/scratch/'+depthDir)

    for filt in filters:
        if not any(x in images for x in ['science', 'bgsub', 'bgmap']):
            continue

        # Get science image for this filter
        sciImg = next((img for img in matchSciImgs if f"UVISTA_{filt}_" in img), None)
        if sciImg is None:
            print(f">>>> No science image found for filter {filt}")
            continue

        if not sciImg.startswith('/'):
            sciImg = os.path.join(dataDir, 'flux_comparison_test/', sciImg)

        command += f" {sciImg} -zscale"

        # Get both bgsub and bgmap images for this filter
        bgImgs = [
    img for img in matchBgImgs
    if f"{filt}bksz" in img and coords in img and (
        img.endswith('_bgsub.fits') or img.endswith('_bgmap.fits') or img.endswith('_seg.fits'))]

        if len(bgImgs) < 3:
            print(f">>>> Warning: expected 3 bg images for {filt}, found {len(bgImgs)}")
            continue

        for bg in sorted(bgImgs):
            command += f" -frame new {bg} -zscale"

    # Finalize and run
    command += " -tile grid layout 4 2 -zoom to fit -frame lock wcs &"

    if verbose:
        print("Running DS9 command:")
        print(command)

    os.system(command)

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

    if size_arcsec is not None:
        show_all = False

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
        if size_arcsec != None or bksz != None or bkfilt != None or ra != None or dec != None:
            cats = []
            #criteria for opening specific cats include: size_arcsec 
            for f in filters:
                for cat in all_cats:
                    full_path = os.path.join(catDir, cat)
                    # match: starts with d{filter}bksz..., ends with {size}arcsec.fits
                    pattern = rf"d{f}bksz{bksz}bkfilt{bkfilt}_{ra}_{dec}_{size_arcsec}arcsec\.fits$"
                    hello = re.search(pattern, full_path)
                    if re.search(pattern, full_path) and full_path not in cats:
                        cats.append(full_path)

    # join all filenames into one space-separated string for shell command
    files_str = " ".join(cats)
    command = f"topcat {files_str} &"
    os.system(command)


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




def detections_fig(cutoutPath, detPath, filters, detCats, size_to_match, cutout_coords_file='random_coords.txt',  detCat_coords_file='random_coords.txt' , subsetCriterion='none'):
    """ 
    Makes a figure showing the science image, detections,
     and subset detections for two filters in a 3x2 formation. Subset creation happens here. 

        cutoutPath(str): /raid/scratch/data/COSMOS_test/random_cutouts/UVISTA_JH_DR6_15007186924767774_26432262759966996.fits"""

    from astropy.visualization import ZScaleInterval

    sciImgs = []
    catPaths = [] # list of paths to cats which include all detections
    catTables = [] # data - all dets
    subsetPaths = [] # list of paths to cats which include a subset of dets  i.e. MAG_APER[1] > 50
    subsetTables = [] # data - only MAG_APER[1] > 50 dets
    subsetNames = [] # name of subset cat
    bad_ratios = [] # all detections : bad dets
    
    ### Get matching coords
    cutoutDir = cutoutPath
    cutout_coords_file = os.path.join(cutoutDir, cutout_coords_file)
    detCat_coords_file = os.path.join(detPath, detCat_coords_file)
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

   #for f, filt in enumerate(filters):
        ## get science data
        # file with matching_coords in filename
    for f in os.listdir(cutoutPath):
        result = extract_filename(f)

        if result is None or len(result) != 6:
            continue  # Skip non-matching files

        prefix, filt, suffix,b,c, size_arcsec = extract_filename(f) # get prefix and suffix to rebuild filename
        filename =  prefix + filt + suffix + matching_coords +'_' + str(size_to_match) +'arcsec.fits' # rebuild filename with matched coords
        filePath = os.path.join(cutoutDir, filename)

        if '_wht_' in filePath:
            continue

        # ensure file exists
        if os.path.isfile(filePath):
            print(f"Found science image: {filePath}")
            sciImgs.append(filePath)
        elif os.path.isfile(filePath)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {filePath}")

    ## detections
    # get det with matching_coords in filename
    for f in detCats:
        detName = f
        if result is None or len(result) != 6:
            continue  # Skip non-matching files

        filt, bksz, bkfilt, ra, dec, size_arcsec = extract_filename(detName)
        detName = 'd'+filt+'bksz'+bksz+'bkfilt'+bkfilt+'_'+matching_coords+'_'+ str(size_arcsec)+'.fits'

        catPath = os.path.join(detPath, detName)
        # ensure file exists
        if os.path.isfile(catPath):
            print(f"Found det cat: {catPath} \n")
            catPaths.append(catPath)
            cat_table = Table.read(catPath) # a subset of detections
            catTables.append(cat_table)
        elif os.path.isfile(catPath)==False:
            print(f">>>>>>>>>>>>>>>>>>>>> MISSING file: {catPath}")

## make subset data
        # file with matching_coords in filename
        subset_file, bad_ratio = make_subset(f, detPath, '_bad_mags.fits', subsetCriterion=subsetCriterion)
        bad_ratios.append(bad_ratio)
        # ensure file exists
        subsetPath = os.path.join(detPath, subset_file)
        if os.path.isfile(subsetPath):
            print(f"Found subset catalogue: {subsetPath}")
            # TODO: check bad_mags file has matching coords
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

    for i, (sci_img, cat_table, subset_table) in enumerate(zip(sciImgs, catTables, subsetTables)):    
        sci_data = fits.getdata(sci_img) ## load sci data
        vmin, vmax = zscale.get_limits(sci_data)
        ax = axes[i, 0]
        ax.imshow(sci_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        plot_name = os.path.basename(sci_img)
        ax.set_title(f"{plot_name}", fontsize=10)
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
        ax_det.set_title(f"All detections: fix title")
        ax_det.invert_xaxis()  # invert RA axis to match sky convention
        ax_det.grid(True)

    ## plot subset
        ra = subset_table['ALPHA_J2000']
        dec = subset_table['DELTA_J2000']
        ax_sub = axes[i, 2]
        ax_sub.scatter(ra, dec, s=1, color='blue', alpha=0.5)
        ax_sub.set_xlabel("RA (deg)")
        ax_sub.set_ylabel("Dec (deg)")
        ax_sub.set_title(f"subset of detections: fix title",)
        ax_sub.text(np.max(ra), np.min(dec)-0.05, f'bad:all = {bad_ratios[i]:.5f}', fontsize=12, color='red')
        ax_sub.invert_xaxis()  # invert RA axis to match sky convention
        ax_sub.grid(True)

    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    plt.show()




