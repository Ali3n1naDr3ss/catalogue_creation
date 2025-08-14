

"""
17/07/25
hu_depths_codes.py
Testing new_depths_codes.py - is SE running correctly? Why aren't all the objects recorded with a measured magnitude?

"""

import numpy as np
from astropy.table import Table
import os
import matplotlib.pyplot as plt


baseDir = '/raid/scratch/'

##############################################################################################################

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
    
    plot_detections(outputCatalogue) # plot all detections
    outputPath = make_subset(outputCatalogue, baseDir+'depths/COSMOS_test/catalogues/')
    plot_detections(outputPath)
    plot_bkg_model([[32,9],[1,2]])
    #unhide_hidden_files()
    return

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
                    
            image_depth(imageDir + imageName, zeropoint, whtName = imageDir + whtName, whtType = whtType, outputDir = outputDir, strips = strips, filterName = filterName, overwrite = overwrite, mask = maskName, gridSepAS = gridSepAS, apDiametersAS = apDiametersAS)

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
    
    return

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


def make_subset(big_cat, outputDir):

    from astropy.table import Table
    filename_root = os.path.splitext(os.path.basename(big_cat))[0]
    subset_cat_name = filename_root + '_bad_mags.fits'
    os.makedirs(filename_root, exist_ok=True) # create dir if it doesn't exist
    table = Table.read(big_cat)
    print(f"Objects in unfiltered catalogue: {len(table)}")
    filtered = table[table["MAG_APER"][:, 1] > 50]  # criteria by which to create a subset
    outputPath = os.path.join(outputDir, subset_cat_name)  
    filtered.write(outputPath, overwrite=True)
    print(f"Objects in filtered sub-catalogue: {len(filtered)}")
    print(f"Saved to {outputPath}")
    return outputPath

def unhide_hidden_files(directory):
    for fname in os.listdir(directory):
        if fname.startswith('.') and not fname.startswith('..'):
            new_name = fname.lstrip('.')
            os.rename(os.path.join(directory, fname), os.path.join(directory, new_name))
            print(f"Renamed hidden file {fname} to {new_name}")

def plot_bkg_model(param_combos, baseDir=baseDir, testDir='depths/COSMOS_test/'):
    """
    what affect do each set of SE bkg modelling parameters have?
    Do visual check with this side-by-side plotter

    param_combos(list): list of lists in the format
                        [[back_size, back_filtersize], [[back_size1, back_filtersize1]]]
    """
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval
    import re
    from collections import defaultdict


    testImgDir = baseDir + testDir + 'images/'
    origImgDir = baseDir + 'data/COSMOS/'

    # get all filenames in Test dir and verify test parameters are present
    all_files = [f for f in os.listdir(testImgDir) if os.path.isfile(os.path.join(testImgDir, f))]
    test_imgs = []    
    for f in all_files:
        for param_combo in param_combos:
            if str(param_combo[0]) and str(param_combo[1]) in f:
                #print(f"Test file found: {f}")
                test_imgs.append(testImgDir+f)
            else:
                 print(f"Test file NOT found: {f}")

    # Compile regex
    pattern = re.compile(r"^(?P<band>.+?)bk_sz(?P<bksz>.+?)bkfisz(?P<bkfsz>.+?)_(?P<type>.+?)\.fits")

    # Group files by (bksz, bkfsz)
    grouped = {}
    bands = []
    for fname in test_imgs:
        match = pattern.match(os.path.basename(fname)) # extracts info from filename if name matches supplied pattern
        if not match:
            print(f"Skipping unmatched filename: {fname}")
            continue  # Skip files that don?t match the pattern

        band = match.group("band")
        bands.append(band)
        bksz = match.group("bksz")
        bkfsz = match.group("bkfsz")
        ftype = match.group("type")

        key = (bksz, bkfsz) # combination of test params

        if key not in grouped:
            grouped[key] = {"J": {}, "JH": {}} # add keys for subdicts
            print(f"Adding: grouped[{key}][{band}][{ftype}]")
        grouped[key][band][ftype] = os.path.join(testImgDir, fname)

    # get all filenames in Orig dir and verify test bands are present
    all_files = [f for f in os.listdir(origImgDir) if os.path.isfile(os.path.join(origImgDir, f))]
    orig_imgs = []   
    orig_pattern = re.compile(r"^UVISTA_(?P<band>.+?)\_DR6.fits") 

    for f in all_files:
        basename = os.path.basename(f)
        match = orig_pattern.match(basename)
        if match:
            matched_band = match.group("band")
            #print(f"Original file found: {f}")
            if matched_band in ['J', 'JH']:
                print(f"Using science image: {f}")
                orig_imgs.append(os.path.join(origImgDir, f))

    # display bgmap, bgsub, original for both filters

    # Find first group with matching files for J and JH
    selected_group = None
    for (bksz, bkfsz), bands in grouped.items():
        if all(ftype in bands["J"] for ftype in ["bgmap", "bgsub"]) and \
           all(ftype in bands["JH"] for ftype in ["bgmap", "bgsub"]):
            selected_group = (bksz, bkfsz)
            break

    if not selected_group:
        raise ValueError("No group found with complete J and JH bgmap, and bgsub images.")

        bksz, bkfsz = selected_group
    # --- Load bgmap and bgsub images ---
    j_images = []
    jh_images = []
    ftypes = ["bgmap", "bgsub"]

    for ftype in ftypes:
        j_path = grouped[(bksz, bkfsz)]["J"][ftype]
        jh_path = grouped[(bksz, bkfsz)]["JH"][ftype]

        j_data = fits.getdata(j_path)
        jh_data = fits.getdata(jh_path)

        j_images.append((j_data, ftype))
        jh_images.append((jh_data, ftype))

    # --- Load original images ---
    orig_pattern = re.compile(r"^UVISTA_(?P<band>.+?)_DR6\.fits")

    for img_path in orig_imgs:
        basename = os.path.basename(img_path)
        match = orig_pattern.match(basename)
        if not match:
            continue

        band = match.group("band")
        data = fits.getdata(img_path)

        if band == "J":
            j_images.append((data, "orig"))
        elif band == "JH":
            jh_images.append((data, "orig"))

    # --- Plot ---
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval()

    # Top row: J-band
    for i, (data, ftype) in enumerate(j_images):
        vmin, vmax = zscale.get_limits(data)
        ax = axes[0, i]
        ax.imshow(data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"J bk_sz{bksz} bkfisz{bkfsz} {ftype}", fontsize=10)
        ax.axis('off')

    # Bottom row: JH-band
    for i, (data, ftype) in enumerate(jh_images):
        vmin, vmax = zscale.get_limits(data)
        ax = axes[1, i]
        ax.imshow(data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"JH bk_sz{bksz} bkfisz{bkfsz} {ftype}", fontsize=10)
        ax.axis('off')

    plt.tight_layout()
    plt.show()
                    





