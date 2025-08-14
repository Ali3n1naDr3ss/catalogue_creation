

"""
17/07/25
hu_depths_codes.py
Testing new_depths_codes.py - is SE running correctly? Why aren't all the objects recorded with a measured magnitude?
why are 100% of objects measured with bad magnitudes!?

"""

import numpy as np
from astropy.table import Table
import os

baseDir = '/raid/scratch/'

##############################################################################################################

def bkg_plotter(param_combos, depthDir='/raid/scratch/depths/COSMOS_test/', dataDir='/raid/scratch/data/COSMOS_test', verbose=True):
    """
    what affect do each set of SE bkg modelling parameters have?
    Do visual check with this side-by-side plotter

    param_combos(list): SE parameters used in bkg modelling test: [BACK_SIZE, BACK_FILTERSIZE]
    depthDir(str):      directory of SE'd find_depth catalogues
    dataDir(str):       directory of science images (full-size and cutouts)
    """

    import re
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from collections import defaultdict
    from astropy.visualization import ZScaleInterval

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
        back_size = 100
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
            print(f"The SEG and/or BG subtracted map exist at: {segName, bgSubName} \n  Moving on...")

    # next step after plotting 17/07/25 ###########################################################
    """
    # Next step is to place apertures down
    aperPhotFile = aperDir + filterName + '_aperPhot.fits'
    overwrite = False
    if os.path.isfile(aperPhotFile) == False or overwrite == True:
        # define the grid separation
        gridSepPixels = np.ceil(gridSepAS/pixScale) # 5'' separation
    #        gridSepPixels = 10.0
        
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
    print("aperPhotFile", aperPhotFile, bgSubName)
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
    print("Output file saved to ", depthFile) """
    # next step after plotting 17/07/25 ##########################################################
    bkg_plotter([32,9])
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
    print("The available filters are:\n", availableFilters)

    if reqFilters[0] != 'all':

        # reduce the list of available filters
        keep = np.zeros(availableFilters.size, dtype=bool)
        for bf, keepfilt in enumerate(reqFilters):
            ii = np.where((availableFilters == keepfilt))
            keep[ii] = True
            
        imagedata = imagedata[keep]
        availableFilters = imagedata['Name']
        
    print("Only running for filters: ",'\n', availableFilters)

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
            print("No queue. Running now. \n")

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
            

            print(apDiametersASstring)
            print("Spawning in the queue...", queue)
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
    
