import os
import re
import warnings
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from datetime import datetime
from astropy.table import Table

# make cutours (optional)
# get_depths
## read_image_lis
# get required filter data, mask file
# image_depth
## define output Dirs
## check image and weight files exist
## get pixel scale from Header
## convert aperture sizes in arcsec to pix using pix scale

# makes different runs easier to find/distinguish
print("hu_depth_codes.py run at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

# base output dir
baseDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/'

######################## file-handling and set-up functions ################################
def get_fits_data(fitsImg, verbose=True):
    """ Returns data, header, and wcs of a fits image. """
    if verbose:
        print(f"Getting data for {fitsImg} \n")
    with fits.open(fitsImg) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

    return data, header, wcs

def read_image_lis(dirHere):
    """ Read .lis file and returns the contents as a Table. 

        dirHere - dir of .lis file"""
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print(f"Reading in images from {inputFile}")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
            
    return Table.read(inputFile, format = 'ascii.commented_header')


def open_cats(cats, open_subset=True, overwrite=True):
    """
        Opens catalogues for filters of interest in TOPCAT.

    INPUT(s)
        cats[list(str)]     list of paths to the cats you want to open. 
                            Typically, the ones you just made by looping through image_depths.  """
    command = f"topcat "
    catsToOpen = []

    for cat in cats:
        catDir = os.path.dirname(cat)
        catName = os.path.basename(cat)[:-len('.fits')]
        subsetBase = os.path.join(catDir, catName)
        if '_subset' in subsetBase: # if a subset already exists, skip
            continue
        else:
            catsToOpen.append(cat)
    files_str = " ".join(catsToOpen)
    files_str = " ".join(catsToOpen)
    command += f"{files_str} "
   
    if open_subset == False:
        command += " &"
        print(command)
        os.system(command)

    elif open_subset == True:
        # assuming subset = MAG_APER[1]>50
        subsets = []
        for cat in cats:
            # if subset name not already in cats, make subset 
            catDir = os.path.dirname(cat)
            catName = os.path.basename(cat)[:-len('.fits')]
            subsetBase = os.path.join(catDir, catName)

            if '_subset' in subsetBase: # if a subset already exists, skip
                continue
            else:
                subsetPath = subsetBase +'_subset.fits'
                subsets.append(subsetPath)
            if os.path.isfile(subsetPath) == False:
                table = Table.read(cat)
                subset = table[table["MAG_APER"][:, 1] > 50 ]  # criterion by which to create a subset
                subset.write(subsetPath, overwrite=overwrite)
                subsets.append(subsetPath)
        
        files_str = " ".join(subsets)
        command += f"{files_str} &"
        print(command)
        os.system(command)

        return subsets

######################## make cutout for testing on ################################
def make_cutout(filters, size_arcsec, centre_ra, centre_dec, dataDir, verbose=True, overwrite=True):
    """ Makes a square cutout of the original data files (image and weight).
        Used in testing scripts. Saves cutouts to dataDir+'cutouts/'
    
    INPUT(s)
        filters(list[str])  Used to find the image and weight file for the UVISTA filters of interest
        size_arcsec(int)    Desired size of the square cutout in arcsecs
        centre_ra(float)    RA coordinate of the centre of the the square cutout in degrees
        centre_dec(float)   As above, Declination.
   
    OUTPUT(s)
        cutoutImg(str)      Filename of the square image file cutout 
        cutoutWht(str)      Filename of the square weight file cutout 
    """

    ### set up warnings feedback
    ### make the cutout
    ## get paths of original images from which to make cutout
    ## retreive original data, header
    ## get pixel scale from original image
    ## define centre of square cutout in pix coords
    ## used Cutout2D to make cutout
    ## save cutout
    ## check header has been preserved

    import astropy.units as u
    from astropy.nddata import Cutout2D
    from astropy.coordinates import SkyCoord

    warnings_triggered = 0
    unexpectedValues = False
    unexpectedKeys = False

    # construct name(s) of original files for each filter
    origFilePaths = [] # files from which to cut out
    outFiles = [] # output files written by function

    for filt in filters:
        imgPath = dataDir + 'UVISTA_' + filt + '_DR6.fits'
        origFilePaths.append(imgPath)
        whtPath = dataDir + 'UVISTA_' + filt + '_DR6_wht.fits'
        origFilePaths.append(whtPath)
    
    for f, file_ in enumerate(origFilePaths):
        print(file_)
        data, header, wcs = get_fits_data(file_, verbose=verbose) # retrieve the original image data, header

        # Convert RA/Dec to pixel coordinates
        skycoord = SkyCoord(centre_ra, centre_dec, unit='deg')
        x, y = wcs.world_to_pixel(skycoord) ## define centre of square cutout in pix coords
        position = (x, y)

        # Get scale in arcsec/pixel
        try:
            scale = abs(header['CD1_1']) * 3600  # deg to arcsec
        except KeyError:
            scale = abs(header['CDELT1']) * 3600
        if verbose: 
            print(f"the pixel scale of {os.path.basename(file_)} is: ", scale, "arcsec/pix")
        size_pix = int(size_arcsec / scale)
        
        # Perform the cutout
        cutout = Cutout2D(data, position, (size_pix, size_pix), wcs=wcs)

        # Copy the original header and update it with the cutout WCS
        new_header = header.copy()
        new_header.update(cutout.wcs.to_header())

        hdu = fits.PrimaryHDU(data=cutout.data, header=new_header) # use cutout object
        hdul_out = fits.HDUList([hdu])
        
        ## save the cutout
        # make outputDir if necessary
        outputDir = dataDir + 'cutouts/'
        if os.path.isdir(outputDir) == False:
            os.system('mkdir ' + outputDir)
        savedCutouts = os.path.join(outputDir, 'cutoutNames.txt')

        # construct outputPath(s)
        ra_filename = str(centre_ra).replace('.', '')
        dec_filename = str(centre_dec).replace('.', '')
        outPrefix = dataDir + 'cutouts/' + 'UVISTA_' 
        outSuffix = '_DR6_' + ra_filename +'_'+ dec_filename +'_'+ 'size' + str(size_arcsec) + '.fits'
        cutoutPath = outPrefix + filt + outSuffix
        filt = os.path.basename(file_).split('_')[1] # get filt from filename so it matches to file_ in loop

        if  '_wht' in fits.open(file_).filename():
            outWhtSuffix = '_DR6_wht_' + ra_filename +'_'+ dec_filename +'_'+ 'size' + str(size_arcsec) + '.fits'
            cutoutWhtPath = outPrefix + filt + outWhtSuffix

            # Read existing cutouts (if file exists)
            if os.path.exists(savedCutouts):
                with open(savedCutouts, 'r') as f:
                    existing = {line.strip() for line in f}
            else:
                existing = set()

            # Append only if not already in file
            if cutoutWhtPath not in existing:
                with open(savedCutouts, 'a') as f:
                    print(f"Writing file: {cutoutWhtPath}")
                    f.write(cutoutWhtPath + '\n')
                    hdul_out.writeto(cutoutWhtPath, overwrite=overwrite)
                    outFiles.append(cutoutWhtPath) # for later in func
            cutouPathHeaderTest = cutoutWhtPath

        elif '_wht' not in fits.open(file_).filename():
            # Read existing cutouts (if file exists)
             if os.path.exists(savedCutouts):
                with open(savedCutouts, 'r') as f:
                    existing = {line.strip() for line in f}

             else:
                existing = set()
                
             if cutoutPath not in existing:
                with open(savedCutouts, 'a') as f:
                    print(f"Writing file: {cutoutPath}")
                    f.write(cutoutPath + '\n')
                    hdul_out.writeto(cutoutPath, overwrite=overwrite)
                    outFiles.append(cutoutPath) # for later in func
             cutouPathHeaderTest = cutoutPath

        print("Cutout file written. Performing Header checks... ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        # ensure header of data file is as expected 
        # i.e. if using a copy/cutout, header should == header of original file

        cutouthdulist = fits.open(cutouPathHeaderTest)
        cutoutHeader = cutouthdulist[0].header
        cutoutHeaderKeys = list(cutoutHeader.keys()) # are all the keys the same?

        originalhdulist = fits.open(file_)
        originalHeader = originalhdulist[0].header
        originalHeaderKeys = list(originalHeader.keys())

        # broard check: are all the keys the same?
        if cutoutHeaderKeys == originalHeaderKeys:
            print("ImageHeader matches expected format. Checking values...")

        elif cutoutHeaderKeys != originalHeaderKeys:
            # get a unique set of the element within the Keys lists
            print("what is it", type(cutoutHeaderKeys))
            cutoutHeaderKeys = set(cutoutHeaderKeys)
            originalHeaderKeys = set(originalHeaderKeys)
            # sort the sets alphabetically to allow comparison
            uniqueCutoutHeaderKeys = sorted(cutoutHeaderKeys - originalHeaderKeys)
            uniqueOrigImageHeaderKeys = sorted(originalHeaderKeys - cutoutHeaderKeys)

            warnings_triggered += 1
            unexpectedKeys = True
            if verbose:
                print("Keys unique to Cutout Header: ",uniqueCutoutHeaderKeys)
                print("Keys unique to Original Header: ",uniqueOrigImageHeaderKeys)

        # ensure values at each Key are as expected
        unmatchedValues = []
        skip_keys = {'NAXIS1', 'NAXIS2', 'CRPIX1', 'CRPIX2'} # these are expected to b different: skip them.
        for key in originalHeaderKeys:
            if key in cutoutHeaderKeys:
                if originalHeader[key] == cutoutHeader[key]:
                    if verbose:
                        print(f"original == cutout value: {originalHeader[key] == cutoutHeader[key]}", key, originalHeader[key], cutoutHeader[key])
                elif key in skip_keys:
                    continue
                else:
                    unexpectedValues = True
                    warnings_triggered += 1
                    unmatchedValues.append(key)
                    if verbose:
                        print(f"original == cutout value: {originalHeader[key] == cutoutHeader[key]}", key, originalHeader[key], cutoutHeader[key])

    print("\nCutout image and weight files have been made \n")
    if verbose:
        print(outFiles)
    if warnings_triggered == 0:
        print("All expected values are present.")
    if warnings_triggered > 0:
        print("\n>>>>> The following warnings were triggered during their creation: \n")
        if unexpectedKeys:
            warnings.warn(">>>>>>> WARNING: cutoutHeader does not match expected Header. Known potential issues: incorrect calculation of pixScale. Proceeding...")
            print("\nKeys unique to Cutout Header: ",uniqueCutoutHeaderKeys, '\n')
            print("Keys unique to Original Header: ",uniqueOrigImageHeaderKeys, '\n')
        if unexpectedValues:
            print("original == cutout value: False ... for these Keys: \n", unmatchedValues)
    
    return outFiles          
######################## do depths ##################################################
def aperture_photometry_blank(imageName, segMap, whtMap, apSize, gridSeparation = 100, pixScale = -99.0, next = 0, clean = False, outputFitsName = 'none', imageDir = '', verbose =False, field = 'NORMAL', overwrite = False):
    print("aperphotblank ", segMap, whtMap, apSize)
    #from photutils import CircularAperture
    #from photutils import aperture_photometry
    from astropy.io import fits
    from astropy.table import Table, hstack, join
    import time
    #import sep # aperture photometry from SExtractor!
    
    print("hellooooo",imageName, outputFitsName, (os.path.isfile(outputFitsName) and (overwrite == False)))
    # first check if output exists
    if os.path.isfile(outputFitsName) and (overwrite == False):

        origTable = Table.read(outputFitsName)
        cols = np.array(origTable.colnames)
        
        # check if all the columns are there
        change = -1
        
        # modify the column names
        for tt, typ in enumerate(['IMAGE', 'SEG', 'WHT']):
            for	ai, apD	in enumerate(apSize):
                
                oldcolname = '{0}_flux_{1}'.format(typ, ai)
                newcolname = '{0}_flux_{1:.1f}as'.format(typ, apD)
                
                if np.any(oldcolname == cols):
                    origTable.rename_column(oldcolname, newcolname)
                    change = 1
                    print('Renaming column from ', oldcolname, newcolname)
        
        # overwrite the file
        if change > 0:
            print(origTable.colnames)
            origTable.write(outputFitsName, overwrite = True)
            print('aperture phot file updated ', outputFitsName)

            
    # check if the column I want exists yet or not
    if os.path.isfile(outputFitsName) and (overwrite == False):
        
        origTable = Table.read(outputFitsName)
        cols = np.array(origTable.colnames)
        missingAps = np.ones(apSize.size, dtype = bool)
        
        for ai, apD in enumerate(apSize):
            reqCol = '{0}_flux_{1:.1f}as'.format('IMAGE', apD)
            kk = np.any(cols == reqCol)
            if kk:
                missingAps[ai] = False
                
        print('Checking the apertures ', missingAps)
        if np.any(missingAps):
            print('I need to re-run adding this aperture', apSize[missingAps])
        else:
            print('All required aps are present, do not need to run again')
            return

        append = True
        apSize = apSize[missingAps]
        
    else:
        append = False
        
    ## Get the pixel scale
    if verbose:
        print(imageName)
    
    hdulist = fits.open(imageName)
    header = hdulist[next].header
    imageData = hdulist[next].data
    
    if pixScale < 0.0:
        # read from header
        if 'CD1_1' in header:
            cdone_o = -3600.0*header['CD1_1']
        else:
            cdone_o = 3600.0*np.abs(header['CDELT1'])
        pixScale = round(cdone_o, 5)
        
    
    ## Get the apertures size, in pixels
    apSizePix = apSize/pixScale
    
    ## First just try a simple grid
    ## grab the dimensions
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    
    #gridSeparation = 20 ## pixels
    
    ## create arrays of the central coordinates
    numberX = int((naxis1-gridSeparation)/gridSeparation)
    numberY = int((naxis2-gridSeparation)/gridSeparation)
    numberApertures = numberX*numberY
    print('The number of apertures is ', numberApertures)
    
    
    #apertureArray = np.zeros([numberApertures, 2])
    xArray = np.zeros(numberApertures)
    yArray = np.zeros(numberApertures)
    halfGrid = gridSeparation/2
        
    numberHere = 0
    for xi in range(numberX):
        for yi in range(numberY):
            #apertureArray[numberHere, :] = [halfGrid+xi*gridSeparation, halfGrid+yi*gridSeparation]
            # print "the coords are ", apertureArray[:, numberHere]
            xArray[numberHere] = halfGrid+xi*gridSeparation
            yArray[numberHere] = halfGrid+yi*gridSeparation
            numberHere = numberHere + 1

    ## now do aperture photometry on both
    ## setup the apertures
    radii = apSizePix/2.0
    if verbose:
        print("Here")

    ## 1) the image
    tic = time.time()
    phot_image = aperture_phot_fast(imageData, xArray, yArray, radii)
    toc = time.time()
    hdulist.close()
    if verbose:
        print("Finished doing the photometry for image in time {0}".format(toc-tic))

    ## 2) the seg
    ## I don't care about interpolation here
    hdulist = fits.open(segMap)
    segData = hdulist[next].data
    phot_seg = aperture_phot_fast(segData, xArray, yArray, np.array(radii), subpix = 1)    
    hdulist.close()
    if verbose:
        print("Finished doing the photometry for seg in time {0}".format(toc-tic))
    
    ## 3) the wht
    ## to exclude pixels off the edge
    if whtMap[-4:].lower() == 'none':
        print("No weight data. ")
        ## Just use the image instead.
        phot_wht = Table(phot_image, copy = True)
        
        ## absolute these
        for ri, r in enumerate(radii):
            name = 'flux_' + str(ri)
            phot_wht[name] = np.abs(phot_wht[name])
        
    else:
        hdulist = fits.open(whtMap)
        whtData = hdulist[next].data
        phot_wht = aperture_phot_fast(whtData, xArray, yArray, np.array(radii), subpix = 1)
        print("phot_wht", phot_wht, "that means sep has run")
    # centre means a pixel is either in or outside the aperture
        hdulist.close()
        
    ## Save these results to a fits file
    ## I can do cuts etc in another code
    ## to speed this up!!

    if outputFitsName == 'none':
        directory = imageDir + 'depths/catalogues/'
        filterName = segMap[0:segMap.rfind('_')]
        outputFitsName = filterName + + '_aperPhot.fits'
        
    # stack the tables
    bigTable = hstack([phot_image, phot_seg, phot_wht], table_names=['IMAGE', 'SEG', 'WHT'], uniq_col_name='{table_name}_{col_name}')
    
    bigTable = Table(bigTable)
    
    # remove columns to make it more streamlined
    bigTable.remove_column('WHT_xcenter')
    bigTable.remove_column('WHT_ycenter')
    bigTable.remove_column('SEG_xcenter')
    bigTable.remove_column('SEG_ycenter')

    if append:
        for ri, r in enumerate(radii):
            bigTable.remove_column('SEG_flux_{0}'.format(ri))
            bigTable.remove_column('WHT_flux_{0}'.format(ri))
    else:
        for ri, r in enumerate(radii):
            if ri != 2:
                bigTable.remove_column('SEG_flux_{0}'.format(ri))
                bigTable.remove_column('WHT_flux_{0}'.format(ri))        

                
        if clean:

            print('Cleaning at the aperture phot level')
            # remove the bad columns here.
            smallNum = 0.0000001
            deltaZero = 1E-13 #0.00000001
            apString = '2'
            
            seg_sum = np.array(bigTable['SEG_flux_' + apString])
            wht_sum = np.array(bigTable['WHT_flux_' + apString])
            ap_sum = np.array(bigTable['IMAGE_flux_' + apString])
            
            if field == 'NIRSPEC':
                good_indicies = (seg_sum < 0.5)  & (wht_sum > -1E-28) & (wht_sum < 1E28)
                
            else:
                good_indicies = (seg_sum < 0.5)  & (wht_sum > smallNum) & ((ap_sum > deltaZero) | (ap_sum < -deltaZero))
            #good_indicies = np.where(good_indicies)
            bigTable = bigTable[good_indicies]

    # rename the columns with the diameter size
    for tt, typ in enumerate(['IMAGE', 'SEG', 'WHT']):
        for ai, apD in enumerate(apSize):
            
            oldcolname = '{0}_flux_{1}'.format(typ, ai)
            newcolname = '{0}_flux_{1:.1f}as'.format(typ, apD)
            if oldcolname in np.array(bigTable.colnames):
                print('Renaming column from ', oldcolname, newcolname)
                bigTable.rename_column(oldcolname, newcolname)
            
    bigTable.info
    
    # remove the wht and seg columns
    #bigTable.remove_column('WHT_flux_' + apString)
    #bigTable.remove_column('SEG_flux_' + apString)

    if append:
        print(bigTable.colnames)
        print(origTable.colnames)
        
        # join with the big table!
        print('Appending to big aper phot table, lengths = {0}, {1}'.format(len(bigTable), len(origTable)))
        bigTable = join(origTable, bigTable, keys = ['IMAGE_xcenter', 'IMAGE_ycenter'])
        print(bigTable.colnames)
        print('After ', bigTable)
    #exit()
        
    #bigTable.meta['aperture_photometry_args'] = ''
    bigTable.write(outputFitsName, overwrite = True)
    print("tired",bigTable)
    print("Aperture table has been saved to ", outputFitsName)
     
    return


def image_depth(imagePath, zeropoint, cutouts=[], size='none', apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtPath='NONE', whtType='NONE', IRACapDiametersAS=np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segPath='NONE', outputDir='none', filterName='NONE', numApertures=300, step=200, overwrite=False, inputSex=baseDir+'data/bertin_config/video_mine.sex', strips=False, bgSub=True, mask='none', gridSepAS=3.0):

    warnings_triggered = 0
    cd1Test = False

    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    
    # Tinsley SE env
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/share/sextractor'

    # define the output files
    if outputDir == 'none':
        outputDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/'
        if os.path.isdir(outputDir) == False:
            os.system('mkdir ' + outputDir)
        print("outputDir will be: ", outputDir)

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
   
    if len(cutouts) > 0:
        Dir = imagePath.split('UVISTA_')[0]
        cutoutDir = Dir + 'cutouts/'

        if '.fits' not in imagePath.split('/')[-1]:
            imagePath = cutoutDir + imagePath.split('/')[-1]+'.fits'
        else:
            imagePath = cutoutDir + imagePath.split('/')[-1]

        print("the imagepath is", imagePath)

        whtDir = whtPath.split('UVISTA_')[0]
        whtDir = whtDir + 'cutouts/'

        if '.fits' not in whtPath.split('/')[-1]:
            whtPath = whtDir + whtPath.split('/')[-1]+'.fits'
        else:
            whtPath = whtDir + whtPath.split('/')[-1]

        print("the whtpath is", whtPath)

    parts = imagePath.split('/') # split image dir by /
    if len(parts) > 1:
        baseName = parts[-1] # image filename string

    # also remove the file extension
    pparts = baseName.split('.')
    baseName = pparts[0] # image name string without file extension
    print("The base name is ", baseName)

    # check all the necessary files exist
    imyes = os.path.isfile(imagePath)
    if imyes == False:
        print("ERROR: the image file does not exist, check the path:\n" + imagePath)
        exit()

    if whtType != "NONE":
        whtyes = os.path.isfile(whtPath) # find correct weight file for image
        if whtyes == False:
            print("ERROR: the weight file does not exist, check the path:\n" + whtPath)            
            exit()

    ###### seg map, SE, and Bkg Subtraction ###########################################

    if filterName == 'NONE':
        filterName = baseName # image name string without file extension

    # get the pixel scale
    hdulist = fits.open(imagePath)
    imageHeader = hdulist[0].header
    imageHeaderKeys = list(imageHeader.keys())

    # check CD1_1 is in Header - if using a new cutout, this has already been checked but good to check again
    cd1Test = ['CD1_1' in imageHeaderKeys]
    print(f">>>> WARNING: Passing file to SE. Header includes CD1_1?: {cd1Test} Is this the right file?", '\n', imagePath, '\n') 
    if cd1Test == True:
        warnings_triggered += 1

    if 'CD1_1' in imageHeader:
        cdone_o = -3600.0*imageHeader['CD1_1']
    else:
        cdone_o = 3600.0*np.abs(imageHeader['CDELT1'])
    pixScale = round(cdone_o, 5)
    print("The pixel scale is {0:.4f}".format(pixScale))
    
    if segPath == 'NONE':
        # run source extractor
        ## convert aperture sizes in arcsec to pix using pix scale
        ## SE needs aperture sizes in pix (PHOT_AP)
        apStringPix = str(apDiametersAS[0]/pixScale) 
        for ai in range(1, apDiametersAS.size):
            apStringPix = apStringPix + ',' + str(apDiametersAS[ai]/pixScale)
            
        ### now run SE 
        inputSex = baseDir + 'data/bertin_config/video_mine.sex'
        param_file = baseDir + 'data/bertin_config/default.param'

        if len(cutouts) > 0: 
            segPath = imageDir + filterName + size +'_cutout_seg.fits'  
            bgSubPath = imageDir  + filterName + size + '_cutout_bgsub.fits'
            bgMapPath = imageDir  + filterName + size + '_cutout_bgmap.fits'
            outputCatalogue = catDir + 'd' + size + filterName + '_cutout.fits'
        else:
            segPath = imageDir + filterName + '_seg.fits'  
            bgSubPath = imageDir  + filterName + '_bgsub.fits'
            bgMapPath = imageDir  + filterName + '_bgmap.fits'
            outputCatalogue = catDir + 'd' + filterName + '.fits'

        keywordsbase = ' -CATALOG_TYPE FITS_1.0 -CATALOG_NAME ' + \
                       outputCatalogue + \
                       ' -MAG_ZEROPOINT '+ str(zeropoint) + \
                       ' -WEIGHT_TYPE ' + whtType + \
                       ' -WEIGHT_IMAGE ' + whtPath
        # -CHECKIMAGE_TYPE "BACKGROUND" for background map -BACKGROUND for bkg subtracted
        keywords = keywordsbase + \
                   ' -CHECKIMAGE_TYPE "BACKGROUND,-BACKGROUND,SEGMENTATION" '\
                   +'-CHECKIMAGE_NAME "' + \
                   bgMapPath + ',' + bgSubPath + ',' + segPath + '" -PHOT_APERTURES ' \
                   + apStringPix # "', '" needed to interpret bgMApPAth etc as str

# with -CATALOG_TYPE FITS_1.0, the input parameter settings are saved in the header.
# NOTE: the fits option can't handle array outputinformation such as MAG APER, FLUX RADIUS if more than one value! - SE for Dummies 
# PHOT_APERTURES specifies the number and size of apertures (in pix)

        command = '/usr/local/sextractor/bin/sex '+ imagePath + ' -c ' + inputSex + ' -PARAMETERS_NAME ' + param_file + ' ' + keywords
        
        if os.path.isfile(bgSubPath) == False or os.path.isfile(segPath) == False or overwrite:
            print("The SEG and BG subtracted map do not exist, or overwrite:True.  Running SE like so: \n")

            print(command)
            os.system(command)
        else:
            print(f"The SEG and/or BG subtracted map exist at: {segPath, bgSubPath} \n")

    ### TODO: get the rest of this function from new_depth_codes i.e. Next step is to place apertures down
    aperPhotFile = aperDir + filterName + '_aperPhot.fits'
    overwrite = False
    if os.path.isfile(aperPhotFile) == False or overwrite == True:
        # define the grid separation
        gridSepPixels = np.ceil(gridSepAS/pixScale) # 5'' separation
        print("test: ",gridSepPixels)
        # gridSepPixels = 10.0
        
        # make this tunable...
        if bgSub == False:
            bgSubPath = imagePath
            print("Not using bg subtracted image.")
        print("Measuring the aperture photometry.")

        ii = segPath.rfind('NIRSPEC')
        if ii > 0:
            field = 'NIRSPEC'
        else:
            field = 'NONE'

        print("aperture_photometry_blank is running")
        # if aperphotfile doesn't exist, call ap phot blank - the output fits is aperphotfile
       # aperture_photometry_blank(bgSubPath, segPath, whtPath, apDiametersAS, gridSeparation = gridSepPixels, clean = True, outputFitsName = aperPhotFile, imageDir = imageDir, field = field, overwrite = overwrite)

    return

def get_depths(fieldName, cutouts, size='none',queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir=baseDir+'data/', outputDir='none', overwrite=False):

    # set the grid seperation in arcsec
    if fieldName == 'NIRSPEC':
        print('Nirspec data, setting grid separation to be much lower!')
        gridSepAS = 0.5
    else:
        gridSepAS = 3.0
        
    # Read in the images file
    dirHere = dataDir + fieldName + '/'
    imagedata = read_image_lis(dirHere)
    availableFilters = np.array(imagedata['Name'])
    print("The available filters are ", availableFilters)

    if reqFilters[0] != 'all':

        # reduce the list down!
        keep = np.zeros(availableFilters.size, dtype=bool)
        for bf, keepfilt in enumerate(reqFilters):
            ii = np.where((availableFilters == keepfilt))
            keep[ii] = True
            
        imagedata = imagedata[keep]
        availableFilters = imagedata['Name']
    print("Only running for filters ", availableFilters)

    # check if we are in UltraVISTA
    stripFilt = np.array(['Y', 'J', 'H', 'Ks', 'JH', 'YJ', 'HKs', 'Y_DR4', 'J_DR4', 'H_DR4', 'Ks_DR4', 'NB118', 'NB118_DR4'])

    # loop through each filter
    # for the queue run, run each as a separate file...

    for fi, filterName in enumerate(availableFilters):
        # define the images etc to send through
        imageName = imagedata['Image'][fi]
        whtName = imagedata['Weight'][fi]
        whtType = imagedata['Wht_type'][fi]
        zeropoint = imagedata['zeropoint'][fi]
        imageDir = imagedata['directory'][fi]
        maskName = '/raid/scratch/data/masks/{0}/{1}'.format(fieldName, imagedata['mask'][fi])
        
        if len(cutouts) > 0:    
            for cutout in set(cutouts): # for unique items in list
                cutoutPath = cutout
                cutoutName = os.path.basename(cutoutPath)
                imageName = imageName.split('.')[0]
                if cutoutName.startswith(imageName) and ('wht' not in cutoutName):
                    imageName = cutoutName
                    print("The imageName is ", 'cutouts/'+imageName)

                if 'wht' in cutout:
                    cutWhtPath = cutout
                    cutWhtName = os.path.basename(cutWhtPath)
                    whtName = whtName.split('.')[0]
                if cutoutName.startswith(whtName):
                    whtName = cutoutName
                    print("The whtName is ", 'cutouts/'+whtName)

        if imageDir == 'here':
            imageDir = dataDir + fieldName + '/'
        
        # Now spawn the depths!
        if queue == 'none':
            strips = True
            if fieldName == 'COSMOS':
                jj = (filterName == stripFilt)
                if np.any(jj):
                    strips = True 
            image_depth(imageDir+imageName, zeropoint, cutouts=cutouts, size=size, whtPath=imageDir+whtName, whtType=whtType, outputDir=outputDir, strips=strips, filterName=filterName, overwrite=overwrite, mask=maskName, gridSepAS=gridSepAS, apDiametersAS=apDiametersAS)

        else:
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

    # after each filter has been SE'd, open the cats
    cats = []
    if outputDir == 'none':
        outputDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/catalogues/'
        for cat in os.listdir(outputDir):
            for filt in reqFilters:
                pattern = (rf"d{size}{filt}_cutout\.fits")
                if re.fullmatch(pattern, cat):
                    catPath = os.path.join(outputDir, cat)
                    cats.append(catPath)

    subsets = open_cats(cats, open_subset=True, overwrite=True)

    for subset in subsets:
        for filt in reqFilters:
            subpattern = (rf"{baseDir}depths/catalogues/"
                        rf"d{size}{filt}_cutout_subset\.fits")
            if re.fullmatch(subpattern, subset):
                print("test!!!!!", subset)
                cutoutPath = f"{baseDir}depths/catalogues/d{size}{filt}_cutout.fits"
                cutoutTable = Table.read(cutoutPath)
                subsetTable = Table.read(subset)
                print(f"{filt}-bad-proportion: ", round(100*(len(subsetTable)/len(cutoutTable)),5), "%")
    return
    

