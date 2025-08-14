import os
import warnings
import numpy as np
from astropy.table import Table

baseDir = '/raid/scratch/'

def image_depth(imageName, zeropoint, apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtName = 'NONE', whtType = 'NONE', IRACapDiametersAS = np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segName = 'NONE', outputDir = 'none', filterName = 'NONE', numApertures = 300, step = 200, overwrite = False, inputSex = baseDir + 'data/bertin_config/video_mine.sex', strips = False, bgSub = True, mask = 'none', gridSepAS = 3.0):

    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    
    from astropy.io import fits
    #os.environ['EXTRACTOR_DIR'] = '/usr/local/shared/sextractor/2.25.0/share/sextractor'
    # updated for tinsley
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/share/sextractor'
    
    #    bgSub = False
    
    if filterName[0:2] == 'ch':
        print('Using IRAC apertures.')
        apDiametersAS = IRACapDiametersAS

    # main code where all the action happens!
    # first define the output files
    if outputDir == 'none':
        # put everything here!
        # make a new directory
        outputDir = 'depths/'
        if os.path.isdir(outputDir) == False:
            os.system('mkdir ' + outputDir)

    # for the seg map
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
        
    # make a sub-directory here?  no but useful
    parts = imageName.split('/')
    if len(parts) > 1:
        baseName = parts[-1]

    # also remove the cat bit
    pparts = baseName.split('.')
    baseName = pparts[0]
    print("The base name is ", baseName)

    # make a directory for the output
    #specificDir = outputDir + baseName + '/'
    #if os.path.isdir(specificDir) == False:
    #    os.system('mkdir ' + specificDir)
        
    # check all the necessary files exist
    imyes = os.path.isfile(imageName)
    if imyes == False:
        print("ERROR: the image file does not exist, check the path:\n" + imageName)
        exit()

    if whtType != "NONE":
        whtyes = os.path.isfile(whtName)
        if whtyes == False:
            print("ERROR: the weight file does not exist, check the path:\n" + whtName)            
            exit()
            
    ####################################################
    # Get the seg map

    if filterName == 'NONE':
        filterName = baseName

    # get the pixel scale
    hdulist = fits.open(imageName)
    imageHeader = hdulist[0].header
    imageHeaderKeys = list(imageHeader.keys())
    # ensure header of data file is as expected 
    # i.e. if using a copy/cutout, header should == header of original file
    #originalImageName = os.path.basename(imageName)[:-len('.fits')]+'_backup.fits'
    originalImageName = os.path.basename(imageName)
    originalHDUlist = fits.open(os.path.dirname(imageName)+'/'+originalImageName)
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
            
        # now run sex
        #HOLLY 29/05/25 - improve bkg sub for VISTA
        # now run sex
        segName = imageDir + filterName + '_seg.fits'
        bgSubName = imageDir  + filterName + '_bgsub.fits'
        outputCatalogue = catDir + 'd' + filterName + '.fits'

        keywordsbase = ' -CATALOG_TYPE FITS_1.0 -CATALOG_NAME ' + \
                       outputCatalogue + \
                       ' -MAG_ZEROPOINT '+ str(zeropoint) + \
                       ' -WEIGHT_TYPE ' + whtType + \
                       ' -WEIGHT_IMAGE ' + whtName
        
        keywords = keywordsbase + \
                   ' -CHECKIMAGE_TYPE "-BACKGROUND,SEGMENTATION" '\
                   +'-CHECKIMAGE_NAME "' + \
                   bgSubName + ',' + segName + '" -PHOT_APERTURES ' \
                   + apStringPix 
    
        command = '/usr/local/sextractor/bin/sex '+ imageName +' -c ' + inputSex + keywords
        print(command)
        
        if os.path.isfile(bgSubName) == False or os.path.isfile(segName) == False or overwrite:
            print("The SEG and BG subtracted map do not exist.  Running...")
            os.system(command)
        else:
            print("The SEG and BG subtracted map exist.  Moving on...")

 
    #######################################################################
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
    print("Output file saved to ", depthFile)
    
    return



# Note from Rohan (01/02/2022): I am unable to run with dataDir='/users/bowlerr/vardygroupshare/data', change to absolute path '/mnt/vardy/vardygroupshare/data'.
def get_depths(fieldName, queue = 'none', reqFilters= ['all'], apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir = baseDir + 'data/', outputDir = 'none', overwrite = False):

    
    # set the grid sep
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

    # check if we are in UltraVISTa
    stripFilt = np.array(['Y', 'J', 'H', 'Ks', 'JH', 'YJ', 'HKs', 'Y_DR4', 'J_DR4', 'H_DR4', 'Ks_DR4', 'NB118', 'NB118_DR4'])

    # loop through each filter
    # for the queue run, run each as a separate file...

    for fi, filterName in enumerate(availableFilters):
        # define the images etc to send through
        imageName = imagedata['Image'][fi]
        imageName = imageName.split('.')[0]+'_backup.fits'
        whtName = imagedata['Weight'][fi]
        whtName = whtName.split('.')[0]+'_backup.fits'
        whtType = imagedata['Wht_type'][fi]
        zeropoint = imagedata['zeropoint'][fi]
        imageDir = imagedata['directory'][fi]
        maskName = '/raid/scratch/data/masks/{0}/{1}'.format(fieldName, imagedata['mask'][fi])
        
        
        if imageDir == 'here':
            imageDir = dataDir + fieldName + '/'
        
        # Now spawn the depths!
        if queue == 'none':
            print("Running here ")

            strips = True
            if fieldName == 'COSMOS':
                jj = (filterName == stripFilt)
                if np.any(jj):
                    strips = True
                    
            image_depth(imageDir + imageName, zeropoint, whtName = imageDir + whtName, whtType = whtType, outputDir = outputDir, strips = strips, filterName = filterName, overwrite = overwrite, mask = maskName, gridSepAS = gridSepAS, apDiametersAS = apDiametersAS)


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
    
    return
    
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

def read_image_lis(dirHere):

    # read in filters
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
            
    return Table.read(inputFile, format = 'ascii.commented_header')

def aperture_phot_fast(imageData, xArray, yArray, radii, subpix = 5):
    
    # Subpix = 5 is what SEXtractor uses.
    import sep
    from astropy.table import Table, Column
    
    data = imageData.byteswap().newbyteorder()    
    for ri, r in enumerate(radii):
        #print(data, r, subpix)
        #print(xArray)
        flux, fluxerr, flag = sep.sum_circle(data, xArray, yArray, r, subpix = subpix)
        
        if ri < 1:
            # create a table of the results
            phot_apertures = Table([xArray, yArray, flux], names = ['xcenter', 'ycenter', 'flux_0'], dtype = ['f4', 'f4', 'f4'])
            
        else:
            newcolumn = Column(name = 'flux_' + str(ri), data = flux, dtype = 'f4')
            phot_apertures.add_column(newcolumn)
    
    return phot_apertures

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
    #magmin = 23
    #magmax = 28
    # remove things close to zero!
    deltaZero = 1E-13
    
    # extract name for local depth results
    #if clean:
    #    basename = inputTableFile[:-21]
    #else:
    basename = inputTableFile[:-13]

    colourArray = ['Blue', 'Green', 'Green', 'Green', 'Green', 'Red', 'Red', 'Red', 'Red', 'Red']
    
    ## Now loop through the different regions of the image
    if strips:
        # for ultravista!
        regions = ['fullimage', 'stripone', 'striptwo', 'stripthree', 'stripfour', 'gap1', 'gap2', 'gap3', 'gap4']
        regions = ['full', 'str1', 'str2', 'str3', 'str4', 'gap1', 'gap2', 'gap3', 'gap4']
        #deepStrips_ra_low = [149.3, 149.65, 150.02, 150.4]
        #deepStrips_ra_high = [149.5, 149.85, 150.25, 150.6]
        # convert to x and y!
        #strips_x_high = [9513, 18626, 27500, 35895]
        #strips_x_low = [4716, 13111, 22703, 31098]
        #x_low_limit = strips_x_low + [0] + strips_x_high
        #x_high_limit = strips_x_high + strips_x_low + [36261]
        
        ## this is conservative, gives only the deepest part
        ## for the gaps I should probably add a few 2000...
        #gaps = gapone | gaptwo | gapthree| gapfour | gapfive
        
    else:
        regions = ['full']
        
    global_depth = np.zeros([len(regions), apDiametersAS.size])
    median_local_depth = np.zeros([len(regions), apDiametersAS.size])
    mode_local_depth = np.zeros([len(regions), apDiametersAS.size])
    
    # define a nice figure
    # save the plot somewhere sensible
    # extract the filtername
    startI = inputTableFile.rfind('/') + 1
    endI = inputTableFile.rfind('_')
    filterName = inputTableFile[startI:endI]
    plotName = plotDir + filterName + '_' + str(numApertures) + '_{0}.pdf'.format(step)
    #plotName = "test.pdf"
    pdf = matplotlib.backends.backend_pdf.PdfPages(plotName)
    
    # Loop through the available apertures
    for ai, apDiAS in enumerate(apDiametersAS):
        
        print("Extracting local and global depths for aperture = ", apDiAS, " as.")
        
        # The apertures are strings _0, _1 etc
        #apString = '_' + str(ai)
        apString = '_{0:.1f}as'.format(apDiAS)
        
        #if clean:
        localDepthsFile = basename + str(apDiAS) + 'as_gridDepths_{0}_{1}.fits'.format(numApertures, step)
        #else:
        #    localDepthsFile = basename + str(apDiAS) + 'as_gridDepths.fits'
            
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
        #print inputTable.colnames
        #print inputTable['IMAGE_xcenter']
        #print inputTable['IMAGE_ycenter']
         ## cut the table depending on apertures
         ## only accept coordinates where the seg map
         ## aperture is blank, and the wht map is > 0
            seg_sum = np.array(inputTable['SEG_flux' + apString])
            wht_sum = np.array(inputTable['WHT_flux' + apString])
            ap_sum = np.array(inputTable['IMAGE_flux' + apString])
            good_indicies = (seg_sum < 0.5)  & (wht_sum > smallNum) & \
                ((ap_sum > deltaZero) | (ap_sum < -deltaZero))
            print("There are ", sum(good_indicies), " good indicies, out of ", len(seg_sum))
            good_indicies = np.where(good_indicies)
            
        ## This table has all the good apertures
            reducedTable = inputTable[good_indicies]
            
        else:
            reducedTable = inputTable
            print("inputTable",inputTableFile)
        ## Get the x and y coordianates of the apertures
        apX = reducedTable['IMAGE_xcenter'] #/u.pix
        apY = reducedTable['IMAGE_ycenter'] #/u.pix
        
        ## Check if I need to run local depths from the
        ## apertures
        
        if (os.path.isfile(localDepthsFile) == False) or recalculate:
            ## Files doesn't exist or I want to recalculate it  
            
            #############################################################
            ################## RUNNING LOCAL DEPTHS ####################
            
           ## Find the min max so I can create array for local depths
            #maxNumx = max(apX)
            #maxNumy = max(apY)
            #minNumx = min(apX)
            #minNumy = min(apY)

            ## Why?!
            #xmax = max(apX) + min(apX)*3.0
            #ymax = max(apY) + min(apY)*3.0
            print("debugggg", type(apX))
            xmax = max(apX)# - min(apX)#*3.0
            ymax = max(apY)# - min(apY)#*3.0
            
            numX = np.ceil(xmax/step)
            numY = np.ceil(ymax/step)

            x = min(apX) + np.arange(numX)*step ## modifed 21/9/2018 to add min(apX).
            y = min(apY) + np.arange(numY)*step
            print("Step = ", step, " numx, y = ", numX, numY)
            print("Max = ", xmax, ymax, max(apX), min(apX))
            
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

            #            x = x[7000:7010]
            #            y = y[7000:7010]
            
            ## Now run local depths at those points
            depthsLocalFull, maskArray = local_depths(reducedTable, apString, x, y, numApertures, zeropoint = zeropoint, mask = True, sigmaClip = 3.0)#, plot = plotDir + filterName + '_' + str(numApertures))
            
            ## remove points that lie off the image
            #good_ind = depthsLocal > 0.0
            #x = x[good_ind]
            #y = y[good_ind]
            #depthsLocalFull = depthsLocal[good_ind]
            
            ## Now save these results for faster calculating/plotting in future
            ## Create a table
            localTable = Table([x, y, depthsLocalFull, maskArray], names = ['x', 'y', 'depths', 'mask'], dtype = ['f4', 'f4', 'f4', 'f4'])
            localTable.write(localDepthsFile, overwrite=True)
            print("Local depths saved to ", localDepthsFile)
            
        else:
            # simply restore the results
            localTable = Table.read(localDepthsFile)
            
                ## extract the depths
            #if 'depth' in keys:
            #        #print "Yes!"
            #    depthsLocalFull = localTable['depth']
            #else:
            #        #print "Calculating depths here. "
            #    depthsLocalFull = localTable['depths']
                

                
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
                        print(regDir)
                        
                        fff = refimage.find('HSC-R')
                        if fff > -1:
                            circlesFile = regDir + 'HSC_circle_cdfs_R_xy.reg'
                        else:
                            circlesFile = regDir + 'HSC_circle_cdfs_xy.reg'
                            
                        good_indicies = mask_column(x, y, maskreg, tokeep = True, hsc = hsc, xy = True, circlesFile = circlesFile)
                    else:
                        print(ra, dec, maskreg, hsc)
                        good_indicies = mask_column(ra, dec, maskreg, tokeep = True, hsc = hsc)
                    
                
                #print "There are ", x.shape, y.shape, " local depth positions..."

            print('There are {0} good indicies'.format(np.sum(good_indicies)))

            good_indicies = good_indicies & np.logical_not(np.isnan(localTable['depths']))
            
            finalTable = localTable[good_indicies]
            apXregion = finalTable['x'] #/u.pix
            apYregion = finalTable['y'] #/u.pix
            depthsLocal = finalTable['depths'] #depthsLocalFull[good_indicies]

            ii = np.logical_not(np.isnan(depthsLocal))
            
            print('Check for NANS', np.isnan(np.sum(depthsLocal)), np.sum(ii), depthsLocal[ii])
            
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
                
#                plt.scatter(x[good_indicies], y[good_indicies], s = 11, linewidth = 0.1, facecolor = 'none', edgecolor = 'k', alpha = 0.5)
                plt.colorbar(sc)
                plt.title('Local depths for filter {0}\n Aperture diameter is {1:.1f}as'.format(filterName, apDiAS))
                
                # now make a histogram to go underneath, do this differently if in strips!
                ax = plt.subplot(gs[1])
                binwidth = 0.01 #0.01
                #low = np.round(magmin*100)/100.0
                low = magmin
                
                bins = np.arange(low+binwidth/2.0, magmax+binwidth, binwidth)
                print('Bins defined here ', bins)
                
                if strips:
                    ## plot subsets of the results
                    ## and colour based on strip/gap 
                    print('In strips', bins.shape)
                    
                    ## 1) get the strips
                    print("allgood")
                    ss = return_instrips(apXregion, apYregion)
                    print("allgood2")
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
                    print('Definiting msoothed here, bins ', smoothed.shape, bins.shape)
                    
                else:

                    #print('Strips = false?')
                    print(depthsLocal)
                    print(magmin, magmax)
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
            print("The median local depth is then: MAD = {0:.2f}, Gaussian = {1:.2f}".format(medianLocalDepth, global_depth[ri, ai]))
            median_local_depth[ri, ai] = medianLocalDepth

            if region == 'full':
                # Mode depth
                #            mod = (histy == np.max(histy))
                #mod = (smoothed == np.max(smoothed))
                mod = np.argmax(smoothed)
                smallbins = bins[0:-1]
                print(smoothed.shape, smallbins.shape, bins.shape, mod.shape)
                print(bins)
                print('smallbinds ', smallbins)
                print(smallbins[mod])
                
                mode = smallbins[mod]+binwidth/2.0
                print(mode)
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
                        print(peaks)
                        print('Smoothed!', smoothx[peaks])
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


def local_depths(cleanTable, apString, x, y, numApertures, zeropoint = -99.0, verbose=False, mask = False, sigmaClip = 3.0, plot = 'none', regFile = 'none', fitGauss = False):
    ''' Code to find the local depth around a given aperture coordinate, or coodinates.
    Using the closest numApertures apertures.'''

    #plot = 'nirspec_test'
    #verbose = False
    
    #from astropy import units as u
    import numpy as np
    import matplotlib.pyplot as plt
    
    ## Check x and y have the same size
    if len(x) != len(y):
        print("Error, the x and y arrays have a different length.")
        print(len(x), len(y))
        exit()
    
    ## extract the required columns from the input table
    apX = np.array(cleanTable['IMAGE_xcenter']) #/u.pix
    apY = np.array(cleanTable['IMAGE_ycenter']) #/u.pix
    number = np.arange(apX.size)
    
    allApertureResults = np.array(cleanTable['IMAGE_flux'+apString])
    
    ## do two lots of sigma clipping
    
    ## Make a nice output array
    localDepths = np.zeros(x.size)
    maskArray = np.zeros(x.size)
    
    print("There are ", x.size, " positions to calculate depths for.")
    #print cleanTable
    #print allApertureResults
    #print localDepths
    #exit()
    
    if regFile != 'none':
        tf = open(regFile, 'w')
        print('Masking with reg file, ', regFile)

    if plot != 'none':
        import matplotlib.backends.backend_pdf
        plotname = plot + '_check_mad.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(plotname)

    # testing
    #x = x[1000:1100]
    
    # get a minimum separation
    diffx = np.min(np.abs(y - np.roll(y, 1)))
    print('the miminum x separation is ', diffx)
    
    
    
    # loop through the positions
    for xi, xpos in enumerate(x):
        #for yi, ypos in enumerate(y):
        ypos = y[xi]
        
        #if yi > 4000:
        #    print "Initialising loop."
            
        ## calculate the radius        
        deltaX = apX - xpos
        deltaY = apY - ypos
        radius = np.sqrt(deltaX*deltaX + deltaY*deltaY)
        
        ## sort this array, and then the table
        #sortedIndicies = np.argsort(radius)
        ## do this faster!
        idx = np.argpartition(radius, numApertures)
        useIndicies = idx[0:numApertures]
        
        ## Check that the radius here is close!!
        #print "The radius is ", radius[useIndicies]
        #print "the pos is ", xpos, ypos
        apRedX= apX[useIndicies]
        apRedY= apY[useIndicies]
        
        #print numApertures
        if regFile != 'none':
            for i in range(apRedX.size):
                tf.write('circle\t{0}\t{1}\t6\n'.format(apRedX[i], apRedY[i]))
        
        # this return the indicies of the lowerest numApertures apertures
        # not necessarily sorted
        
        #if yi > 4000:
        #    print "Is it the sort?", sortedIndicies.size
  
        ## Take the closest XX
        #useIndicies = sortedIndicies[0:numApertures]
        
        #sortedRadius = radius[useIndicies]
        #sortedTable = cleanTable[useIndicies]
        #print "the size is ", sortedRadius.size
        
        ## extract the data here
        apertureResults = allApertureResults[useIndicies]
        #print "Now working with {:4f} results.".format(apertureResults.size)
        smallRadius = radius[useIndicies]
        sortedRadius = smallRadius[np.argsort(smallRadius)]
        if verbose:
            print("The largest radii", sortedRadius[-4:])
        
        
        ## do a check to see if there are actually apertures there
        ## or if we are at the edge of the image
        
        #############################################
        ## now calculate the depth!
        ## plot the results
        
            
        ## First clip the data to remove outliers
        ## that will skew the SD
        medianFlux = np.median(apertureResults)
        mad = np.median(abs(apertureResults - medianFlux))
        sigma_mad = 1.4826*mad

        
        
        if (plot != 'none'):
            fig = plt.figure()
            bins = np.arange(-3.0*sigma_mad, 3.0*sigma_mad, sigma_mad/5.0)
            #print bins
            #exit()
            n,bin,patches = plt.hist(apertureResults, bins = bins, facecolor = 'green', alpha = 0.75)
            
            # split by region
            north = (apRedX > np.median(apRedX)) & (apRedY > np.median(apRedY))
            #n,bin,patches = plt.hist(apertureResults[north], bins = bins, facecolor = 'yellow', alpha = 0.75)
            south = (apRedX < np.median(apRedX)) & (apRedY < np.median(apRedY))
            #n,bin,patches = plt.hist(apertureResults[south], bins = bins, facecolor = 'blue', alpha = 0.75)
            
            # plot the median etc
            plt.plot([medianFlux, medianFlux], [0, max(n)])
            plt.plot([medianFlux-sigma_mad, medianFlux-sigma_mad], [0, max(n)], color = 'k')
            plt.plot([medianFlux+sigma_mad, medianFlux+sigma_mad], [0, max(n)], color = 'k')
            plt.xlim([medianFlux -5.0*sigma_mad, medianFlux + 5.0*sigma_mad])
            
        
        if verbose:
            print("The mad sigma is ", sigma_mad, " after first run, mag = {0:.2f}".format(-2.5*np.log10(5.0*sigma_mad) + zeropoint))
            
            #if yi > 4000:
            #        print "Searching for apertures in this sigma range. "
            
            #print "the mad = ", sigma_mad

        if sigma_mad > 1E-15:
            
            if sigmaClip > 0.0:
                # First sigma clip
                good_indicies = (apertureResults > medianFlux - sigmaClip*sigma_mad) & \
                    (apertureResults < medianFlux + sigmaClip*sigma_mad)
                
                clippedResults = apertureResults[good_indicies]
                if verbose:
                    print("After the first clipping we have ", np.sum(good_indicies), " good indicies.")
                
        ## Second sigma clip
                medianFlux = np.median(clippedResults)
                mad = np.median(abs(clippedResults - medianFlux))
                sigma_mad = 1.4826*mad
                
                good_indicies_two = (clippedResults > medianFlux - sigmaClip*sigma_mad) & \
                    (clippedResults < medianFlux + sigmaClip*sigma_mad)
                
                if verbose:
                    print("After the second clipping we have ", np.sum(good_indicies_two), " good indicies.")
                    
                finalClippedResults = clippedResults[good_indicies_two]
                if (plot != 'none'):
                    n,bin,patches = plt.hist(finalClippedResults, bins = bins, facecolor = 'red', alpha = 0.75)
                
                    if fitGauss:
                        from scipy.stats import norm
                        import matplotlib.mlab as mlab
                        
                    # fit a gaussian to these points!
                        
                        (mu, sigma) = norm.fit(clippedResults)                
                        yll = mlab.normpdf(bins, mu, sigma)
                        
                        # normalise this
                        plt.plot(bins, np.sum(n)*yll/np.sum(yll), 'k', linewidth=3)
                        print ("The gauss sigma is ", sigma, " after first run, mag = {0:.2f}".format(-2.5*np.log10(5.0*sigma) + zeropoint))

            
            ## Now calculate the final mad sigma value
                medianFlux = np.median(finalClippedResults)
                mad = np.median(abs(finalClippedResults - medianFlux))
                sigma_mad = 1.4826*mad
                if verbose:
                    print("After second clipping... = ", sigma_mad)
                    print("The depth = {0:.2f}".format(return_mag(sigma_mad, zeropoint, sigma = 5.0)))
                    
            
            
        #######################################
        # if not clip then just use this.
        #######################################
        
        if float(zeropoint) > -1.0:
                ## convert to magnitudes
            #print sigma_mad, xpos, ypos
            localDepths[xi] = return_mag(sigma_mad, zeropoint, sigma = 5.0)
                #print "The depth is then: {0:.2f} from {1} apertures.".format(return_mag(sigma_mad, zeropoint, sigma = 5.0), finalClippedResults.size)
            
        else:
                ## keep as fluxes
                #print "The depths are going to be in fluxes, not magnitudes"
            localDepths[xi] = sigma_mad
                #if yi > 4000:


        if np.isnan(localDepths[xi]):
            print('NAN here', sigma_mad)
            exit()

        if sortedRadius[0] < diffx*0.1:
            ## There are no apertures nearby...
            ## Set a masked value (for plotting)
            maskArray[xi] = 1
            #print('Masking this as good ', diffx*10, sortedRadius[0])
            
## set the depth to -99
         #   localDepths[xi] = -np.abs(localDepths[xi])
            #print "Here ", localDepths
            #    print "I have added to the final array."
            #exit()
            
        if xi % 1000 == 0:
            print("At position ", xi)

        if (plot != 'none'):
            pdf.savefig(fig)
            plt.close()
    
    if regFile != 'none':
        tf.close()

    if plot != 'none':
        pdf.close()
        print('Plot at ', plotname)

#    print(localDepths)
    # get rid of nans
    bad_indicies = np.isnan(localDepths)
    if np.any(bad_indicies):
        print('Fixing NANs in the local_depth code')
        localDepths[bad_indicies] = -99.0
        #print localDepths[xi], xpos, ypos, sigma_mad

    if mask: 
        return localDepths, maskArray
    else:
        print("CONSIDER updating code to use masked depths.")
        return localDepths

def return_mag(flux_counts, zeropoint, sigma = 1.0):
    
    import math
    
    if flux_counts < 1E-15:
        return -99.0
    else:
        return -2.5*math.log(sigma*flux_counts, 10.0) + zeropoint
    
def grid_depths(gridTable, x, y, faster = True, verbose = False, nearby = False):
    
   ''' Code to find the closest depth measurement from my previous analysis. Faster than truely local depths '''
   
   import numpy as np
   
   xgrid = gridTable['x']
   ygrid = gridTable['y']
   keys = gridTable.colnames
   #print keys
   #print "Grid params."
   #print len(xgrid), len(ygrid)
   #print xgrid[0], ygrid[0], xgrid[1], ygrid[1]
   
   depthsOverField = gridTable['depths']
   
   ## Make an output array
   depthArray = np.zeros(x.size)
   depthArray[:] = -99.0
   
   if faster:
       if verbose:
           print("Using faster method.")
           print("Input array size is ", x.size)
       deltay = np.min(ygrid)
       deltax = np.min(xgrid)
      
    #print "The delta is ", deltax, deltay
       #print deltax, np.max(xgrid)
       #print deltay, np.max(ygrid)
       #print xgrid[0:10]
       #print ygrid[0:10]
    ## loop through the grid instead of each object
       for xi in range(xgrid.size):
           
           xmin = xgrid[xi] - deltax
           xmax = xgrid[xi] + deltax
           ymin = ygrid[xi] - deltay
           ymax = ygrid[xi] + deltay
           #print xmin, xmax, ymin, ymax
           #exit()
           

           ii = (x > xmin) & (x <= xmax) & (y > ymin) & (y <= ymax)
           
           depthArray[ii] = depthsOverField[xi]
               
           # use this to average over nearby pixels.
           
               
   else:
       
   ## Find the closest point to the objects x and y positions
   ## Loop!
       for xi in range(x.size):
           
       ## make a radius array
           deltax = (xgrid - x[xi])
           deltay = (ygrid - y[xi])
           radius = np.sqrt(deltax*deltax + deltay*deltay)
           mini = np.argmin(radius)
           
       ## try using argpartition
           numpoints = 10
           idx = np.argpartition(radius, numpoints)
           
           if nearby:
               
               mini = idx[0:numpoints]
               print("The nearby depths are = ", depthsOverField[mini])
               print("Before = ", depthsOverField[mini][0])
           
           
       #print "The closest point is ", xgrid[mini], ygrid[mini], " to ", x[xi], y[xi]
           depthArray[xi] = depthsOverField[mini][0]
       #exit()
   
   return depthArray
       
       
   
def grid_psf(gridTable, x, y, faster = True, verbose = False, nearby = False):
    
   ''' Code to find the closest psf measurement from PSFEx map '''
   '''Note from Rohan: copy and pasted above grid_depths, changed a couple things to adapt to PSFEx map.'''
   import numpy as np
   
   xgrid = gridTable['x']
   ygrid = gridTable['y']
   keys = gridTable.colnames
   #print keys
   #print "Grid params."
   #print len(xgrid), len(ygrid)
   #print xgrid[0], ygrid[0], xgrid[1], ygrid[1]
   
   PSFsOverField = gridTable['ef_2.0']
   
   ## Make an output array
   psfArray = np.zeros(x.size)
   psfArray[:] = -99.0
   
   if faster:
       if verbose:
           print("Using faster method.")
           print("Input array size is ", x.size)
       deltay = np.min(ygrid)
       deltax = np.min(xgrid)
      
    #print "The delta is ", deltax, deltay
       #print deltax, np.max(xgrid)
       #print deltay, np.max(ygrid)
       #print xgrid[0:10]
       #print ygrid[0:10]
    ## loop through the grid instead of each object
       for xi in range(xgrid.size):
           
           xmin = xgrid[xi] - deltax
           xmax = xgrid[xi] + deltax
           ymin = ygrid[xi] - deltay
           ymax = ygrid[xi] + deltay
           #print xmin, xmax, ymin, ymax
           #exit()
           

           ii = (x > xmin) & (x <= xmax) & (y > ymin) & (y <= ymax)
           
           psfArray[ii] = PSFsOverField[xi]
               
           # use this to average over nearby pixels.
           
               
   else:
       
   ## Find the closest point to the objects x and y positions
   ## Loop!
       for xi in range(x.size):
           
       ## make a radius array
           deltax = (xgrid - x[xi])
           deltay = (ygrid - y[xi])
           radius = np.sqrt(deltax*deltax + deltay*deltay)
           mini = np.argmin(radius)
           
       ## try using argpartition
           numpoints = 10
           idx = np.argpartition(radius, numpoints)
           
           if nearby:
               
               mini = idx[0:numpoints]
               print("The nearby PSFs are = ", PSFsOverField[mini])
               print("Before = ", PSFsOverField[mini][0])
           
           
       #print "The closest point is ", xgrid[mini], ygrid[mini], " to ", x[xi], y[xi]
           depthArray[xi] = PSFsOverField[mini][0]
       #exit()
   
   return psfArray
       
       
   
