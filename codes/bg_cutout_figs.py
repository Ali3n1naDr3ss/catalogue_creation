
# Global stuff
baseDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/'
import os
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from matplotlib.ticker import ScalarFormatter
from astropy.visualization import ZScaleInterval

def get_fits_data(fitsImg, verbose=True):
    """ Returns data, header, and wcs of a fits image. """
    if verbose:
        print(f"Getting data for {fitsImg} \n")
    with fits.open(fitsImg) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

    return data, header, wcs

def make_cutout(filters, size_arcsec, centre_ra, centre_dec, dataDir, verbose=False, overwrite=True):
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
        if 'HSC_' in filt:
            imgPath = dataDir + filt + '_DR3.fits'
            origFilePaths.append(imgPath)
            whtPath = dataDir +  filt + '_DR3_wht.fits'
            origFilePaths.append(whtPath)
        else:
            imgPath = dataDir + 'UVISTA_' + filt + '_DR6.fits'
            origFilePaths.append(imgPath)
            whtPath = dataDir + 'UVISTA_' + filt + '_DR6_wht.fits'
            origFilePaths.append(whtPath)
    
    for f, file_ in enumerate(origFilePaths):
        #print(file_)
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

        vistafilt = os.path.basename(file_).split('_')[1] # get filt from filename so it matches to file_ in loop
        hscfilt = os.path.basename(file_).split('_DR')[0]

        if 'HSC_' in file_:
            filt = hscfilt
            outPrefix = dataDir + 'cutouts/'
            outSuffix = '_DR3_' + ra_filename +'_'+ dec_filename +'_'+ 'size' + str(size_arcsec) + '.fits'
            cutoutPath = outPrefix + filt + outSuffix
        elif 'HSC_' not in file_:
            filt = vistafilt
            outPrefix = dataDir + 'cutouts/' + 'UVISTA_' 
            outSuffix = '_DR6_' + ra_filename +'_'+ dec_filename +'_'+ 'size' + str(size_arcsec) + '.fits'
            cutoutPath = outPrefix + filt + outSuffix

        if  '_wht' in file_:
            if 'HSC_' in file_:
                filt = hscfilt
                outWhtPrefix = dataDir + 'cutouts/'
                outWhtSuffix =  '_DR3_wht_' + ra_filename +'_'+ dec_filename +'_'+ 'size' + str(size_arcsec) + '.fits'
                cutoutWhtPath = outWhtPrefix + filt + outWhtSuffix
            elif 'HSC_' not in file_:
                filt = vistafilt
                outWhtPrefix = dataDir + 'cutouts/' + 'UVISTA_' 
                outWhtSuffix = '_DR6_wht_' + ra_filename +'_'+ dec_filename +'_'+ 'size' + str(size_arcsec) + '.fits'
                cutoutWhtPath = outWhtPrefix + filt + outWhtSuffix

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

    if len(outFiles) != 0:
        print("\nThe following cutout image and weight files have been made: \n", outFiles)
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

def prep(size_arcsec=100, centre_ra=149.5, centre_dec=2.2, filters=['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y'], dataDir=baseDir+'data/COSMOS/'):
    import re

    cutoutDir = os.path.join(dataDir, f'cutouts/')
    if os.path.isdir(cutoutDir) == False:
        os.system('mkdir ' + cutoutDir)
    cutoutFile = os.path.join(cutoutDir, 'cutoutNames.txt')
    if not os.path.isfile(cutoutFile):
        with open(cutoutFile, "w") as f:
            pass  # just create and close, leaving it empty

    imagePaths = []
    existingCutouts = []
    necessaryCuts = []

    centre_ra_str = str(centre_ra).replace('.', '')
    centre_dec_str = str(centre_dec).replace('.', '')

    for filt in filters:
        for file_ in os.listdir(dataDir):
            pattern = f"UVISTA_{filt}_DR6.fits"
            HSCpattern = f"{filt}_DR3.fits"
            if file_ == pattern:
                imagePath = os.path.join(dataDir, file_)
                imagePaths.append(imagePath)
            elif file_ == HSCpattern:
                imagePath = os.path.join(dataDir, file_)
                imagePaths.append(imagePath)

    # check if the files have been made already
    # determine which files are necessary for desired run
    for filt in filters:
        if 'HSC_' in filt:
            necessaryCut = f"{cutoutDir}{filt}_DR3_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
            necessaryWht = f"{cutoutDir}{filt}_DR3_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
            if necessaryCut not in necessaryCuts:
                necessaryCuts.append(necessaryCut)
            if necessaryWht not in necessaryCuts:
                necessaryCuts.append(necessaryWht)
        else:
            necessaryCut = f"{cutoutDir}UVISTA_{filt}_DR6_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
            necessaryWht = f"{cutoutDir}UVISTA_{filt}_DR6_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}.fits"
            if necessaryCut not in necessaryCuts:
                necessaryCuts.append(necessaryCut)
            if necessaryWht not in necessaryCuts:
                necessaryCuts.append(necessaryWht)

    # check if the necessary files have already been made
    if os.path.isfile(cutoutFile):
        with open(cutoutFile, 'r') as c:
            print("cutoutFile", cutoutFile)
            if os.path.getsize(cutoutFile) == 0:
                print("Cutout file is empty")
            for line in c:
                line = line.strip()
                for filt in filters:
                    if 'HSC_' in filt:
                        pattern = (
                            rf"{cutoutDir}"
                            rf"{filt}_DR3_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                        if re.fullmatch(pattern, line):
                            existingCutouts.append(line)
                        whtpattern = (
                            rf"{cutoutDir}"
                            rf"{filt}_DR3_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                        if re.fullmatch(whtpattern, line):
                            existingCutouts.append(line)

                    elif 'HSC_' not in filt:
                        pattern = (
                            rf"{cutoutDir}"
                            rf"UVISTA_{filt}_DR6_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                        if re.fullmatch(pattern, line):
                            existingCutouts.append(line)
                        whtpattern = (
                            rf"{cutoutDir}"
                            rf"UVISTA_{filt}_DR6_wht_{centre_ra_str}_{centre_dec_str}_size{size_arcsec}\.fits")
                        if re.fullmatch(whtpattern, line):
                            existingCutouts.append(line)

    elif os.path.isfile(cutoutFile) != True:
        print(f">>>>>>>>>>>>>>>>>> WARNING: Could not find {cutoutFile}...")

    if existingCutouts != necessaryCuts: # if the correct cutouts do not exist, make them
        unmatchedCutouts = make_cutout(filters, size_arcsec, centre_ra, centre_dec, dataDir) # cutout Paths
    else:
        return necessaryCuts

def bg_subtraction(zeropoint, cutouts=[], size='none', back_size=32, back_filtersize=9, apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtPath='NONE', whtType='NONE', IRACapDiametersAS=np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segPath='NONE', outputDir='none', numApertures=300, step=200, overwrite=False, inputSex=baseDir+'data/bertin_config/video_mine.sex', strips=False, bgSub=True, mask='none', gridSepAS=3.0):

    """
    Modified from hu_depth_codes which uses Bowler's image_depths()

    Perform background modelling and subtraction using SourceExtractor. 
    If blank-sky aperture photometry results file doesn't exist, call aperture_photometry_blank() to
    estimate the blank-sky noise and write the file. Call extract_local_ddepths() to give local and 
    global depths. Returns the paths to the bg-map and bg-subbed images and seg map.
    """

    warnings_triggered = 0
    cd1Test = False

    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    
    # Tinsley SE env
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/share/sextractor'

    # define the output files
    if outputDir == 'none':
        outputDir = os.path.join(baseDir, 'depths/COSMOS_bgFigs/')
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

    if cutouts == []:
        print(">>>>>>>>>> WARNING: No cutout file names inc cutout list...")

    ###### seg map, SE, and Bkg Subtraction ###########################################
    bgSubDict = {} # name, back_size, back_sizefilter
    
    for cutout in cutouts:
        if 'HSC_' in cutout:
            filterName = os.path.basename(cutout).split('_DR')[0]
        elif 'HSC_' not in cutout:
            filterName = os.path.basename(cutout).split('_')[1]
        # get the pixel scale
        hdulist = fits.open(cutout)
        imageHeader = hdulist[0].header
        imageHeaderKeys = list(imageHeader.keys())

        # check CD1_1 is in Header - if using a new cutout, this has already been checked but good to check again
        cd1Test = ['CD1_1' in imageHeaderKeys]
        print(f">>>> WARNING: Passing file to SE. Header includes CD1_1?: {cd1Test} Is this the right file?", '\n', cutout, '\n') 
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

            size = str(size)
            back_size = str(back_size)
            back_filtersize = str(back_filtersize)
            segPath = imageDir + filterName + size + '_'+ back_size + '_' + back_filtersize +'_cutout_seg.fits'  
            bgSubPath = imageDir  + filterName + size + '_'+ back_size + '_' + back_filtersize + '_cutout_bgsub.fits'
            bgMapPath = imageDir  + filterName + size + '_'+ back_size + '_' + back_filtersize + '_cutout_bgmap.fits'
            outputCatalogue = catDir + 'd' + filterName + '_'+ back_size + '_' + back_filtersize + size + '_cutout.fits'
            # dict used in bg_plotter
            bgSubDict['Path'] = bgSubPath
            bgSubDict['back_size'] = back_size
            bgSubDict['back_filtersize'] = back_filtersize

        keywordsbase = ' -CATALOG_TYPE FITS_1.0 -CATALOG_NAME ' + \
                       outputCatalogue + \
                       ' -MAG_ZEROPOINT '+ str(zeropoint) + \
                       ' -WEIGHT_TYPE ' + whtType + \
                       ' -WEIGHT_IMAGE ' + whtPath
        # -CHECKIMAGE_TYPE "BACKGROUND" for background map -BACKGROUND for bkg subtracted
        keywords = keywordsbase + \
                   ' -CHECKIMAGE_TYPE "BACKGROUND,-BACKGROUND,SEGMENTATION" '\
                   +'-CHECKIMAGE_NAME "' + \
                   bgMapPath + ',' + bgSubPath + ',' + segPath + '" ' + \
                   '-PHOT_APERTURES ' + apStringPix + \
                   ' -BACK_SIZE ' + str(back_size) + \
                   ' -BACK_FILTERSIZE ' +str(back_filtersize)
                   # "', '" needed to interpret bgMApPAth etc as str

# with -CATALOG_TYPE FITS_1.0, the input parameter settings are saved in the header.
# NOTE: the fits option can't handle array outputinformation such as MAG APER, FLUX RADIUS if more than one value! - SE for Dummies 
# PHOT_APERTURES specifies the number and size of apertures (in pix)

        command = '/usr/local/sextractor/bin/sex '+ cutout + ' -c ' + inputSex + ' -PARAMETERS_NAME ' + param_file + ' ' + keywords
        
        if os.path.isfile(bgSubPath) == False or os.path.isfile(segPath) == False or overwrite:
            
            print("The SEG and BG subtracted map do not exist, or overwrite:True.  Running SE like so: \n")
            print("Running image_depths() SE at: ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            print(command)
            os.system(command)

        else:
            print(f"The SEG and/or BG subtracted map exist at: {segPath, bgSubPath} \n")

def plot_bgmap(bgMapPath, savePath=None, cmap='binary', vmin=None, vmax=None):
    """
    Make a figure of the background map produced by Source Extractor.
    
    Parameters
    ----------
    bgMapPath : str
        Path to the bgmap FITS file.
    savePath : str, optional
        Path to save the figure (if None, just shows it).
    cmap : str
        Colormap for imshow (default 'viridis').
    vmin, vmax : float, optional
        Limits for color scale. If None, auto-scaled.
    """
    # load the fits bgmap
    with fits.open(bgMapPath) as hdul:
        data = hdul[0].data

    # Apply z-scaling
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(data)

    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    plt.figure(figsize=(8, 6))

    plt.imshow(data, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar(label="Background level")
    plt.title(f"Background map: {os.path.basename(bgMapPath)}")
    plt.xlabel("X (pixels)")
    plt.ylabel("Y (pixels)")
    
    if savePath:
        plt.savefig(savePath, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved background map figure to {savePath}")
    else:
        plt.show()

filt = 'Y'
size_arcsec = 500 # TODO: testing on different sizes
back_size = 4
back_filtersize = 1
cutouts = prep(size_arcsec=size_arcsec, centre_ra=149.5, centre_dec=2.2, filters=[f'{filt}']) # TODO: doesn't work for len(filter lists)>1
bg_subtraction(zeropoint=30, cutouts=cutouts, size=size_arcsec, back_size=back_size, back_filtersize=back_filtersize)

bgMapDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/COSMOS_bgFigs/images/'
size_arcsec = str(size_arcsec)
back_size = str(back_size)
back_filtersize = str(back_filtersize)

bgMapName = f'{filt}{size_arcsec}_cutout_bgmap.fits'
bgMapPath = bgMapDir + bgMapName
savePath = savePath=f"/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/COSMOS_bgFigs/bg_maps/{filt}_{back_size}_{back_filtersize}_bgmap.png"

plot_bgmap(bgMapPath, savePath=savePath)



