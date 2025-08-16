import os
import re
import warnings
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from datetime import datetime
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from astropy.visualization import ZScaleInterval

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
        Also creates the subsets and returns them for future use.

    INPUT(s)
        cats[list(str)]     list of paths to the cats you want to open. 
                            Typically, the ones you just made by looping through image_depths.
    OUTPUT(s)
        subsets(list)       List of subset paths  """

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
            catDir = os.path.dirname(cat)
            catName = os.path.basename(cat)[:-len('.fits')]
            subsetPath = os.path.join(catDir, catName + '_subset.fits')

            # skip if subset already exists
            if os.path.isfile(subsetPath):
                subsets.append(subsetPath)
                continue

            # create subset
            table = Table.read(cat)
            subset = table[table["MAG_APER"][:, 1] > 50]  # criterion
            subset.write(subsetPath, overwrite=overwrite)
            subsets.append(subsetPath)

        files_str = " ".join(subsets)
        command += f"{files_str} &"
        print(command)
        #os.system(command) #TODO: turn me back on when needed! 

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

######################## Plotting and Figures ##################################################

def detections_fig_3by2(cutoutPaths, catPaths, subsetPaths, verbose=True):
    """ 
    Makes a figure showing the science image, detections,
    and subset detections for two filters in a 3x2 formation. 
    Ready-made subsets should be passed to this function. 

    cutoutPaths(list[str]):     A list of paths to the science cutouts you want to show e.g. 
                                raid/scratch/hullyott/cataloguing/DepthsTestDir/data/COSMOS/cutouts/
                                UVISTA_J_DR6_1495_22_size500.fits

    catPaths(list[str]):        A list of paths to the detections catalogues made in the cutout regions
                                e.g. /raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/catalogues/
                                d500JH_cutout.fits

    subsetPaths(list[str]):     A list of paths to the subsets of detections made in the cutout regions
                                e.g. /raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/catalogues/
                                d500JH_cutout_subset.fits
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    from astropy.visualization import ZScaleInterval

    ### initialise fig
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval() # z-scale all images
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    for i, cutoutPath in enumerate(cutoutPaths):
        cutPlotTitle = os.path.basename(cutoutPath)
        cutoutData, _, _ = get_fits_data(cutoutPath, verbose=verbose)
        vmin, vmax = zscale.get_limits(cutoutData)
        ax = axes[i, 0]
        ax.imshow(cutoutData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"{cutPlotTitle}", fontsize=10)
        ax.xaxis.set_major_formatter(formatter)
        ax.axis('off')

    for i, catPath in enumerate(catPaths):
        catPlotTitle = os.path.basename(catPath)
        catTable = Table.read(catPath)
        ra = catTable['ALPHA_J2000']
        dec = catTable['DELTA_J2000']
        ax_cat = axes[i, 1]
        ax_cat.scatter(ra, dec, color='blue', alpha=0.5, s=1)
        ax_cat.set_title(f"{catPlotTitle}", fontsize=10)
        ax_cat.xaxis.set_major_formatter(formatter)
        ax_cat.invert_xaxis()  # invert RA axis to match sky convention
        ax_cat.grid(True)

    for i, subsetPath in enumerate(subsetPaths):
        subsetPlotTitle = os.path.basename(subsetPath)
        subsetTable = Table.read(subsetPath)
        ra = subsetTable['ALPHA_J2000']
        dec = subsetTable['DELTA_J2000']
        ax_sub = axes[i, 2]
        ax_sub.scatter(ra, dec, color='green', alpha=0.5, s=1)
        ax_sub.set_title(f"{subsetPlotTitle}", fontsize=10)
        ax_sub.xaxis.set_major_formatter(formatter)
        ax_sub.invert_xaxis()  # invert RA axis to match sky convention
        ax_sub.grid(True)

    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    #plt.show()

def detections_fig(cutoutPaths, catPaths, subsetPaths, badMags,saveFig=True, overwrite=True, verbose=True):
    """ 
    Makes a figure showing the science imagea and measurements,
    overlayed by subset detections for two filters in a 2x2 formation. 
    Ready-made subsets should be passed to this function. 

    cutoutPaths(list[str]):     A list of paths to the science cutouts you want to show e.g. 
                                raid/scratch/hullyott/cataloguing/DepthsTestDir/data/COSMOS/cutouts/
                                UVISTA_J_DR6_1495_22_size500.fits

    catPaths(list[str]):        A list of paths to the detections catalogues made in the cutout regions
                                e.g. /raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/catalogues/
                                d500JH_cutout.fits

    subsetPaths(list[str]):     A list of paths to the subsets of detections made in the cutout regions
                                e.g. /raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/catalogues/
                                d500JH_cutout_subset.fits

    badMags(float):             Percentage of bag-mag detections within region. bad-mag == MAG_APER[1]>50
    """

    ### initialise fig
    fig, axes = plt.subplots(nrows=len(cutoutPaths), ncols=2, figsize=(12, 8))
    axes = np.atleast_2d(axes)  # forces 2D indexing
    zscale = ZScaleInterval() # z-scale all images
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    for i, (cutoutPath, subsetPath, badMag) in enumerate(zip(cutoutPaths, subsetPaths, badMags)):
        cutPlotTitle = os.path.basename(cutoutPath).replace('.fits', '')
        subsetPlotTitle = os.path.basename(subsetPath).replace('_cutout_subset.fits', 'Mag>50')
        cutoutData, _, wcs = get_fits_data(cutoutPath, verbose=verbose)
        subsetTable = Table.read(subsetPath)
        subRA = subsetTable['ALPHA_J2000']
        subDec = subsetTable['DELTA_J2000']
        # Convert RA/Dec -> pixel coordinates
        x_pix, y_pix = wcs.world_to_pixel_values(subRA, subDec)

        vmin, vmax = zscale.get_limits(cutoutData)
        ax = axes[i, 0]
        ax.imshow(cutoutData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"{cutPlotTitle}", fontsize=10)
        ax.scatter(x_pix, y_pix, color='red', alpha=0.5, s=10, label=subsetPlotTitle)
        ax.legend()
        ax.grid(True)
        ax.text(0.0,-50, f"MAG_APER[1]>50: {badMag}%", color='red')
        ax.xaxis.set_major_formatter(formatter)
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("Dec (deg)")
        ax.axis('off')

    for i, (catPath, subsetPath) in enumerate(zip(catPaths, subsetPaths)):
        catPlotTitle = os.path.basename(catPath)
        size = re.findall(r'\d+', catPlotTitle)[0]
        subsetPlotTitle = os.path.basename(subsetPath)
        catTable = Table.read(catPath)
        subsetTable = Table.read(subsetPath)
        ra = catTable['ALPHA_J2000']
        dec = catTable['DELTA_J2000']
        subRA = subsetTable['ALPHA_J2000']
        subDec = subsetTable['DELTA_J2000']

        ax_cat = axes[i, 1]
        ax_cat.scatter(ra, dec, color='blue', alpha=0.5, s=1, label='all measurements')
        ax_cat.scatter(subRA, subDec, color='red', alpha=0.5, s=1, label=subsetPlotTitle)
        ax_cat.legend()
        ax_cat.set_title(f"{catPlotTitle}", fontsize=10)
        ax_cat.xaxis.set_major_formatter(formatter)
        ax_cat.invert_xaxis()  # invert RA axis to match sky convention
        ax_cat.grid(True)

    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    if saveFig and overwrite==True:
        outDir = '/raid/scratch/hullyott/cataloguing/plots/'
        outputPath = outDir + size + "arcsec" + "detections_badMags.pdf"
        try:
            plt.savefig(outputPath, format='pdf')
            print(f"Figure saved as {outputPath}")
        except PermissionError:
            print(f"Error: Permission denied. Could not save the file {outputPath}.")
        except Exception as e:
            print(f"Error: Could not save the file {outputPath}. Reason: {e}")

    #plt.show()

def bg_plotter(cutoutPaths='none', bgMapPaths='none', listBgSubDicts='none', badMags='none', segPaths='none', whtPaths='none', show_bad_mags=False, saveFig=True, overwrite=True, verbose=True):
    """
        
    """

    ### initialise fig
    ncols = 0
    figTitleSubstr = []
    for arg in [cutoutPaths, bgMapPaths, listBgSubDicts, segPaths, whtPaths]:
        if arg != 'none':
            ncols += 1

    for substrlist in [cutoutPaths, bgMapPaths, listBgSubDicts, badMags, segPaths, whtPaths]:
        for substr in substrlist:
            if type(substr) == dict or type(substr) == float or substrlist=='none': # if the list is listBgSubDicts or bad mags or not set
                continue
            substr = os.path.basename(substr)
            substr = substr.replace('.fits', '')
            substr = substr.replace('UVISTA_', '')
            figTitleSubstr.append(substr)
    figTitleSubstr = "_".join(figTitleSubstr).replace(' ', '_')
    figTitleSubstr = list(dict.fromkeys(figTitleSubstr))
    #print(figTitleSubstr)
    breakpoint()
    fig, axes = plt.subplots(nrows=len(listBgSubDicts), ncols=ncols, figsize=(12, 8))
    axes = np.atleast_2d(axes)  # forces 2D indexing

    zscale = ZScaleInterval() # z-scale all images
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)


    for i in range(len(listBgSubDicts)):  # one row per BgSubDict
        col_idx = 0

        # plot cutout region
        if cutoutPaths != 'none':
            cutoutPath = cutoutPaths[i]
            cutPlotTitle = os.path.basename(cutoutPath).replace('.fits', '')
            cutoutData, _, _ = get_fits_data(cutoutPath, verbose=verbose)
            vmin, vmax = zscale.get_limits(cutoutData)
            ax = axes[i, col_idx]
            ax.imshow(cutoutData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
            ax.set_title(f"{cutPlotTitle}", fontsize=10)
            ax.xaxis.set_major_formatter(formatter)
            ax.axis('off')
            col_idx += 1

        # plot bgmap
        if bgMapPaths != 'none':
            bgMapPath = bgMapPaths[i]
            bgMapPlotTitle = os.path.basename(bgMapPath).replace('.fits', '').replace('_cutout_', '')
            bgMapData, header, _ = get_fits_data(bgMapPath, verbose=verbose)
            vmin, vmax = zscale.get_limits(bgMapData)
            ax_bgmap = axes[i, col_idx]
            ax_bgmap.imshow(bgMapData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
            ax_bgmap.set_title(f"{bgMapPlotTitle}", fontsize=10)
            ax_bgmap.xaxis.set_major_formatter(formatter)
            ax_bgmap.axis('off')
            col_idx += 1

        # plot bgSub
        if listBgSubDicts != 'none':
            BgSubDict = listBgSubDicts[i]
            bgSubPath = BgSubDict['Path']
            bgSubName = os.path.basename(bgSubPath).replace('.fits', '').replace('_cutout_', ' ')
            bgSubParams = f"BKSZ: {BgSubDict['back_size']} BKFILT: {BgSubDict['back_filtersize']}"
            bgSubPlotTitle = bgSubName + bgSubParams
            bgSubData, header, _ = get_fits_data(bgSubPath, verbose=verbose)
            vmin, vmax = zscale.get_limits(bgSubData)
            ax_bgsub = axes[i, col_idx]
            ax_bgsub.imshow(bgSubData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
            ax_bgsub.set_title(f"{bgSubPlotTitle}", fontsize=10)
            ax_bgsub.xaxis.set_major_formatter(formatter)
            ax_bgsub.axis('off')
            col_idx += 1

        # plot segmentation maps
        if segPaths != 'none':
            segPath = segPaths[i]
            segPlotTitle = os.path.basename(segPath).replace('.fits', '').replace('_cutout_', '')
            segData, header, _ = get_fits_data(segPath, verbose=verbose)
            vmin, vmax = zscale.get_limits(segData)
            ax_seg = axes[i, col_idx]
            ax_seg.imshow(segData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
            ax_seg.set_title(f"{segPlotTitle}", fontsize=10)
            ax_seg.xaxis.set_major_formatter(formatter)
            ax_seg.axis('off')
            col_idx += 1

        # plot weight maps
        if whtPaths != 'none':
            whtPath = whtPaths[i]
            whtPlotTitle = os.path.basename(whtPath).replace('.fits', '').replace('_cutout_', '')
            whtData, header, _ = get_fits_data(whtPath, verbose=verbose)
            vmin, vmax = zscale.get_limits(whtData)
            ax_wht = axes[i, col_idx]
            ax_wht.imshow(whtData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
            ax_wht.set_title(f"{whtPlotTitle}", fontsize=10)
            ax_wht.xaxis.set_major_formatter(formatter)
            ax_wht.axis('off')
            col_idx += 1

    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    if saveFig and overwrite==True:
        outDir = '/raid/scratch/hullyott/cataloguing/plots/'
        #outputPath = outDir+"background_comparison_"+bgSubParams.replace(":", '').replace(' ', '_')+ ".pdf"
        outputPath = outDir +  figTitleSubstr + bgSubParams.replace(":", '').replace(' ', '_') + ".pdf"
        try:
            plt.savefig(outputPath, format='pdf')
            print(f"Figure saved as {outputPath}")
        except PermissionError:
            print(f"Error: Permission denied. Could not save the file {outputPath}.")
        except Exception as e:
            print(f"Error: Could not save the file {outputPath}. Reason: {e}")

    plt.show()

       
######################## do depths ##################################################
def aperture_photometry_blank(imageName, segMap, whtMap, apSize, gridSeparation=100, pixScale=-99.0, next=0, clean=False, outputFitsName='none', imageDir='', verbose=False, field='NORMAL', overwrite=False):
    """
    Performs aperture photometry across the image at fixed intervals (defined by gridSeperation), 
    without regard for the positions of real astronomical sources.
    This allows estimation of the background noise.

    INPUT(s)
        imageName(str)       FITS image to measure.
	    segMap(str)          Segmentation map (mask for objects).
	    whtMap(str)          Weight map (for noise or exposure info).
	    apSize(arr)          Array of aperture diameters (arcseconds).
	    gridSeparation(int)  Pixel separation between apertures.
	    pixScale(float)      Pixel scale in arcseconds/pixel (reads from header if negative).
	?	next(int)            HDU index in FITS file.
	?	clean(Boolean)       Whether to filter out "bad" apertures.
	    outputFitsName(str)  Path for output table (creates if doesn?t exist).
	    imageDir(str)        Default directory for outputs.
	?	field(str)           Special cleaning rules for certain fields (e.g., NIRSPEC).



     """
    
    print("aperphotblank ", segMap, whtMap, apSize)

    from astropy.io import fits
    from astropy.table import Table, hstack, join
    import time
    #import sep # aperture photometry from SExtractor!
    
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


def aperture_phot_fast(imageData, xArray, yArray, radii, subpix=5):
    
    # Subpix = 5 is what SEXtractor uses.
    import sep
    from astropy.table import Table, Column
    
    data = imageData.byteswap().newbyteorder()    
    for ri, r in enumerate(radii):
        flux, fluxerr, flag = sep.sum_circle(data, xArray, yArray, r, subpix = subpix)

        if ri < 1:
            # create a table of the results
            phot_apertures = Table([xArray, yArray, flux], names = ['xcenter', 'ycenter', 'flux_0'], dtype = ['f4', 'f4', 'f4'])
            
        else:
            newcolumn = Column(name = 'flux_' + str(ri), data = flux, dtype = 'f4')
            phot_apertures.add_column(newcolumn)
    
    return phot_apertures



def image_depth(imagePath, zeropoint, cutouts=[], size='none', back_size=32, back_filtersize=9, apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtPath='NONE', whtType='NONE', IRACapDiametersAS=np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segPath='NONE', outputDir='none', filterName='NONE', numApertures=300, step=200, overwrite=False, inputSex=baseDir+'data/bertin_config/video_mine.sex', strips=False, bgSub=True, mask='none', gridSepAS=3.0):

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
    bgSubDict = {} # name, back_size, back_sizefilter

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
            # dict used in bg_plotter
            bgSubDict['Path'] = bgSubPath
            bgSubDict['back_size'] = back_size
            bgSubDict['back_filtersize'] = back_filtersize
        else:
            outputCatalogue = catDir + 'd' + filterName + '.fits'
            segPath = imageDir + filterName + '_seg.fits'  
            bgMapPath = imageDir  + filterName + '_bgmap.fits'
            bgSubPath = imageDir  + filterName + '_bgsub.fits'
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
        # if aperphotfile doesn't exist, call ap phot blank - outputFitsName is aperPhotFile
        aperture_photometry_blank(bgSubPath, segPath, whtPath, apDiametersAS, gridSeparation = gridSepPixels, clean = True, outputFitsName = aperPhotFile, imageDir = imageDir, field = field, overwrite = overwrite)

    return bgMapPath, bgSubDict, segPath

def get_depths(fieldName, cutouts, size='none', back_size=32, back_filtersize=9, queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir=baseDir+'data/', outputDir='none', overwrite=False, ra_str='none', dec_str='none', verbose=True, saveFig=True):

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
    
    bgMapPaths = [] # used in bg_plotter
    listBgSubDicts = [] # used in bg_plotter
    segPaths = [] # used in bg_plotter

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
            bgMapPath, bgSubDict, segPath = image_depth(imageDir+imageName, zeropoint, cutouts=cutouts, size=size, back_size=back_size, back_filtersize=back_filtersize, whtPath=imageDir+whtName, whtType=whtType, outputDir=outputDir, strips=strips, filterName=filterName, overwrite=overwrite, mask=maskName, gridSepAS=gridSepAS, apDiametersAS=apDiametersAS)

            # used in bg_plotter
            bgMapPaths.append(bgMapPath) 
            listBgSubDicts.append(bgSubDict)
            segPaths.append(segPath)

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
    catPaths = []
    if outputDir == 'none':
        outputDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/catalogues/'
        for catName in os.listdir(outputDir):
            for filt in reqFilters:
                pattern = (rf"d{size}{filt}_cutout\.fits")
                if re.fullmatch(pattern, catName):
                    catPath = os.path.join(outputDir, catName)
                    catPaths.append(catPath)

    subsetPaths = open_cats(catPaths, open_subset=True, overwrite=True)
    cutoutPaths = [] # selection of files that will be passed to detections_fig() (i.e. excl wht files)
    whtPaths = []    
    badMags = [] # all detections/bad_mag measurements

    for subsetPath in subsetPaths:
        for filt in reqFilters:
            subpattern = (rf"{baseDir}depths/catalogues/"
                        rf"d{size}{filt}_cutout_subset\.fits")
            if re.fullmatch(subpattern, subsetPath):
                catPath = f"{baseDir}depths/catalogues/d{size}{filt}_cutout.fits"
                cutoutPath = f"{baseDir}data/COSMOS/cutouts/UVISTA_{filt}_DR6_{ra_str}_{dec_str}_size{size}.fits" # corresponds to the catalogue and subset catalogue
                cutoutPaths.append(cutoutPath)
                catTable = Table.read(catPath)
                subsetTable = Table.read(subsetPath)
                badMag = round(100*(len(subsetTable)/len(catTable)),5)
                badMags.append(badMag)
                print(f"{filt}-bad-proportion: ", badMag, "%")

    for subsetPath in subsetPaths:
        for filt in reqFilters:
            subpattern = (rf"{baseDir}depths/catalogues/d{size}{filt}_cutout_subset\.fits")
            if re.fullmatch(subpattern, subsetPath):
                catPath = f"{baseDir}depths/catalogues/d{size}{filt}_cutout.fits"
                whtPath = f"{baseDir}data/COSMOS/cutouts/UVISTA_{filt}_DR6_wht_{ra_str}_{dec_str}_size{size}.fits" # corresponds to the catalogue and subset catalogue
                whtPaths.append(whtPath)

    detections_fig(cutoutPaths=cutoutPaths, catPaths=catPaths, subsetPaths=subsetPaths, badMags=badMags, saveFig=saveFig, overwrite=overwrite, verbose=verbose)

    bg_plotter(cutoutPaths=cutoutPaths, bgMapPaths=bgMapPaths, listBgSubDicts=listBgSubDicts, badMags=badMags, segPaths=segPaths, whtPaths=whtPaths, saveFig=True, overwrite=True, verbose=verbose)


    return
    

