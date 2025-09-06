import os
import re
import warnings
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from datetime import datetime
from astropy.table import Table
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
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
        #os.system(command) #TODO: turn on

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

######################## Plotting and Figures ##################################################

def detections_fig(imagePaths, catPaths, subsetPaths, badMags, saveFig=True, overwrite=True, show_fig=False, verbose=True):
    """ 
    Makes a figure showing the science imagea and measurements,
    overlayed by subset detections for two filters in a 2x2 formation. 
    Ready-made subsets should be passed to this function. 

    imagePaths(list[str]):     A list of paths to the science full-size/cutouts you want to show e.g. 
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
    import math
    import matplotlib
    matplotlib.use("TkAgg")  # or "QtAgg"
    import matplotlib.pyplot as plt

    if any('cutout' in path for path in imagePaths):
        full_size = False
    else:
        full_size = True
        step = 1

    rows_per_fig = 3  

    # split imagePaths etc. into chunks of 3
    n_batches = math.ceil(len(imagePaths) / rows_per_fig)

    for batch_idx in range(n_batches):
        start = batch_idx * rows_per_fig
        end = start + rows_per_fig
        
        # slice the relevant chunk of data
        image_chunk   = imagePaths[start:end]
        subset_chunk  = subsetPaths[start:end]
        badMag_chunk  = badMags[start:end]
        cat_chunk     = catPaths[start:end]
        subsetC_chunk = subsetPaths[start:end]

        # initialise fig for this batch
        fig, axes = plt.subplots(
            nrows=len(image_chunk), ncols=2, figsize=(12, 4*len(image_chunk))
        )
        axes = np.atleast_2d(axes)   # always 2D for indexing
        axes = axes.reshape(-1, 2)

        #axes = np.atleast_2d(axes)  # forces 2D indexing
        axes = axes.reshape(-1, 2)  # force (nrows, 2) shape
        zscale = ZScaleInterval() # z-scale all images
        formatter = ScalarFormatter(useMathText=False)
        formatter.set_scientific(False)
        formatter.set_useOffset(False)

        print("\n>>>>>>> Getting data for detections figures. May take a few mins. Started at: ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        for i, (imagePath, subsetPath, badMag) in enumerate(zip(image_chunk, subset_chunk, badMag_chunk)):
            imagePlotTitle = os.path.basename(imagePath).replace('.fits', '')
            if 'cutout' in os.path.basename(subsetPath):
                subsetPlotTitle = os.path.basename(subsetPath).replace('_cutout_subset.fits', 'Mag>50')
            else:
                subsetPlotTitle = os.path.basename(subsetPath).replace('_subset.fits', 'Mag>50')

            imageData, _, wcs = get_fits_data(imagePath, verbose=verbose)
            if full_size:
                imageData = imageData[::step, ::step] # make number of points manageable!
            subsetTable = Table.read(subsetPath)
            subRA = subsetTable['ALPHA_J2000']
            subDec = subsetTable['DELTA_J2000']
            # Convert RA/Dec -> pixel coordinates
            x_pix, y_pix = wcs.world_to_pixel_values(subRA, subDec)

            vmin, vmax = zscale.get_limits(imageData)
            ax = axes[i, 0]
            ax.imshow(imageData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
            ax.set_title(f"{imagePlotTitle}", fontsize=10)
            ax.scatter(x_pix, y_pix, color='red', alpha=0.5, s=1, label=subsetPlotTitle)
            ax.legend(loc='lower left')
            ax.grid(True)
            ax.text(0.0,-50, f"MAG_APER[1]>50: {badMag}%", color='red')
            ax.xaxis.set_major_formatter(formatter)
            ax.set_xlabel("RA (deg)")
            ax.set_ylabel("Dec (deg)")
            ax.axis('off')

        for i, (catPath, subsetPath) in enumerate(zip(cat_chunk, subsetC_chunk)):
            catPlotTitle = os.path.basename(catPath)
            if 'cutout' in catPath:
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
            ax_cat.legend(loc='lower left')
            ax_cat.set_title(f"{catPlotTitle}", fontsize=10)
            ax_cat.xaxis.set_major_formatter(formatter)
            ax_cat.invert_xaxis()  # invert RA axis to match sky convention
            ax_cat.grid(True)

        plt.subplots_adjust(hspace=0.4, wspace=0.3)
        print("\n>>>>>>>> Saving detections figure... ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        if saveFig and overwrite:
            outDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/plots/'
            if any('cutout' in path for path in subset_chunk):
                outputPath = outDir + f"{batch_idx}_cutout_detections_badMags.png"
            else:
                outputPath = outDir + f"{batch_idx}_fullsize_detections_badMags.png"
            try:
                plt.savefig(outputPath, format='png', dpi=300, bbox_inches='tight')
                plt.close()  # closes current figure
                print(f"Figure saved as {outputPath}")
            except PermissionError:
                print(f"Error: Permission denied. Could not save the file {outputPath}.")
            except Exception as e:
                print(f"Error: Could not save the file {outputPath}. Reason: {e}")
        if show_fig:   
            print("\n>>>>>>>> Opening detections figure. ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            plt.show()

def bg_plotter(imagePaths='none', bgMapPaths='none', listBgSubDicts='none', badMags='none', segPaths='none', whtPaths='none', show_bad_mags=False, saveFig=True, overwrite=True, show_fig=False, verbose=True):
    """
        Makes figure showing data, bg mpa, bg subtracted, and seg map.
    """

    ### initialise fig
    ## file-save name, count number of cols needed 
    import math
    import matplotlib
    matplotlib.use("TkAgg")  # or "QtAgg"
    import matplotlib.pyplot as plt

    rows_per_fig = 3
    n_batches = math.ceil(len(listBgSubDicts) / rows_per_fig)

    print("\n>>>>>>>>> Making background comparison figures. For large images, this can take a few mins. Started at: ",
          datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    for batch_idx in range(n_batches):
        start = batch_idx * rows_per_fig
        end = start + rows_per_fig

        # slice chunks for this batch
        bgSub_chunk  = listBgSubDicts[start:end]
        image_chunk  = [] if imagePaths == 'none' else imagePaths[start:end]
        bgMap_chunk  = [] if bgMapPaths == 'none' else bgMapPaths[start:end]
        seg_chunk    = [] if segPaths == 'none' else segPaths[start:end]
        wht_chunk    = [] if whtPaths == 'none' else whtPaths[start:end]

        # determine number of rows in this figure
        nrows = len(bgSub_chunk)

        fig, axes = plt.subplots(nrows=nrows, ncols=5, figsize=(20, 4*nrows))
        axes = np.atleast_2d(axes)  # always 2D
        axes = axes.reshape(-1, 5)

        zscale = ZScaleInterval()
        formatter = ScalarFormatter(useMathText=False)
        formatter.set_scientific(False)
        formatter.set_useOffset(False)

        for i in range(nrows):  # one row per BgSubDict
            col_idx = 0

            # plot cutout region
            if imagePaths != 'none':
                imagePath = image_chunk[i]
                imagePlotTitle = os.path.basename(imagePath).replace('.fits', '')
                imageData, _, _ = get_fits_data(imagePath, verbose=verbose)
                vmin, vmax = zscale.get_limits(imageData)
                ax = axes[i, col_idx]
                ax.imshow(imageData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                ax.set_title(imagePlotTitle, fontsize=10)
                ax.xaxis.set_major_formatter(formatter)
                ax.axis('off')
                col_idx += 1

            # plot bgmap
            if bgMapPaths != 'none':
                bgMapPath = bgMap_chunk[i]
                bgMapPlotTitle = os.path.basename(bgMapPath).replace('.fits', '').replace('_cutout_', '')
                bgMapData, header, _ = get_fits_data(bgMapPath, verbose=verbose)
                vmin, vmax = zscale.get_limits(bgMapData)
                ax_bgmap = axes[i, col_idx]
                ax_bgmap.imshow(bgMapData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                ax_bgmap.set_title(bgMapPlotTitle, fontsize=10)
                ax_bgmap.xaxis.set_major_formatter(formatter)
                ax_bgmap.axis('off')
                col_idx += 1

            # plot bgSub
            if listBgSubDicts != 'none':
                BgSubDict = bgSub_chunk[i]
                bgSubPath = BgSubDict['Path']
                bgSubName = os.path.basename(bgSubPath).replace('.fits', '').replace('_cutout_', ' ')
                bgSubParams = f" BKSZ:{BgSubDict['back_size']} BKFILT:{BgSubDict['back_filtersize']}"
                bgSubPlotTitle = bgSubName + bgSubParams
                bgSubData, header, _ = get_fits_data(bgSubPath, verbose=verbose)
                vmin, vmax = zscale.get_limits(bgSubData)
                ax_bgsub = axes[i, col_idx]
                ax_bgsub.imshow(bgSubData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                ax_bgsub.set_title(bgSubPlotTitle, fontsize=10)
                ax_bgsub.xaxis.set_major_formatter(formatter)
                ax_bgsub.axis('off')
                col_idx += 1

            # plot segmentation maps
            if segPaths != 'none':
                segPath = seg_chunk[i]
                segPlotTitle = os.path.basename(segPath).replace('.fits', '').replace('_cutout_', '')
                segData, header, _ = get_fits_data(segPath, verbose=verbose)
                vmin, vmax = zscale.get_limits(segData)
                ax_seg = axes[i, col_idx]
                ax_seg.imshow(segData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                ax_seg.set_title(segPlotTitle, fontsize=10)
                ax_seg.xaxis.set_major_formatter(formatter)
                ax_seg.axis('off')
                col_idx += 1

            # plot weight maps
            if whtPaths != 'none':
                whtPath = wht_chunk[i]
                whtPlotTitle = os.path.basename(whtPath).replace('.fits', '').replace('_cutout_', '')
                whtData, header, _ = get_fits_data(whtPath, verbose=verbose)
                vmin, vmax = zscale.get_limits(whtData)
                ax_wht = axes[i, col_idx]
                ax_wht.imshow(whtData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
                ax_wht.set_title(whtPlotTitle, fontsize=10)
                ax_wht.xaxis.set_major_formatter(formatter)
                ax_wht.axis('off')
                col_idx += 1

    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    # --- save each batch as its own file ---
    if saveFig and overwrite:
        outDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/plots/'
        outputPath = outDir + f"{figTitleSubstr}_batch{batch_idx}.png"
        try:
            plt.savefig(outputPath, format='png', dpi=300)
            print(f"Figure saved as {outputPath}")
            plt.close()
        except PermissionError:
            print(f"Error: Permission denied. Could not save the file {outputPath}.")
        except Exception as e:
            print(f"Error: Could not save the file {outputPath}. Reason: {e}")

    print("\n>>>>>>>>> Making background comparison figure. For large images, this can take a few mins. Started at: ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    if saveFig and overwrite==True:
        outDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/plots/'
        outputPath = outDir +  figTitleSubstr + bgSubParams.replace(":", '').replace(' ', '_') + ".pdf"
        try:
            plt.savefig(outputPath, format='png', dpi=300)
            print(f"Figure saved as {outputPath}")
            plt.close(outputPath)
        except PermissionError:
            print(f"Error: Permission denied. Could not save the file {outputPath}.")
        except Exception as e:
            print(f"Error: Could not save the file {outputPath}. Reason: {e}")

    if show_fig: 
        print("\n>>>>>>>> Opening background figure. ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        plt.show()
       
######################## do photometry and measure depths ##############################################
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
            #print(origTable.colnames)
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
                #print('Renaming column from ', oldcolname, newcolname)
                bigTable.rename_column(oldcolname, newcolname)
            
    bigTable.info
    
    # remove the wht and seg columns
    #bigTable.remove_column('WHT_flux_' + apString)
    #bigTable.remove_column('SEG_flux_' + apString)

    if append:
        #print(bigTable.colnames)
        #print(origTable.colnames)
        
        # join with the big table!
        print('Appending to big aper phot table, lengths = {0}, {1}'.format(len(bigTable), len(origTable)))
        bigTable = join(origTable, bigTable, keys = ['IMAGE_xcenter', 'IMAGE_ycenter'])
        #print(bigTable.colnames)
        #print('After ', bigTable)
    #exit()
        
    #bigTable.meta['aperture_photometry_args'] = ''
    bigTable.write(outputFitsName, overwrite = True)
    #print("tired",bigTable)
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

def local_depths(cleanTable, apString, x, y, numApertures, zeropoint=-99.0, verbose=False, mask=False, sigmaClip=3.0, plot='none', regFile='none', fitGauss=False):
    ''' 
    Measures the local depth around a given aperture coordinate, or coodinates.
    Uses the closest numApertures apertures.

    Defines positions around which to make measurements > create .reg file if needed > do sigma clipping (twice) > fitr a gaussian to sigma-clipped distribution > normalise >> if zeropoint<-1.0, localDepths == sigma_MAD elif zeropoint>-1. call return_mag: i.e. if sigma_MAD is small, localDepth == -99 >> if there are no apertures nearby to x,y, it's a masked region. Record the position in a MaskArr as ==1 > check for NANs and set depth here to -99 > Return localDepths!

    cleanTable()        Reduced blank ap phot table which contain s only phot from  blank 
                        regions (according to seg map) which have a weight of greater than 0
    apString(str)       size of aperutre of interest as a string; for file naming purpose, etc.
    x(arr)              pixel coordinates at which to make depth measurement
    y(arr)              pixel coordinates at which to make depth measurement
    numApertures(int)   number of apertures you want to place around each object
    zeropint(float)     zeropoint for magnitude calculation
    mask(Boolean)       True for recorded masked regions to be returned (not near a good aperture)
    sigmaClip(float)    Number of times to perform sigma-clipping
    plot(str)           name of plot
    regFile(str)        Path to regions file containing known masked regions
    fitGauss(Boolean)   fit Gaussian to sigma-clipped distribution?
    '''
    
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
    
    ## do two lots of sigma clipping to improve accuracy
    
    ## define output array
    localDepths = np.zeros(x.size)
    maskArray = np.zeros(x.size)
    
    print("There are ", x.size, " positions to calculate depths for.")

    if regFile != 'none':
        tf = open(regFile, 'w')
        print('Masking with reg file, ', regFile)

    if plot != 'none':
        import matplotlib.backends.backend_pdf
        plotname = plot + '_check_mad.pdf' # Median Absolute Deviation 
        pdf = matplotlib.backends.backend_pdf.PdfPages(plotname)
    
    # get a minimum separation between depth measurements
    diffx = np.min(np.abs(y - np.roll(y, 1)))
    print('the miminum x separation is ', diffx)
    
    # loop through the positions
    for xi, xpos in enumerate(x):
        #for yi, ypos in enumerate(y):
        ypos = y[xi]

        ## calculate the radius        
        deltaX = apX - xpos
        deltaY = apY - ypos
        radius = np.sqrt(deltaX*deltaX + deltaY*deltaY)
        
        ## sort this array, and then the table
        #sortedIndicies = np.argsort(radius)
        ## do this faster!
        idx = np.argpartition(radius, numApertures)
        useIndicies = idx[0:numApertures] # use these indices for local measurements: meet min seperation requirement, uses correct step size
        
        ## Check that the radius here is close!!
        #print "The radius is ", radius[useIndicies]
        #print "the pos is ", xpos, ypos
        apRedX= apX[useIndicies]
        apRedY= apY[useIndicies]
        
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
        ## that will skew the StdDev
        medianFlux = np.median(apertureResults)
        mad = np.median(abs(apertureResults - medianFlux)) # Median Absolute Deviation
        sigma_mad = 1.4826*mad

        if (plot != 'none'):
            fig = plt.figure()
            bins = np.arange(-3.0*sigma_mad, 3.0*sigma_mad, sigma_mad/5.0)
            n,bin,patches = plt.hist(apertureResults, bins = bins, facecolor = 'green', alpha = 0.75)
            # split by region
            north = (apRedX > np.median(apRedX)) & (apRedY > np.median(apRedY))
            south = (apRedX < np.median(apRedX)) & (apRedY < np.median(apRedY))
            
            # plot the median etc
            plt.plot([medianFlux, medianFlux], [0, max(n)])
            plt.plot([medianFlux-sigma_mad, medianFlux-sigma_mad], [0, max(n)], color = 'k')
            plt.plot([medianFlux+sigma_mad, medianFlux+sigma_mad], [0, max(n)], color = 'k')
            plt.xlim([medianFlux -5.0*sigma_mad, medianFlux + 5.0*sigma_mad])

        if verbose:
            print("The mad sigma is ", sigma_mad, " after first run, mag = {0:.2f}".format(-2.5*np.log10(5.0*sigma_mad) + zeropoint))

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
                    if fitGauss: # fit a gaussian to these points!
                        from scipy.stats import norm
                        import matplotlib.mlab as mlab
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
        
        # give the user an update every 1000 positions about the progess of calculating the local depths     
        if xi % 1000 == 0:
            print("Calculating local depth at position ", xi, "out of ",x.size, " positions", datetime.now().strftime("%H:%M:%S"))

        if (plot != 'none'):
            pdf.savefig(fig)
            plt.close()
    
    if regFile != 'none':
        tf.close()

    if plot != 'none':
        pdf.close()
        print('Plot at ', plotname)

    # get rid of nans
    bad_indicies = np.isnan(localDepths)
    if np.any(bad_indicies):
        print('Fixing NANs in the local_depth code')
        localDepths[bad_indicies] = -99.0

    if mask: 
        return localDepths, maskArray
    else:
        print("CONSIDER updating code to use masked depths.")
        return localDepths

def return_mag(flux_counts, zeropoint, sigma = 1.0):
    """ If there is a very low number of flux counts, the magnitude is set to -99. 
         Otherwise, the magnitude is calculated as below. """
    import math
    if flux_counts < 1E-15:
        return -99.0
    else:
        return -2.5*math.log(sigma*flux_counts, 10.0) + zeropoint

def grid_depths(gridTable, x, y, faster=True, verbose=False, nearby=False, n_jobs=8):
    """
    Find the position of the closest depth measurement from a grid.

    Parameters
    ----------
    gridTable(astropy Tb)   Must have columns 'x', 'y', 'depths'.
    x, y(ndarray)           Coordinates of objects.
    faster(bool)            If True, uses grid-based binning.
                            If False, uses nearest-neighbour search.
    verbose(bool)           Print debug info.
    nearby(bool)            If True, averages over nearby pixels (numpoints=10).
    n_jobs(int)             Number of parallel jobs for joblib (-1 = all available cores).
    """

    xgrid = np.array(gridTable['x'])
    ygrid = np.array(gridTable['y'])
    depthsOverField = np.array(gridTable['depths'])

    depthArray = np.full(x.size, -99.0)

    if faster:
        if verbose:
            print("Using faster method.")
            print("Input array size is ", x.size)

        deltay = np.min(ygrid)
        deltax = np.min(xgrid)

        print("Doing grid_depths...", datetime.now().strftime("%d %H:%M:%S"))

        def process_gridpoint(xi):
            xmin = xgrid[xi] - deltax
            xmax = xgrid[xi] + deltax
            ymin = ygrid[xi] - deltay
            ymax = ygrid[xi] + deltay
            ii = (x > xmin) & (x <= xmax) & (y > ymin) & (y <= ymax)
            return ii, depthsOverField[xi]

        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(process_gridpoint)(xi) for xi in range(xgrid.size)
        )

        for ii, depth in results:
            depthArray[ii] = depth

    else:
        def process_object(xi):
            deltax = (xgrid - x[xi])
            deltay = (ygrid - y[xi])
            radius = np.sqrt(deltax*deltax + deltay*deltay)
            idx = np.argpartition(radius, 10)
            if nearby:
                mini = idx[:10]
            else:
                mini = [np.argmin(radius)]
            return xi, depthsOverField[mini][0]

        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(process_object)(xi) for xi in range(x.size)
        )

        for xi, depth in results:
            depthArray[xi] = depth

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
  
       print("The delta is ", deltax, deltay)
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
               #print("The nearby PSFs are = ", PSFsOverField[mini])
               #print("Before = ", PSFsOverField[mini][0])
           
           
       #print "The closest point is ", xgrid[mini], ygrid[mini], " to ", x[xi], y[xi]
           depthArray[xi] = PSFsOverField[mini][0]
       #exit()
   
   return psfArray


def extract_local_depths(inputTableFile, apDiametersAS, zeropoint, step=500, numApertures=200, strips=False, plot=False, local=True, recalculate=False, globalplot=False, clean=True, plotDir='', maskreg='none', refimage='none'):
    ''' extractDepths.py
    
    For each aperture diameter supplied in apDiametersAS, 
    extract the local depths from the blank aperture photometry table (inputTableFile) by first defining the suitable blank-sky positions to measeure the depths in by requiring the seg map to be empty and the weight at the position to be greater than zero. If the local depths have not already been calculated and stored in localDepthsFile, call local_depths() to make it (and optionally the maskArr). If localDepthsFile already exists, use this. 

    Remove any negative values from localDepthsFile, finally defining the local depths table for the image. Next, find the local depths for any regions of interest e.g. stripes. Call return_instrips() to do this. If maskreg does not exist/is not supplied, crudely buffer around the edges of the image where depth is expected to be lower (VIRCAM/VISTA User Manual, 2009, Fig 17). Define the local depths tables for these regions of interest. Make plots. 

    GLOBAL DEPTHS

    INPUT(s)
        inputTableFile(str)     Path to blank aperture photometry file i.e. phot/H_aperPhot.fits
        apDiametersAS(arr)      Array of aperture diameters you want to do phot with i.e. 1.8, 2.0, 3.0
        zeropoint
        step(int)
        numApertures(int)       number of apertures you want to place around each object
        strips(Boolean)         True for UVISTA. Observiing stategy means stripes must be accounted for 
        plot(Boolean)           Do you want to make a plot of ...
        local(Boolean)          ?
        recalculate(Boolean)    Do you want to overwrite the depths/phot/localDepthsFile ?
        globalplot(Boolean)     ?
        clean(Boolean?          ?
        plotDir(str)            Where do you want to save optional plot?
        maskreg(str)            ? directory to mask.reg file
        refimage(str)           ? reference image

    OUTPUT(s)
        Writes localDepthsFile and optionally MaskArr

        regions, global_depth, median_local_depth, mode_local_depth
    '''

    import matplotlib
    #matplotlib.use('pdf')
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    from matplotlib import gridspec
    from scipy.stats import norm
    from astropy import units as u
    import matplotlib.backends.backend_pdf
    from hu_catalogue_codes import return_instrips, mask_column
    
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
        #deepStrips_ra_high = [149.5, 149.85, 150.25, 150.6] #TODO - did RB have the same issue?/think bout the bad mag patches?
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

    # create empty arrays for results - "bad" measurements are set to 0 here
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
            
        ## Get the x and y coordianates of the apertures
        apX = reducedTable['IMAGE_xcenter'] #/u.pix pixel coords
        apY = reducedTable['IMAGE_ycenter'] #/u.pix
        
        ## Check if I need to run local depths from the apertures
        
        if not os.path.isfile(localDepthsFile) or recalculate:
            ## Files doesn't exist or I want to recalculate it  
            
            #############################################################
            ################## RUNNING LOCAL DEPTHS ####################
            
            ## Create array of positions at which to measure the local depths based on the xy coords in the reduced blankapPhot table
            xmax = max(apX)
            ymax = max(apY)
            
            numX = np.ceil(xmax/step) # number of measurements per dimension i.e. grid
            numY = np.ceil(ymax/step)

            x = min(apX) + np.arange(numX)*step ## modifed 21/9/2018 to add min(apX).
            y = min(apY) + np.arange(numY)*step
            #print("Step = ", step, " numx, y = ", numX, numY)
            #print("Max = ", xmax, ymax, max(apX), min(apX))
            
            # create x, y arrays
            x = np.zeros(1) #.value
            y = np.zeros(1) #.value                
            
            #populate the zero arrays
            for xi in np.arange(step/2.0, numX*step, step):
                for yi in np.arange(step/2.0, numY*step, step):
                    x = np.append(x, xi)
                    y = np.append(y, yi)
                    
            # I want a constant grid over the image.
            
            # remove the first elements
            x = x[1:]
            y = y[1:]

            ## Now run local depths at those points

            depthsLocalFull, maskArray = local_depths(reducedTable, apString, x, y, numApertures, zeropoint=zeropoint, mask=True, sigmaClip=3.0)#, plot = plotDir + filterName + '_' + str(numApertures))
            
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
            #print("Calculating depths in region ", region)
            if region != 'full':
                #print('Splitting by strip')
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
                    magmax = 27.0 #TODO why hard-coded?
                    
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
                #print('Bins defined here ', bins)
                
                if strips:
                    ## plot subsets of the results
                    ## and colour based on strip/gap 
                    #print('In strips', bins.shape)
                    
                    ## 1) get the strips
                    ss = return_instrips(apXregion, apYregion)
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
                    #print('Definiting msoothed here, bins ', smoothed.shape, bins.shape)
                    
                else:
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
    #for x, xi in enumerate(regions): print('{0}, {1:.2f}'.format(xi, global_depth[x,0]))
    
    return regions, global_depth, median_local_depth, mode_local_depth

def image_depth(imagePath, zeropoint, cutouts=[], size='none', back_size=32, back_filtersize=9, apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtPath='NONE', whtType='NONE', IRACapDiametersAS=np.array([2.8, 3.8, 5.8, 9.8, 11.6]), segPath='NONE', outputDir='none', filterName='NONE', numApertures=300, step=200, overwrite=False, inputSex=baseDir+'data/bertin_config/video_mine.sex', strips=False, bgSub=True, mask='none', gridSepAS=3.0):
    """
    Perform background modelling and subtraction using SourceExtractor. If blank-sky aperture photometry results file doesn't exist, call aperture_photometry_blank() to estimate the blank-sky noise and write the file. Call extract_local_ddepths() to give local and global depths. Returns the paths to the bg-map and bg-subbed images and seg map.

    """

    warnings_triggered = 0
    cd1Test = False

    # with 300 apertures, the radius is roughly 200 due to blending etc
    # hence a step of 200 means that these don't overlap
    
    # Tinsley SE env
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/share/sextractor'

    # define the output files
    if outputDir == 'none':
        outputDir = baseDir
        if os.path.isdir(outputDir) == False:
            os.system('mkdir ' + outputDir)
        print("outputDir will be: ", outputDir)

    # seg map, bkg map, bgsub images
    imageDir = outputDir + 'images/'
    if os.path.isdir(imageDir) == False:
        breakpoint()
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

    if cutouts != 'none':
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

        if cutouts != 'none': 
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
            print("Running image_depths() SE at: ", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            print(command)
            os.system(command)

        else:
            print(f"The SEG and/or BG subtracted map exist at: {segPath, bgSubPath} \n")
    '''
    aperPhotFile = aperDir + filterName + '_aperPhot.fits' # blank aperture phot measuremnts
    overwrite = False
    if os.path.isfile(aperPhotFile) == False or overwrite == True:
        # define the grid separation
        gridSepPixels = np.ceil(gridSepAS/pixScale) # 5'' separation
        #print("test: ",gridSepPixels)
        # gridSepPixels = 10.0
        
        # make this tunable...
        if bgSub == False:
            bgSubPath = imagePath
            print("Not using bg subtracted image.")
        print(">>>>>>> Measuring the aperture photometry...", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        ii = segPath.rfind('NIRSPEC')
        if ii > 0:
            field = 'NIRSPEC'
        else:
            field = 'NONE'

        # if aperphotfile doesn't exist, call ap phot blank - outputFitsName is aperPhotFile
        aperture_photometry_blank(bgSubPath, segPath, whtPath, apDiametersAS, gridSeparation = gridSepPixels, clean = True, outputFitsName = aperPhotFile, imageDir = imageDir, field = field, overwrite = overwrite)

    #######################################################################
    # Then calculate the local depths, and make a nice plot
    # if COSMOS, I need to run in strips too
    recalculate = False

    # mask
    regions, globaldepths, meddepths, modedepths = extract_local_depths(aperPhotFile, apDiametersAS, zeropoint, recalculate=recalculate, numApertures=numApertures, step=step, plotDir=plotDir, strips=strips, maskreg=mask, refimage=bgSubPath) #, plot = True)
    
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
    print("Output file saved to ", depthFile)'''
    
    return bgMapPath, bgSubDict, segPath

def get_depths(fieldName, fullsizeimages='none', cutouts='none', size='none', back_size=32, back_filtersize=9, queue='none', reqFilters=['all'], apDiametersAS=np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir=baseDir+'data/', outputDir='none', overwrite=False, ra_str='none', dec_str='none', verbose=True, saveFig=True):
    """
    Calls necessary functions to measure and save the local and global depths for given filters within a supplied .lis file. Opens the detection catalogues for each filter (and a "bad" - MAG_APER[1]>50 - subset of each catalogue) in TOPCAT. Shows and saves comparitive images of the background models etc. Shows and saves the detections and the bad subset on top of the science image.
    """

    # set the grid seperation in arcsec
    if fieldName == 'NIRSPEC':
        print('Nirspec data, setting grid separation to be much lower!')
        gridSepAS = 0.5
    else:
        gridSepAS = 3.0
        
    # Read in the images file
    dirHere = dataDir 
    print("dirHere = dataDir ", dataDir )
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
        
        if cutouts != 'none':
            for cutout in set(cutouts): # for unique items in list
                cutoutPath = cutout
                cutoutName = os.path.basename(cutoutPath)
                imageName = imageName.split('.')[0]
                #if cutoutName.startswith(imageName) and ('wht' not in cutoutName):
                if cutoutName.startswith(imageName):
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
            imageDir = dataDir
        
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

            # make shell script
            tmpName = "tmp_{1}_{0}.sh".format(filterName, fieldName)
            f = open(tmpName, 'w')
            f.write('#!/bin/bash\n')
            f.write('python3 stupid.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(imageDir + imageName, imageDir +  whtName, whtType, zeropoint, outputDir, strips, filterName, overwrite, maskName, gridSepAS, apDiametersASstring))
            f.close()
            
            # now execute this
            command = "addqueue -c 'tmp_{0}' -m 9 -q {0} -d ./{1}".format(queue, tmpName)
            print(command)
            os.system('chmod +x {0}'.format(tmpName))
            os.system(command)

    # after each filter has been SE'd, open the cats
    catPaths = []
    if outputDir == 'none':
        outputDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/catalogues/'

        if cutouts != 'none':
            for catName in os.listdir(outputDir):
                for filt in reqFilters:
                    pattern = (rf"d{size}{filt}_cutout\.fits")
                    if re.fullmatch(pattern, catName):
                        catPath = os.path.join(outputDir, catName)
                        catPaths.append(catPath)
        elif fullsizeimages != 'none':
            for catName in os.listdir(outputDir):
                for filt in reqFilters:
                    pattern = (rf"d{filt}\.fits")
                    if re.fullmatch(pattern, catName):
                        catPath = os.path.join(outputDir, catName)
                        catPaths.append(catPath)
        else:
            print("No bingo!")

    subsetPaths = open_cats(catPaths, open_subset=True, overwrite=True)

    cutoutPaths = [] # selection of files that will be passed to detections_fig() (i.e. excl wht files)
    whtPaths = []    
    badMags = [] # all detections/bad_mag measurements

    for subsetPath in subsetPaths:
        if cutouts != 'none' and fullsizeimages == 'none':
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

        elif cutouts == 'none' and fullsizeimages != 'none':
            for filt in reqFilters:
                subpattern = (rf"{baseDir}depths/catalogues/"
                            rf"d{filt}_subset\.fits")
                if re.fullmatch(subpattern, subsetPath):
                    catPath = f"{baseDir}depths/catalogues/d{filt}.fits"
                    imagePath = f"{baseDir}data/COSMOS/UVISTA_{filt}_DR6.fits" # corresponds to the catalogue and subset catalogue
                    HSCimagePath = f"{baseDir}data/COSMOS/{filt}_DR3.fits" 
                    catTable = Table.read(catPath)
                    subsetTable = Table.read(subsetPath)
                    badMag = round(100*(len(subsetTable)/len(catTable)),5)
                    badMags.append(badMag)
                    print(f"{filt}-bad-proportion: ", badMag, "%")

    for subsetPath in subsetPaths:
        for filt in reqFilters:
            subpattern = f"{baseDir}depths/catalogues/d{size}{filt}_cutout_subset.fits"
            fullsize_subpattern = f"{baseDir}depths/catalogues/d{filt}_subset.fits"
            if subpattern == subsetPath:
                catPath = f"{baseDir}depths/catalogues/d{size}{filt}_cutout.fits"
                whtPath = f"{baseDir}data/COSMOS/cutouts/UVISTA_{filt}_DR6_wht_{ra_str}_{dec_str}_size{size}.fits" # corresponds to the catalogue and subset catalogue
                HSCwhtPath = f"{baseDir}data/COSMOS/cutouts/{filt}_DR3_wht_{ra_str}_{dec_str}_size{size}.fits" 
                if os.path.isfile(whtPath):
                    whtPaths.append(whtPath)
                elif os.path.isfile(HSCwhtPath):
                    whtPaths.append(HSCwhtPath)
            elif fullsize_subpattern == subsetPath:
                catPath = f"{baseDir}depths/catalogues/d{filt}.fits"
                whtPath = f"{baseDir}data/COSMOS/UVISTA_{filt}_DR6_wht.fits" # corresponds to the catalogue and subset catalogue
                HSCwhtPath = f"{baseDir}data/COSMOS/{filt}_DR3_wht.fits"
                if os.path.isfile(whtPath):
                    whtPaths.append(whtPath)
                elif os.path.isfile(HSCwhtPath):
                    whtPaths.append(HSCwhtPath)

    if cutouts != 'none':
        #detections_fig(cutoutPaths=cutoutPaths, catPaths=catPaths, subsetPaths=subsetPaths, badMags=badMags, saveFig=saveFig, overwrite=saveFig, verbose=verbose)

        bg_plotter(cutoutPaths=cutoutPaths, bgMapPaths=bgMapPaths, listBgSubDicts=listBgSubDicts, badMags=badMags, segPaths=segPaths, whtPaths=whtPaths, saveFig=saveFig, overwrite=saveFig, verbose=verbose)

    #else:
        #detections_fig(fullsizeimages, catPaths=catPaths, subsetPaths=subsetPaths, badMags=badMags, saveFig=saveFig, overwrite=saveFig, verbose=verbose)

        #bg_plotter(fullsizeimages, bgMapPaths=bgMapPaths, listBgSubDicts=listBgSubDicts, badMags=badMags, segPaths=segPaths, whtPaths=whtPaths, saveFig=saveFig, overwrite=saveFig, verbose=verbose)
        

    return
    
