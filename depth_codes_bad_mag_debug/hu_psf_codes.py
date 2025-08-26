import numpy as np
import os
from astropy.table import Table,Column, hstack

baseDir = '/raid/scratch/'
baseDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/'# TODO change back

def psfex(imageName, filterName, fieldName, zeropoint, depthDir, apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), whtName = 'NONE', whtType = 'NONE', segName = 'NONE', outputDir = 'none', overwrite = False, overwritePSF = False, inputSex = baseDir + 'data/bertin_config/video_psfex.sex', inputPSFEx = baseDir + 'data/bertin_config/default_cube.psfex', strips = False, bzk = False, starsOnly = False, useCat = 'NONE'):
    print("Running psfex(): ", '\n')
# Note: changed to my config_files dir in data/HSC_SSP_DR3

    from astropy.io import fits
    import sep
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import matplotlib.backends.backend_pdf

    # updated for tinsley
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/share/sextractor'
    
    plot = True
    
    ##############################################
    psfSize = 75
    psfDegree = 5
    psfNSNAP = 10

    assocRad = 1.0/0.15 

    rstepAS = 0.05 ## arcsec
    rmaxAS = 9.0 ## arcsec
    numRad = int(rmaxAS/rstepAS)
    rAS = np.arange(numRad)*rstepAS + rstepAS
    normRas = 5.0
    
    ##############################################
    
    # main code where all the action happens!
    # first define the output files
    if outputDir == 'none':
        # put everything here!
        # make a new directory
        outputDir = 'psf/'
        if os.path.isdir(outputDir) == False:
            os.system('mkdir ' + outputDir)

    # for the seg map
    #imageDir = outputDir + 'images/'
    #if os.path.isdir(imageDir) == False:
    #    os.system('mkdir ' + imageDir)

    plotDir = outputDir + 'plots/'
    if os.path.isdir(plotDir) == False:
        os.system('mkdir ' + plotDir)

    catDir = outputDir + 'catalogues/'
    if os.path.isdir(catDir) == False:
        os.system('mkdir ' + catDir)

    efDir = outputDir + 'enclosedflux/'
    if os.path.isdir(efDir) == False:
        os.system('mkdir ' + efDir)

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
        
    # check all the necessary files exist
    imyes = os.path.isfile(imageName)
    if imyes == False:
        print("ERROR: the image file does not exist, check the path:\n" + imageName)
        exit()
        
    if whtName != 'NONE':
        whtyes = os.path.isfile(whtName)
        if whtyes == False:
            print("ERROR: the weight file does not exist, check the path:\n" + whtName)            
            exit()

    
    ####################################################
    # Get a star catalogue
    if filterName == 'NONE':
        filterName = baseName

    # get the pixel scale

    hdulist = fits.open(imageName)
    imageHeader = hdulist[0].header
    cdone_o = -3600.0*imageHeader['CD1_1']
    pixScale = round(cdone_o, 5)
    print("The pixel scale is {0:.4f}".format(pixScale))
    naxis1 = imageHeader['NAXIS1']
    naxis2 = imageHeader['NAXIS2']
    
    rpix = rAS/pixScale

    #####################################################
    # Define the point sources for the run
    depthCat = depthDir + '{0}/catalogues/d{1}.fits'.format(fieldName, filterName)
    # get the five sigma limit
    fiveSigFile = depthDir + '{0}/results/{1}_300.txt'.format(fieldName, filterName)
    fiveSig = Table.read(fiveSigFile, format = 'ascii.commented_header')
    # HOLLY 10/03/25 10:23 - ch1cds/ch2cds don't have 2.0as; 2.8as is next closest size
    if filterName == "ch1cds" or filterName == "ch2cds":
        fiveSigHere = fiveSig['2.8as'][0]
    else:
        fiveSigHere = fiveSig['2.0as'][0]
    
    if bzk:
        print('BZK stars...')
    else:
        # use the depth catalogue
        print('Using catalogue', depthCat)
        plotName = plotDir + '{0}_sizemag.pdf'.format(filterName)
        starName = sizemag_stars(depthCat, fiveSigHere, filterName, fieldName, pixScale, plotDir = plotDir, starDir = catDir, plotName = plotName, sizeType = 'FLUX_RADIUS')
        #starName = catDir + '{0}_stars.ascii'.format(filterName)
    #print(starName)

    if starsOnly:
        return

    if useCat != 'NONE':
        starName = starName.replace(filterName, useCat)
        print('Using star cat name ', starName)

    #############################################################
    # Run a catalogue in the required format for PSFEx
    seCatalogue = catDir  + filterName + '.fits'
    print("seCatalogue: ", seCatalogue)
    
    # FITS_LDAC
    keywords = ' -CATALOG_TYPE FITS_LDAC -CATALOG_NAME ' + seCatalogue + \
               ' -MAG_ZEROPOINT '+ str(zeropoint)
    
    keywordsWeight = ' -WEIGHT_TYPE ' + whtType + \
                     ' -WEIGHT_IMAGE ' +  whtName
    
    #starCatalogue = 'sky.list'
    associateKeywords = ' -ASSOC_DATA 1,2,3 -ASSOC_NAME ' + starName + \
        ' -ASSOC_PARAMS 4,5 -ASSOC_RADIUS ' + str(assocRad) + \
        ' -ASSOC_TYPE NEAREST' + \
        ' -ASSOCSELEC_TYPE MATCHED'
    
    ## for Spitzer I need to alter this a bit _LDAC
    #if filterName == 'ch1' or filterName == 'ch2':
    #    extraKeywords = ' -SATUR_KEY DUMMY -SATUR_LEVEL 35.0 -SEEING_FWHM 2.0 -PHOT_APERTURES ' + str(10.0/pixScale) + ' -DETECT_MINAREA 10 -DETECT_THRESH 200 -ANALYSIS_THRESH 10 -DEBLEND_MINCONT 0.01 -STARNNW_NAME none'
    #extraKeywords = ' -SEEING_FWHM 2.0 -PHOT_APERTURES ' + str(12.0/pixScale) + ' -DETECT_MINAREA 10 -DETECT_THRESH 200 -ANALYSIS_THRESH 10 -DEBLEND_MINCONT 0.01'
    #else:
    #    extraKeywords = ''
        
    extraKeywords = associateKeywords
          
    if os.path.isfile(seCatalogue) and overwrite == False:
        print("Catalogue has been created previous.  Skipping this step.")
    else:
        print('sex '+  imageName +' -c ' + inputSex + extraKeywords + keywords + keywordsWeight )
        os.system('/usr/local/sextractor/bin/sex '+  imageName +' -c ' + inputSex + extraKeywords + keywords + keywordsWeight) # Note from Rohan: need to run SExtractor differently.
        print("Source extractor catalogue has been saved to ", seCatalogue)
        
        ## Read this catalogue in?
        #yy = Table.read(seCatalogue, format = 'fits', hdu = 'LDAC_IMHEAD')
        #xx = Table.read(seCatalogue, format = 'fits', hdu = 'LDAC_OBJECTS')
        if (filterName[0:2] == 'ch'): #fieldName == 'COSMOS' or
            hdulist = fits.open(seCatalogue)
            datahere = hdulist[2].data
            datahere['FLUX_RADIUS'] = datahere['FLUX_RADIUS']/3.0
            #print datahere['FLUX_RADIUS'].size
            #goodi = np.where(datahere['FLUX_RADIUS'] < 3.4)
            #datahere.remove_row(10)# = datahere[goodi]
            
            hdulist.writeto(seCatalogue, clobber = True)
            print("Editing the se catalogue for my ch results.")
            
    #############################################################
    # Run PSFEx on this catalogue!
    
    ## keywords
    keywords = ' -PSF_DIR ' + resultsDir + ' -PSF_SIZE ' + str(psfSize) + \
               ',' + str(psfSize) + ' -PSFVAR_DEGREES ' + str(psfDegree) + \
               ' -PSFVAR_NSNAP ' + str(psfNSNAP) + ' -CHECKIMAGE_TYPE ' + \
               ' SNAPSHOTS,SAMPLES,RESIDUALS ' + ' -CHECKIMAGE_NAME ' + \
               '{0}snap.fits,{0}samples.fits,{0}resi.fits '.format(resultsDir) + \
               ' -SAMPLE_AUTOSELECT N -PSF_SAMPLING 1.0 -SAMPLE_FWHMRANGE 1,10' + \
               ' -CHECKPLOT_TYPE FWHM -CHECKPLOT_DEV PSC' + \
               ' -CHECKPLOT_NAME {0}fwhm'.format(plotDir) # -BASIS_TYPE GAUSS-LAGUERRE -SAMPLE_FWHMRANGE 2.0,8.0 -VERBOSE_TYPE FULL'.format(plotDir, filterName) #,{0}ellipticity,{0}residuals,{0}chis'.format(resultsDir) + \
               #'-SAMPLE_FWHMRANGE 2.0,8.0 ' -PSF_ACCURACY 0.1'#+ ' -CHECKIMAGE_TYPE ' #,ELLIPTICITY,RESIDUALS,CHI2 ' + \
               #' PROTOTYPES,SNAPSHOTS,RESIDUALS ' + ' -CHECKIMAGE_NAME ' + \
               #'proto.fits,snap.fits,resi.fits'
    
    # Define the important psf file
    psfOutFile = resultsDir + 'snap_{0}.fits'.format(filterName)
    
    # Note from Rohan: I want to try some things, so will create different snap files.
    #psfOutFile = resultsDir + 'snap_{0}_rohan.fits'.format(filterName)
    
    print('psfex '+ seCatalogue +' -c ' + inputPSFEx + keywords)
    overwritePSF = False
    
    if os.path.isfile(psfOutFile) and overwritePSF == False:
        print("PSFEx has already been run, do not overwrite.")
    else:
        print("Running psfex", seCatalogue)
        os.system('psfex '+ seCatalogue +' -c ' + inputPSFEx + keywords)    

    #exit()
        
    #############################################################
    # Extract the PSFs over the full FOV
    print("Now extracting the PSF snapshots over the field. ")
    hdulist = fits.open(psfOutFile)
    hdulist.info()
    fullPSFdata = hdulist[0].data
    fullPSFheader = hdulist[0].header
    shapePSF = fullPSFdata.shape
    print("The PSF dimensions are ", shapePSF)
    
    psf_samp = fullPSFheader['PSF_SAMP']
    
    halfwidth = (shapePSF[2]-1)/2
    print("The central pixel is (", halfwidth, ',', halfwidth, ')')
    
    stepsizex = naxis1/psfNSNAP
    stepsizey = naxis2/psfNSNAP
    
    print("There are ", psfNSNAP*psfNSNAP, " psfs to extract.")
    #psfNSNAP = 1
    
    ## Make an output array
    enclosedFluxArray = np.zeros((len(apDiametersAS), psfNSNAP*psfNSNAP))
    x = np.zeros(psfNSNAP*psfNSNAP)
    y = np.zeros(psfNSNAP*psfNSNAP)

    for xi in range(psfNSNAP):
        for yi in range(psfNSNAP):
            
            posx = stepsizex*(xi+0.5)
            posy = stepsizey*(yi+0.5)
            index = yi + xi*psfNSNAP
            
            ## This is correct!!  Checked with UltarVISTA - looking at the corners
            ## not covered by HSC
            x[index] = naxis1-posx
            y[index] = naxis2 -posy
            
            #print "Index = ", index
            
            ## Fix the byte problem

            data = fullPSFdata[xi, yi, :, :]
            psfHereFull = data.byteswap().newbyteorder()    
            psfHere = psfHereFull/np.sum(psfHereFull)

            ## plot the results
            if plot:
                if (index > 5) & (index < 7):
                    cogplot = resultsDir + '{0}_COG.pdf'.format(filterName)
                    fig, ax1 = plt.subplots()    
            
            xArray = halfwidth
            yArray = halfwidth
            
            #xArray = np.array([halfwidth, halfwidth+0.5, halfwidth+0.5, halfwidth-0.5, halfwidth-0.5, halfwidth+1.0, halfwidth+1.0, halfwidth-1.0, halfwidth-1.0]) 
            #yArray = np.array([halfwidth, halfwidth+0.5, halfwidth-0.5, halfwidth-0.5, halfwidth+0.5, halfwidth-1.0, halfwidth+1.0, halfwidth-1.0, halfwidth+1.0]) 
            
            ##############################################
            ## Run my enclosed flux code on this cut-out
            
            ## IMPORTANT - I need to scale by PSF_SAMP
            rpixHere = rpix/psf_samp
            #print(len(psfHere), len(xArray), len(yArray), rlen(pixHere))
            flux, fluxerr, flag = sep.sum_circle(psfHere, [xArray], [yArray], rpixHere, subpix = 5)

            # normalise this
            jj = (rAS == normRas)
            fluxnorm = flux/flux[jj]
            
            ## Plot the results
            if plot:
                if (index > 5) & (index < 7):
                    line1, = plt.plot(rAS, flux, label = 'Global bg', linestyle = '--', color = 'deepskyblue', )


                    plt.plot(rAS, fluxnorm, linestyle = ':', color = 'red')
                    
                    ax1.set_ylim([0.4, 1.1])
                    plt.plot([0, 10], [1.0, 1.0], 'k:', linewidth = 2)
                    
                    plt.xlabel('Aperture Radius (arsec)')
                    plt.title('Enclosed flux for filter ' + filterName)
                    
                    ## Add an image
                    left, bottom, width, height = [0.55, 0.3, 0.3, 0.3]
                    ax2 = fig.add_axes([left, bottom, width, height])
                    ax2.imshow(psfHere, cmap = 'GnBu_r', origin = 'lower', interpolation = 'none')
                
            ## Find the values of interest
            for ri, diameterRequired in enumerate(apDiametersAS):

                radiusRequired = diameterRequired/2.0
                cond = (rAS > radiusRequired-0.0001) & (rAS < radiusRequired+0.0001)
                fluxEnclosed = np.extract(cond, flux)
                fluxEnclosedNorm = np.extract(cond, fluxnorm)
                
                #print  "Aperture = ", radiusRequired, " AS, fenc = ", fluxEnclosed
                enclosedFluxArray[ri, index] = fluxEnclosed[0]
                if plot:
                    if (index > 5) & (index < 7):
                        ax1.text(2.0, 0.5+ri*0.05, '{0:.2f} ({2:.2f}), {1:.1f}as'.format(fluxEnclosed[0], radiusRequired, fluxEnclosedNorm[0]))
                
            if plot:
                if (index > 5) & (index < 7):
                    fig.savefig(cogplot)
                    print(cogplot)
                               
    ############################################################
    ## Tabulate these results, so they can be extracted
    ## at will e.g. in my flux_errors code
    localTable = Table([x, y], names = ['x', 'y'], dtype = ['f4', 'f4'])
    
    ##############################################################
    ## Make plots of the enclosed flux and check these make sense
    ## Create pdf with multiple pages
    plotName = plotDir + filterName + '.pdf'
    pdf = matplotlib.backends.backend_pdf.PdfPages(plotName)
    print("making enclosed flux table")
    outputTable = efDir + filterName + '.fits'

    # make a txt file of the median ef values
    medianFile = efDir + filterName + '_peak.txt'
    efmedian = np.zeros(apDiametersAS.size)
    
    for ri, diameterRequired in enumerate(apDiametersAS):
        
    ## plot the results
        fig = plt.figure(figsize=(6,8))
        gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
        plt.subplot(gs[0])
        plt.axis('equal')
    ## get colours
        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(x, y, s = 30, c = enclosedFluxArray[ri, :], cmap = cm, \
                             linewidth = 0.0) #, vmin = magmin, vmax = magmax)
        plt.colorbar(sc)
        apDstring = "{0:.1f}".format(diameterRequired)
        plt.title('Enclosed flux for filter = ' +filterName + " apD = " + apDstring)
        
        ## now make a histogram to go underneath
        ax = plt.subplot(gs[1])
        binwidth = 0.01
        efhere = enclosedFluxArray[ri,:]
        median = np.median(efhere)
        mad = 1.58*np.median(np.abs(efhere-median))
        #print('med = {0:.2f}, mad = {1:.2f}'.format(median, mad))
        #trimmed = (efhere < 1.0) & (efhere > 0.0)
        #minef = np.round(np.min(efhere)*100.0)/100.0
        #maxef = np.max(efhere)
        minef = median - 3.0*mad
        maxef = median + 3.0*mad

        if (maxef - minef) < 0.05:
            minef = minef-0.02
            maxef = maxef+ 0.02

        print('Min = {0}, max = {1}'.format(minef, maxef))
        bins = np.arange(minef+binwidth/2.0, maxef, binwidth)
        histy, histx, _ = plt.hist(enclosedFluxArray[ri, :], facecolor= 'k', alpha = 0.8, bins = bins,  density = True) #, bins = 100) #, range = [0.0, 1.0])#, \
                                       # bins = np.arange(magmin, magmax+binwidth, binwidth))
        
        ## calculate the median
#        print(histy, bins)
        smallbins = bins[:-1]
        print(histy)
        mm = (histy == np.max(histy))
 #       print('The max is at ', smallbins[mm]+binwidth/2.0)
        enclosedFlux = smallbins[mm]+binwidth/2.0
        #        enclosedFlux = np.median(enclosedFluxArray[ri, :])
        #        enclosedFlux = np.median(efhere[trimmed])
        efmedian[ri] = enclosedFlux[0]
        
        plt.plot([enclosedFlux, enclosedFlux], ax.get_ylim())
        
        pdf.savefig(fig)
        
        ## add the enclosed flux values to the table
        newcolumn = Column(name = 'ef_' + apDstring, data = enclosedFluxArray[ri, :], dtype = 'f4')
        localTable.add_column(newcolumn)     
        
    pdf.close()
    plt.close()
    print("Plot saved to ", plotName)

    # Write out the median
    newtb = Table([apDiametersAS, efmedian], names = ['apD', 'ef'])
    newtb['ef'].format = '7.2f'
    newtb.write(medianFile, overwrite = True, format = 'ascii.commented_header')
    print(medianFile)
    
    # Write out table
    localTable.write(outputTable, overwrite=True)
    print("Table saved to ", outputTable)

    return

def get_psf(fieldName, queue = 'none', reqFilters= ['all'], apDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0]), dataDir = baseDir +'data/', outputDir = 'none', overwrite = False, starsOnly = False, useCat = 'NONE'):
    
    print("Running get_psf(): ", '\n')
    from hu_depth_codes import read_image_lis
    
    # Read in the images file
    dirHere = dataDir + fieldName + '/'
    imagedata = read_image_lis(dirHere)
    print('\n', '\n', "imagedata: ", imagedata['directory'], '\n', '\n')
    availableFilters = np.array(imagedata['Name'])
    print("The available filters are ", availableFilters)

    bzk = False
    
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
    stripFilt = np.array(['Y', 'J', 'H', 'Ks', 'JH', 'YJ', 'HKs']) # 'HKs' changed to 'HK' 22/05/25

    # loop through each filter
    # for the queue run, run each as a separate file...

    for fi, filterName in enumerate(availableFilters):
        # define the images etc to send through
        imageName = imagedata['Image'][fi]
        print('\n', '\n', "imageName: ", imageName, '\n', '\n')
        whtName = imagedata['Weight'][fi]
        whtType = imagedata['Wht_type'][fi]
        zeropoint = imagedata['zeropoint'][fi]
        imageDir = imagedata['directory'][fi]
        print('\n', '\n', "imageDir: ", imageDir, '\n', '\n')


        #depthCat = '{0}depths/{1}/catalogues/d{2}.fits'.format(dataDir, fieldName, filterName)
        depthDir = dataDir + '../depths/'
        print("test: ",depthDir, os.path.isdir(depthDir))

        strips = "False"
        if fieldName == 'COSMOS':
            jj = (filterName == stripFilt)
            if np.any(jj):
                strips = "True"
        
        if imageDir == 'here':
            imageDir = dataDir + fieldName + '/'
        
        # Now spawn the depths!
        if queue == 'none':
            print("Running here ")
            
            psfex(imageDir + imageName, filterName, fieldName, zeropoint, depthDir, whtName = imageDir + whtName, whtType = whtType, outputDir = outputDir, strips = strips, overwrite = overwrite, starsOnly = starsOnly, overwritePSF = False, useCat= useCat)

        else:
            print("Spawning in the queue...", queue)
            # make shell script
            tmpName = "tmp_{1}_{0}.sh".format(filterName, fieldName)
            f = open(tmpName, 'w')
            f.write('#!/bin/bash\n')
            f.write('python3 stupid_psf.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(imageDir + imageName, imageDir +  whtName, whtType, zeropoint, outputDir, strips, filterName, overwrite, starsOnly, fieldName, depthDir))
            f.close()
            
            # now execute this
            command = "addqueue -c 'tmp_{0}' -m 8 -q {0} -d ./{1}".format(queue, tmpName)
            #print(command)
            os.system('chmod +x {0}'.format(tmpName))
            os.system(command)
    
    return

def sizemag_stars(inputCat, fiveSig, filterName, field, pixScale, plotDir = 'plots/', starDir = '', sizeType = 'FWHM_IMAGE', extra= '', plotName = 'sizemag.pdf'):
    print("Running sizemag_stars(): ", '\n')
    
    ############################################
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LogNorm
    import matplotlib.backends.backend_pdf
    ##############################################

    # read in table
    tbdata = Table.read(inputCat)

    if filterName == 'u' or filterName == 'HSC-R' or  filterName == 'JH':
        classstarcut = 0.5
    else:
        classstarcut = 0.95

    
    ps = (tbdata['CLASS_STAR'] > classstarcut)
    pssize = np.array(tbdata[sizeType][ps])
    expectedSize = np.median(pssize)
    print('Expected size is ', expectedSize*pixScale)
    
    if  filterName == 'JH' or  filterName == 'HKs':
        expectedSize = 0.5/pixScale # HOLLY 22/05/25
    
    ###if filterName[0:2] == 'ch':
    #    expectedSize = 2.2 # arcsec
    #elif filterName[-1] == 'u':
    #    expectedSize = 1.0
    #else:
    #    expectedSize = 0.7
        
    minSize = 0.5*expectedSize*pixScale
    maxSize = 1.5*expectedSize*pixScale # arcseconds

    #if sizeType == 'FLUX_RADIUS':
    #    minSize = 0.4*expectedSize
    #    maxSize = 1.0*expectedSize 
    #else:
    #    minSize = 0.75*expectedSize
    #    maxSize = 1.5*expectedSize # arcseconds

    # cut it down a bit
    #print(tbdata.colnames)
    
    apMag = tbdata['MAG_APER'][:,0]
    apMag = tbdata['MAG_AUTO']
    faintMag = fiveSig - 2.0 # 1.5
    brightMag = faintMag -7.5
    
    bright = (apMag < faintMag) & (apMag > brightMag) & (tbdata[sizeType]*pixScale < maxSize) & (tbdata[sizeType]*pixScale > minSize) # cut at 20 sigma, and a reasonable size /minsize
    newtb = tbdata[bright]

    print('Cutting catalogue at five sig = {0:.1f}, {1:.2f} < m < {2:.2f}, {3:.2f} < fwhm < {4:.2f}, {5}'.format(fiveSig, faintMag, brightMag, minSize, maxSize, sizeType))
    print("There are {0}/{1} positions left after inital cuts.".format(len(newtb), len(tbdata)))

    ##############################################################
    # Now make a sanity plot
    pdf = matplotlib.backends.backend_pdf.PdfPages(plotName)
    apMagR = apMag[bright]
    size = newtb[sizeType]*pixScale
    fig = plt.figure()
    plt.hexbin(apMagR, pixScale*newtb[sizeType], bins = 'log')
    
    plt.xlabel('mAB/mag')
    plt.ylabel('{0} (/arcsec)'.format(sizeType))
    
    pdf.savefig()
    
    fig = plt.figure()
    plt.scatter(apMagR, pixScale*newtb[sizeType],s=0.2, rasterized = True)
    
    plt.xlabel('mAB/mag')
    plt.ylabel('{0} (/arcsec)'.format(sizeType))
    
    
    #############################################
    # Try and fit iteratively
    # with minimal intervention.
    from scipy.optimize import curve_fit
    
    # read in the cut from a nice file
    cuts = Table.read('star_param_{0}.txt'.format(field), format = 'ascii.commented_header')
    print(cuts)
    ff = (cuts['name'] == filterName)
    print("ff: ", ff)
    print("cuts['name']", cuts['name'])
    print("filterName: ", filterName)
    brMag = cuts['bright'][ff]
    faMag = cuts['faint'][ff]
    maxsize = cuts['upperfwhm'][ff]
    satLevel = cuts['saturate'][ff]
    nosigma = 2.0
    nosigmaup = 1.5
    faintLevel = cuts['faintest'][ff]
    
    # first get a rough estimate from the brighter stars
    #breakpoint()
    br = (apMagR < faMag) & (apMagR > brMag)  & (size < maxsize)
    firstmedian = np.median(size[br])
    mad = np.median(np.abs(size[br]-firstmedian))
    firstsigma = 1.58*mad

    # do a sigma cut here
    sc = (apMagR < faMag) & (apMagR > brMag)  & (size < firstmedian + 2.0*firstsigma) & (size > firstmedian - 2.0*firstsigma)
    
    # and remeasure
    estmedian = np.median(size[sc])
    mad = np.median(np.abs(size[sc]-firstmedian))
    estsigma = 1.58*mad
    
    #fwhmcut = estlocus + 1.58*mad
    xline = np.arange(brightMag, faintMag, 0.1)
    yline = xline*0.0 + estmedian
    plt.plot(xline, yline, color= 'k', label = 'median of bright stars')
    
    plt.plot(xline, yline+estsigma, color= 'green', label = 'one sigma')
    plt.plot(xline, yline-estsigma, color= 'green')
    plt.legend()
    
    #################################################
    # now fit these points to get a line
    keep = (apMagR < faintLevel) & (apMagR > satLevel) & (size < estmedian+nosigmaup*estsigma) & (size > estmedian-nosigma*estsigma)
    
    x = apMagR[keep]
    y = pixScale*newtb[sizeType]
    y = y[keep]
    
    A,B = curve_fit(line, x, y)[0]
    
    plt.scatter(x, y, color = 'green', s= 0.2, rasterized = True)
    
    # plot the result!
    yline = xline*A + B
    #print(yline)
    #print("Grad etc = ", A, B)
    plt.plot(xline, yline, color= 'red') #, label = 'fit to bright stars')
    
    # cut according to this fittted line
    midpoint = apMagR*A + B
    
    tofit = (apMagR < faintLevel) & (apMagR > satLevel) & (size < midpoint+nosigmaup*estsigma) & (size > midpoint-nosigma*estsigma) & (size < maxsize)

    # this is the final data
    
    locusdata = newtb[tofit]
    plt.scatter(apMagR[tofit], size[tofit], color = 'red', s= 0.2)
    
    pdf.savefig()
    
    # look at the cross section
    fig = plt.figure()
    ax = plt.subplot()
    plt.hist(size[br], 50, label = 'Bright')
    plt.hist(size[keep], 50, label = 'Keep')
    bottom, top = ax.get_ylim()
    plt.xlim([0.9*estmedian, 1.1*estmedian])
    plt.plot([estmedian, estmedian], [0, top], color = 'k')
    plt.plot([estmedian+estsigma, estmedian+estsigma], [0, top], color = 'green')
    plt.legend()
    pdf.savefig()
    
    pdf.close()
    
    # great now return the coordinates or whatever it needs
    tablename = starDir + '{0}_stars{1}.fits'.format(filterName, extra)
    locusdata.rename_column('NUMBER', 'ID')
    locusdata.rename_column('ALPHA_J2000', 'RA')
    locusdata.rename_column('DELTA_J2000', 'DEC')
    
    locusdata.write(tablename, overwrite=True) #, format = 'ascii')
    print("Table saved to ", tablename, " with {0} objects.".format(len(locusdata)))
    
    smalldata = locusdata[['ID', 'RA', 'DEC', 'X_IMAGE', 'Y_IMAGE', 'CLASS_STAR', 'FLUX_RADIUS']]
    starName = starDir + '{0}_stars{1}.ascii'.format(filterName, extra)
    smalldata.write(starName, overwrite=True, format = 'ascii')
    print("Stars written to ", starName)

    print("Inspect ", plotName, " and adjust parameters.")
    
    return starName

def line(x, A, B):
    return A*x + B

#kfilt = 'K' was probabluy meant to be kfilt = 'Ks'
def bzk_stars(inputCatalogue, bfilt = 'B', zfilt = 'Z', kfilt = 'Ks', brange = [18.0, 25.0], zrange = [19.0, 25.0], krange = [17.0, 24.5], radiusCut = [0.0001, 2.5], interceptDelta = [0.0, 0.0]):
    ''' Searching for BZK stars... to stop PSFEx breaking... '''
    print("Running bzk__stars(): ", '\n')
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LogNorm
    import matplotlib.backends.backend_pdf
    
    ## Read in the full catalogue
        
    tbdata = Table.read(inputCatalogue)
    columnNames = np.array(tbdata.colnames)
    #print tbdata.info
    #print tbdata
    
    ## First check if the filternames exist
    if np.any(columnNames == bfilt) and np.any(columnNames == zfilt) and np.any(columnNames == kfilt):
        print("I have found the filters B = ", bfilt, " Z = ", zfilt, " K = ", kfilt)
    else:
        print("Missing filters for BZK selection ",  np.any(columnNames == bfilt), np.any(columnNames == zfilt), np.any(columnNames == kfilt), bfilt, zfilt, kfilt)
        exit()
    
    gg = inputCatalogue.rfind('.fits')
    bb = inputCatalogue.rfind('/')
    basename = inputCatalogue[bb+1:gg]
    
    ## Set up the plot
    
    plotName = 'psf/' + basename + '_bzk_starsel.pdf'
    pdf = matplotlib.backends.backend_pdf.PdfPages(plotName)
               
    ## First plot the r1/2 to mag relation...
    fig = plt.figure()
    gg = (tbdata[kfilt] < krange[1]) #& (tbdata[bfilt] < brange[1]) & (tbdata[zfilt] < zrange[1]) \
        #& (tbdata[bfilt] > brange[0]) & (tbdata[zfilt] > zrange[0])
    cutdata = tbdata[gg]
    
    plt.scatter(cutdata[kfilt], cutdata['FLUX_RADIUS'],s=2, c=cutdata['CLASS_STAR'], rasterized = True)
    plt.xlim([15.0, 26.0])
    plt.ylim([0.0, 6])
    plt.xlabel(kfilt)
    
    ## plot a horizontal line corresponding to the cut
    plt.plot([10.0, 30.0], [radiusCut[0], radiusCut[0]], linewidth = 2.0)
    plt.plot([10.0, 30.0], [radiusCut[1], radiusCut[1]], linewidth = 2.0)

    pdf.savefig(fig)

## Figure Inputs
    maxx = 6.0 # for the x axis of c-o-g
    maxy = 2.0
    
## cut the catalogue a bit
## to define the stellar locus
    gg = (tbdata[kfilt] < krange[1]) & (tbdata[bfilt] < brange[1]) & (tbdata[zfilt] < zrange[1]) \
        & (tbdata[kfilt] > krange[0]) & (tbdata[bfilt] > brange[0]) & (tbdata[zfilt] > zrange[0]) \
        & (tbdata['FLUX_RADIUS'] < radiusCut[1]) & (tbdata['FLUX_RADIUS'] > radiusCut[0])
       # & (tbdata['CLASS_STAR'] > 0.90)
    
    brighttbdata = tbdata[gg]
    print("Before ", len(tbdata), " and now ", len(brighttbdata))
    
    BminZ = brighttbdata[bfilt] - brighttbdata[zfilt]
    ZminK = brighttbdata[zfilt] - brighttbdata[kfilt]
    
    fig = plt.figure()

    bins = [100, 100]
    xyrange = [[-1.0, maxx], [-1.0, 4.0]]
    #hh, locx, locy = scipy.histogram2d(BminZ, ZminK, bins = bins, range = xyrange)
    #plt.imshow(np.flipud(hh.T), cmap = 'jet', extent=np.array(xyrange).flatten(), interpolation = 'none', origin = 'upper')
    plt.scatter(BminZ, ZminK, s=1, c = brighttbdata['CLASS_STAR'])
    plt.xlim([-1.0, maxx])
    plt.ylim([-1.0, maxy])
    plt.xlabel(bfilt + ' - ' + zfilt)
    plt.ylabel(zfilt + ' - ' + kfilt)
    
    ## create cuts to get the stellar locus
    x = maxx*np.arange(100)/100.0
    #print x
    m1 = 0.75
    m2 = 0.25
    c1 = -0.95 + interceptDelta[0]
    c2 = -0.4 + interceptDelta[1]
    y1 = m1*x + c1
    y2 = m2*x + c2
    
    plt.plot(x, y1, linewidth = 2.0)
    plt.plot(x, y2, linewidth = 2.0)
    pdf.savefig(fig)
    
    ## If these cuts are adequate, cut the catalogue and make a nice
    ## BzK sample
    #BminZfull = tbdata[bfilt] - tbdata[zfilt]
    #ZminKfull = tbdata[zfilt] - tbdata[kfilt]
    BminZfull = BminZ
    ZminKfull = ZminK
    
    
    ll = (ZminKfull < BminZfull*m1 + c1) & (ZminKfull < BminZfull*m2 + c2)
    print("There are ", BminZ.size, " objects here :) and ", np.sum(ll), " stars. ")
    
    fig = plt.figure()

    ## Show these!
    locusdata = brighttbdata[ll]
    BminZlocus = locusdata[bfilt] - locusdata[zfilt]
    ZminKlocus = locusdata[zfilt] - locusdata[kfilt]
    
    plt.xlim([-1.0, maxx])
    plt.ylim([-1.0, maxy])
    plt.xlabel(bfilt + ' - ' + zfilt)
    plt.ylabel(zfilt + ' - ' + kfilt)
    plt.plot(x, y1, linewidth = 2.0)
    plt.plot(x, y2, linewidth = 2.0)
    plt.scatter(BminZlocus, ZminKlocus, s=1, c = locusdata['CLASS_STAR'])
    #plt.show()
    pdf.savefig(fig)
    
    ## Finally plot the x and y coordinates just as a check
    fig = plt.figure()
    plt.scatter(locusdata['X_IMAGE'], locusdata['Y_IMAGE'], s = 1)
    pdf.savefig(fig)
    
    #plt.savefig('bzk_test.pdf')
    pdf.close()
    plt.close()
    print("Plot saved to ", plotName)

    ## Save these stars ra and dec to a nice file!
    print(locusdata)
    #smalldata = locusdata[['ID', 'RA', 'DEC', 'X_IMAGE', 'Y_IMAGE', 'CLASS_STAR', 'FLUX_RADIUS']]
    tablename = 'psf/' + basename + '_bzk_stars.fits'

    locusdata.write(tablename, overwrite=True) #, format = 'ascii')
    print("Table saved to ", tablename)
    
    return tablename
