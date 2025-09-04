import os
import numpy as np
from datetime import datetime
from joblib import Parallel, delayed
from astropy.table import Table, Column, hstack, join, vstack, MaskedColumn

baseDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/'

def run_all_source_extractors(fieldName, detectionFilters, apDiametersAS, 
                              queue='none', memory=7, reqFilters=None, 
                              IRACapDiametersAS=None, overwrite=True, 
                              se=True, n_jobs=-1):
    """
    Run source extractor in parallel across detection filters.
    
    Parameters
    ----------
    fieldName : str
        The field name being processed.
    detectionFilters : list
        List of detection filter names.
    apDiametersAS : list
        Aperture diameters.
    queue : str
        Queue to submit jobs to.
    memory : int
        Memory requirement for job submission.
    reqFilters : list
        Required filters.
    IRACapDiametersAS : list
        IRAC aperture diameters.
    overwrite : bool
        Whether to overwrite existing output.
    se : bool
        Whether to actually run SE.
    n_jobs : int
        Number of jobs to run in parallel (-1 = all cores).
    """

    def process_filter(detectionFilt):
        if se:
            print(f"Beginning se = True at {datetime.now().strftime('%d/%m %H:%M:%S')}, filter={detectionFilt}")
            run_source_extractor(fieldName, detectionFilt, apDiametersAS,
                                 queue=queue, memory=memory, reqFilters=reqFilters,
                                 IRACapDiametersAS=IRACapDiametersAS,
                                 overwrite=overwrite)
    Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(process_filter)(detectionFilt) for detectionFilt in detectionFilters)

def mask_column(ra, dec, regFile, ax = None, verbose = False, tokeep = False, plot = False, hsc = True, xy = False, circlesFile= 'NONE'):

    """
    (Holly 03/03/25)

        Creates masks based on supplied region (reg) files. Uses the specified shape to 
        create a region for masking around positions, given by the regFile. For each shape,
        the mask is created by checking if the ra and dec values are within the shape's
        region. Pixels outside of the region are stored in Boolean mask array (bad bits).  
        """
    
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # read in the appropriate file
    ## Read in the reg file
    if os.path.isfile(regFile):
        lineNo = 0
        f = open(regFile, 'r')
        lines = f.readlines()[3:]
    else:
        lines = ''
        
    rasterize = True
    
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        ms = 0.05
        #fig = plt.figure()
        #ax = plt.subplot()
        ax.scatter(ra, dec, s = ms)
        
    # make a mask array!
    maskArray = np.ones(ra.size)
    
    if hsc:
        # apply this first
        if (ra[0] > 10) & (ra[0] < 50) & (dec[0] < 0.0) & (dec[0] > -15):
            hscfile = '/raid/scratch/data/masks/HSC_circle_xmm.reg'
        elif ra[0] > 140:
            hscfile = '/raid/scratch/data/masks/HSC_circle_cosmos.reg'
        else:
            hscfile = '/raid/scratch/data/masks/HSC_circle_cdfs.reg'
            
        if circlesFile != 'NONE':
            hscfile = circlesFile
            
        hsc = open(hscfile, 'r')
        hsclines = hsc.readlines()[3:]
        keephsc = np.ones(ra.size, dtype = bool)
        
        for lli, lline in enumerate(hsclines):
            sline = lline.strip()
            shapeType = sline[:sline.find('(')]
         #   print('line = ', sline)
            
            # split the string
            split = sline[sline.find('(') + 1:-1].split(',')
            
            # the central coordinates
            if xy:
                racentre = float(split[0])
                deccentre = float(split[1])
            else:
                c = SkyCoord(split[0], split[1], unit = (u.hourangle, u.deg))
                racentre = c.ra.degree
                deccentre = c.dec.degree
                
            radiusAS = split[2]

            if xy:
                # in pixels
                radius = float(radiusAS[:-1])
                print('Radius in pixels is ', radius)
            else:
                # weight according to the declination
                radius = float(radiusAS[:-1])/3600.0
                cosfactor = np.cos(np.pi*deccentre/180.0)
                print('Factor is ', cosfactor)
                radius = radius/cosfactor
            
            # cut!
            sep = np.sqrt((ra - racentre)**2 + (dec - deccentre)**2)

            if (ra[0] < 100.0) or xy:
                keep = sep < radius
                
                if lli < 1:
                    bb = np.copy(keep)
                else:
                    bb = keep | bb
            else:
                bb = sep > radius
                
            # plot the circle
            phis = np.arange(0,6.28,0.01)
            xcircle = racentre + radius*np.cos(phis)
            ycircle = deccentre + radius*np.sin(phis) 
            if plot:
                ax.plot(xcircle, ycircle, c='red',ls=':')
            
       # print(bb[0:5])
        if (ra[0] < 100.0) or xy:
            bb = np.invert(bb)
            #    print(bb[0:5])
            #print('Hello, inverting, removing {0}/{1} points'.format(np.sum(bb), ra.size))
            
        maskArray[bb] = 0.0
        if plot:
            ax.scatter(ra[bb], dec[bb], color = 'green', s = ms)
                    
    
    ##############################################
    for li, line in enumerate(lines):
        
        sline = line.strip()
        shapeType = sline[:sline.find('(')]
        
        ## split the string
        split = sline[sline.find('(') + 1:-1].split(',')
        
        if len(split) < 2:
            break

        if xy:
            racentre = float(split[0])
            deccentre = float(split[1])
        else:
            # the central coordinates
            c = SkyCoord(split[0], split[1], unit = (u.hourangle, u.deg))
            racentre = c.ra.degree
            deccentre = c.dec.degree
        if verbose:
            print(c, ra, dec)
        
        if shapeType == 'box':

            # mask in the box!
            raWidthArcsec = split[2]
            decWidthArcsec = split[3]
            if xy:
                raWidth = float(raWidthArcsec[:-1])/2.0
                decWidth = float(decWidthArcsec[:-1])/2.0
            else:
                raWidth = float(raWidthArcsec[:-1])/7200.0
                decWidth = float(decWidthArcsec[:-1])/7200.0 # convert to degrees

            # for the cosine
            #cosfactor = np.cos(np.pi*dec/180.0)
            #raWidth = raWidth/cosfactor
            
            # get the ra and dec dimentions
            decmin = deccentre - decWidth
            decmax = deccentre + decWidth
            ramin = racentre - raWidth
            ramax = racentre + raWidth

            #print(ramin, ramax, decmin, decmax, raWidth, decWidth)
            bb = (ra > ramin) & (ra < ramax) & (dec < decmax) & (dec > decmin)

            # plot the box!
            if plot:
                xbox = [ramin, ramin, ramax, ramax, ramin]
                ybox = [decmin, decmax, decmax, decmin, decmin]
                plt.plot(xbox, ybox, color = 'red')
            
        elif shapeType == 'ellipse':
            print("Ellipse.")
            
        elif shapeType == 'circle':
            if verbose:
                print("Circle.", split)
            
            radiusAS = split[2]
            if xy:
                radius = float(radiusAS[:-1])
            else:
                radius = float(radiusAS[:-1])/3600.0 # in degrees

            #cosfactor = np.cos(np.pi*deccentre/180.0)
            #radius = radius/cosfactor
            
            #print("The radius is ", radius, line)
            #exit()
            
            # cut!
            sep = np.sqrt((ra - racentre)**2 + (dec - deccentre)**2)

            bb = (sep < radius)
            
            # plot the circle
            if plot:
                phis=np.arange(0,6.28,0.01)
                xcircle = racentre + radius*np.cos(phis)
                ycircle = deccentre + radius*np.sin(phis)
                ax.plot(xcircle, ycircle, c='red',ls='-')
            
        elif shapeType == 'polygon':
            
            
            # For video (masked by Boris), and me...
            # First get all the corners into an array
            racorn = np.array(split[0::2])
            deccorn = np.array(split[1::2])

            if xy:
                x = racorn 
                y = deccorn
            else:
                
                c = SkyCoord(racorn, deccorn, unit = (u.hourangle, u.deg))
                x = c.ra.degree
                y = c.dec.degree

            #            print("Polygon", y)
            
            verbose = False
#            if (np.max(y) > -4.3) & (np.min(y) < -4.3):
#                verbose = True
#                print('Check, original ', x,y)

            
            # Make sure the order is right
            if np.argmax(y) > 0:
                value = np.argmax(y)
                x = np.roll(x, -value)
                y = np.roll(y, -value)

                if verbose:
                    print('new x, y = ', x, y)
                                
            if x[1] < x[-1]:
                # change direction to clockwise
                x = np.flip(x, axis = 0)
                y = np.flip(y, axis = 0)
                x = np.roll(x, 1)
                y = np.roll(y,1)
                if verbose:
                    print('new new x, y = ', x, y)
            
            # Section 1
            # define the y range for this section
            ymax = y[0]
            if y[1] > y[-1]:
                ymin = y[1]
            else:
                ymin = y[-1]
            
            # now get gradients of the lines
            mAB = (y[1]-y[0])/(x[1]-x[0])
            mAD = (y[0]-y[-1])/(x[0]-x[-1])
            mBC = (y[2]-y[1])/(x[2]-x[1])
            mCD = (y[3]-y[2])/(x[3]-x[2])
            
            cAD = y[0]-mAD*x[0]
            cAB = y[1]-mAB*x[1]
            cBC = y[2]-mBC*x[2]
            cCD = y[3]-mCD*x[3]

            lineAD = (dec-cAD)/mAD
            lineAB = (dec-cAB)/mAB
            lineCD = (dec-cCD)/mCD
            lineBC =  (dec-cBC)/mBC
            
            # first test
            condy = (dec > ymin) & (dec <= ymax)
            cond = condy & (ra < lineAB) & (ra > lineAD)
            
            # second section
            if y[1] > y[-1]:
                ymax = y[1]

                if y[2] > y[-1]:
                    ymin = y[2]
                else:
                    ymin = y[-1]
                
                # need BC gradient
                lineleft = lineAD
                lineright =lineBC
            else:
                ymax = y[-1]
                ymin = y[1]
                lineleft = lineCD
                lineright = lineAB

            condy = (dec > ymin) & (dec <= ymax)
            condtwo = condy & (ra < lineright) & (ra > lineleft)
            
            # final section
            # final section, three options?
            # C bottom:
            if (y[2] < y[-1]) & (y[2] < y[1]):
                ymin = y[2]
                lineleft = lineCD
                lineright = lineBC
                if y[1] > y[-1]:
                    ymax = y[-1]
                else:
                    ymax = y[1]
                    
            elif (y[1] < y[2]) & (y[1] < y[-1]):
                # B bottom
                ymin = y[1]
                lineleft = lineBC
                lineright = lineAB
                ymax = y[2]
            elif (y[-1] < y[2]) & (y[-1] < y[1]):
                # D bottom
                ymin = y[-1]
                lineleft = lineAD
                lineright = lineCD
                ymax = y[2]
                
            
            condy = (dec > ymin) & (dec <= ymax)
            condthree = condy & (ra < lineright) & (ra > lineleft)
            
            bb = cond | condtwo | condthree
            #if np.sum(bb) < 1:
            #    print('Nothing inside...')
                
            if verbose:
                print('Exciting.')
                exit()
            
            # this gives inside of the polygon
            if plot:
                xfull = np.append(x, x[0])
                yfull = np.append(y, y[0])
                ax.plot(xfull, yfull, color = 'red', zorder = 10)
            
            
        else:
            print("Shape not recognised: ", shapeType)
    
        # mask for this shape
        if plot:
            ax.scatter(ra[bb], dec[bb], color = 'green', s= ms)    
        maskArray[bb] = 0

    print('Created mask with {0} bad bits.'.format(maskArray.size -np.sum(maskArray)))
    if plot:
        ax.set_rasterized(rasterize)
        ax.set_ylim([np.min(dec), np.max(dec)])
        ax.set_xlim([np.min(ra), np.max(ra)])
        #fig.savefig(plotname)
        #print('check plot, ', plotname)
        #exit()
    # invert to mask in next code
    if tokeep:
        return maskArray.astype(bool)
    else:
        return np.invert(maskArray.astype(bool))


def set_negativeflux_zero(catName):
    print("running set_negativeflux_zero: " ,'\n')
    newName = catName[:-5] + '_noneg.fits'
    print('new name = ', newName)

    tb = Table.read(catName)
    
    colnames = np.array(tb.colnames)

    for c, col in enumerate(colnames):

        print(col)
        if col[0:4] == 'flux':

            print('fixing column ', col)
            
            neg = (tb[col] < 0.0)
            if np.any(neg):
                print('setting {0} bad'.format(np.sum(neg)))
                tb[col][neg] = 0.0

    # also cut
    notgarbage = ((tb['HSC-I'] > -90.0) & (tb['HSC-I'] < 80.0)) | ((tb['CFHT-iy'] < 80.0) & (tb['CFHT-iy'] > -90.0))
    print('There are {0}/{1} remaining.'.format(np.sum(notgarbage), len(tb)))
                                                                
                                                                  
    tb[notgarbage].write(newName, overwrite =True)
    print('Fluxes below zero fixed in file ', newName)
    print("set_negative_flux line 35 finished", '\n')
    return

    

def combine_xmmcdfs(detFilter, basename, reqCats, fieldName, outDir = '', requiredSpit = ['2.8', '2.8', '2.8'], individual = True, cuts = [-99.0, -99.0, -99.0], remove_counts = True):
    """(Holly 03/03/25)
        Creates a fits file with catalogues for available filters combined into one.

    """
    print("running combine_xmmcdfs: " ,'\n')
    # if individual, make one catalogue for each reqCat
    # otherwise combine into one
    hsc = True
    if fieldName == 'XMM':
        fields = ['XMM1', 'XMM2', 'XMM3']
        finalstring = 'XMMFULL'
    elif fieldName == 'CDFS':
        fields = ['CDFS1', 'CDFS2', 'CDFS3']
        finalstring = 'CDFSFULL'
        hsc = False
    else:
        print('Error, field {0} not recognised'.format(fieldName))
        exit()
        
    if detFilter == 'CFHT-iy':
        fields = ['XMM3']
    
    from datetime import date
    today = date.today()
    today = today.strftime("%Y_%m_%d")

    if individual == False:
        
        # make a final filename
        finalFilename = outDir + finalstring + basename + '{1}_{0}.fits'.format(today, detFilter) 

    for r, reqCat in enumerate(reqCats):
        
        if requiredSpit[r] == 'NONE':
            spitstring = ''
            optstring = reqCat
        else:
            spitstring = '_IRAC{0}as'.format(requiredSpit[r])
            optstring = '{0}as'.format(reqCat)
        
        if individual:
            finalFilename = outDir + finalstring + basename + '{1}_{2}{3}_{0}.fits'.format(today, detFilter, optstring, spitstring)
            finalCutFilename = outDir + finalstring + basename + '{1}_cut5sig_{2}{3}_{0}.fits'.format(today, detFilter, optstring, spitstring)
            
        # for each field first combine the different catalogues, with sensible names!
        for fi, field in enumerate(fields):
            
            fieldDir = '../catalogues/{0}/det_{1}/'.format(field, detFilter)
            
            # read in the base catalogue
            filename = '{0}{1}{2}{5}_{3}{4}.fits'.format(fieldDir, field, basename, optstring, spitstring, detFilter)
            mags = Table.read(filename)
            print('Including cat ', filename)
            
            # get the flux cat
            fluxname =  '{0}{1}{2}{5}_{3}{4}_cgs.fits'.format(fieldDir, field, basename, optstring, spitstring, detFilter)
            if os.path.isfile(fluxname):
                flux = Table.read(fluxname)
                
                
                #check for infinities
                check = np.isinf(flux['flux_' + detFilter])
                if np.any(check):
                    print('There are infinities ', np.sum(check), ' for filter ', detFilter)
                    flux[check].write('testcat.fits', overwrite = True)
                    print('Check the testcat.fits file')
                    exit()
                else:
                    print('No infinities ')
                    
            else:
                print('Error no flux file ', fluxname)
            
            # rename columns?
            # join these together
            #       fulltb = join(base, baseflux, keys = ['ID'])#table_names = ['MAG', 'RAD']
            #fulltable = base

            yess = False
            yesss = False
            yessss = False

            if ('HSC_circle' in mags.colnames):
                
                mags.remove_column('HSC_circle')
            #mags.remove_column('VISTA_tile')
            #mags.remove_column('VISTA_tile_used')
            
            # add label columns
            colnamearray = np.array(mags.colnames)
            exe = (colnamearray == 'VISTA_tile')
            if np.any(exe) == False:
                
                vistatile = which_field(mags['RA'], mags['DEC'], overlap =True)
                newcol = Column(vistatile, name = 'VISTA_tile')
                mags.add_column(newcol)
                yess = True
                
            exe = (colnamearray == 'VISTA_tile_used')
            if np.any(exe) == False:
                # also add final column
                finaltile = which_field(mags['RA'], mags['DEC'])
                newcol = Column(finaltile, name = 'VISTA_tile_used')
                mags.add_column(newcol)
                yesss = True


            if hsc:
                exe = (colnamearray == 'HSC_circle')
                if np.any(exe) == False:
                    # and add a column of the HSC field
                    finaltile, hsctile = which_field(mags['RA'], mags['DEC'], hsc = True)
                    newcol = Column(hsctile.astype(str), name = 'HSC_circle')
                    mags.add_column(newcol)
                    yessss = True
            else:
                yessss = True
                
            if yess | yesss | yessss:
                # now save this
                mags.write(filename, overwrite = True)
                
            # join mags and flux together
            
            # rename the mags columns
            for ii, colhere in enumerate(np.array(mags.colnames)):
                if colhere[0:4] == 'flux':
                    #                print(colhere, 'counts ', colhere[5:])

                    if remove_counts:
                        mags.remove_column(colhere)
                    else:
                        mags.rename_column(colhere, 'counts_' + colhere[5:])
                    
            
            #       print(mags.colnames)
            if os.path.isfile(fluxname):
                fulltb = join(mags, flux, keys = ['ID', 'X_IMAGE', 'Y_IMAGE', 'RA', 'DEC'])
            else:
                fulltb = mags
                
            # remove sources that do not have the right VISTA flag
            keep = (fulltb['VISTA_tile_used'] == field)
            print('Keeping {0}/{1} that are in {2}'.format(np.sum(keep), len(fulltb), field))
            
            fulltb = fulltb[keep]
            
            if fi < 1:
                finaltb = fulltb
                
            else:
                
                finaltb = vstack([finaltb, fulltb])
                
            
            # make a plot of the different regions

        finaltb.sort('RA')
        
        # change the id number, make new column
        finaltb.rename_column('ID', 'ID_se')
        # move this column to the end?
        
        # make a new column
        newid = Column(np.arange(len(finaltb), dtype = int), name = 'ID')
        finaltb.add_column(newid, index = 0)
        
        # fix masked values
        for ii, colname in enumerate(np.array(finaltb.colnames)):
            finaltb[colname].fill_value = -99.0
            
        finaltb = finaltb.filled()
        
        finaltb.write(finalFilename, overwrite = True)
        
        print('The final file has {0} objects'.format(len(finaltb)))
        print('The final filename = ', finalFilename)

        if cuts[r] > 0.0:
            # also cut at five sigma
            keep = (finaltb[detFilter] <= cuts[r])
            print('Cutting at five sigma of {0} in {1}, leaves {2}/{3}'.format(cuts[r], finalCutFilename, np.sum(keep), len(finaltb)))
            finaltb[keep].write(finalCutFilename, overwrite = True)
            print('Saved to ', finalCutFilename)
            exit()
            
    return


def which_field(ra, dec, hsc= False, verbose = True, overlap = False):
    print("Running which_field: ", '\n')
    ## load in the different edges
#    xmm3ra = 35.995833 ### 36.125
#    xmm3dec = -5.1831944
#    xmm23ralow = 36.125 ## bit near the cutout crap edge in XMM3
    
    #xmm12ra_high = 35.048750
    #xmm12ra_low = 34.922500
#    xmm12rahigh = 34.927917
#    xmm12ra = 35.038029
#    decide_dec = -4.4309639 ##-4.84 ## to decide which ra cut here

    ###### New definition for DR2 catalogues
    # for just labelling the different overlaps
    lab_xmm2min = 34.925
    lab_xmm2minNotch = -5.285

    lab_xmm3min = 35.97
    lab_xmm3minNotch = -5.195
    
    # and for inclusion in final catalogue
    xmm2min = 34.935
    xmm2minNotch = -5.285 #-5.275
    xmm2minNotch = -5.275
    xmm3min = 35.98
    xmm3minNotch = -5.195 #-5.19 to make the line straight!
    xmm3minNotch = -5.19
    xmm1max = 35.048
    xmm1maxNotch =-4.436
    xmm2max = 36.13

    # and for cdfs
    # for inclusion in final catalogue 
    cdfs_deccut = -28.16
    cdfs_racut = 53.31

    # for labelling the different regions
    cdfs12_up = -28.03
    cdfs12_down = -28.17
    cdfs3_left = 53.255
    cdfs3_right = 53.363

    field = np.chararray(ra.shape, itemsize = 20)
    field[:] = 'UNKNOWN'
    
    # Check XMM
    xmm = (dec < -3.0) & (dec > -7.0)
    field[xmm] = 'XMM'

    cdfs = (ra > 40) & (ra < 60) & (dec < -20.0) & (dec > -30)
    field[cdfs] = 'CDFS'
    
    #########################
    # first just label the fields, including overlap
    if np.any(xmm):
        if overlap:
            one = (ra < xmm1max) 
            if np.any(one):
                field[one] = 'XMM1'
                
            ot = (ra < xmm1max) & (ra > lab_xmm2min)
            if np.any(ot):
                field[ot] = 'XMM12'
        
            two = (ra > xmm1max) & (ra < lab_xmm3min)
            if np.any(two):
                field[two] = 'XMM2'
                
            tt = (ra > lab_xmm3min) & (ra < xmm2max)
            if np.any(tt):
                field[tt] =	'XMM23'
                
            three = (ra > xmm2max)
            if np.any(three):
                field[three] = 'XMM3'
                
        else:
            # properly label
            # one
            
            field[:] = 'XMM2'
            
            one = ((dec <xmm2minNotch) & (ra < xmm1max)) | ((dec >= xmm2minNotch) & (ra < xmm2min))
            if np.any(one):
                field[one] = 'XMM1'
            
            three = ((dec <xmm3minNotch) & (ra > xmm2max)) | ((dec >= xmm3minNotch) & (ra > xmm3min))
            if np.any(three):
                field[three] = 'XMM3'
        
    if np.any(cdfs):
        print('In CDFS')

        if overlap:

            one = (ra < cdfs3_right) & (dec > cdfs12_down)
            if np.any(one):
                field[one] = 'CDFS1'

            ot = (ra < cdfs3_right) & (dec > cdfs12_down) & (dec < cdfs12_up)
            if np.any(ot):
                field[ot] = 'CDFS12'

            two = (ra < cdfs3_right) & (dec < cdfs12_down)
            if np.any(two):
                field[two] = 'CDFS2'

            three = (ra >= cdfs3_right)
            print('The sources in three are ', np.sum(three))
            if np.any(three):
                field[three] = 'CDFS3'

            oth = (ra < cdfs3_right) & (ra > cdfs3_left) & (dec > cdfs12_up)
            if np.any(oth):
                field[oth] = 'CDFS13'
            
            twothree = (ra < cdfs3_right) & (ra > cdfs3_left) & (dec < cdfs12_down)
            if np.any(twothree):
                field[twothree] = 'CDFS23'
            
            ott = (ra < cdfs3_right) & (ra > cdfs3_left) & (dec > cdfs12_down) & (dec < cdfs12_up)
            if np.any(ott):
                field[ott] = 'CDFS123'


        else:
            print('HERE looking at overlap ')
            # now for splitting the catalogues
            field[:] = 'CDFS2'

            one = (ra < cdfs_racut) & (dec > cdfs_deccut)
            if np.any(one):
                field[one] = 'CDFS1'

            three = (ra >= cdfs_racut)
            print('There are {0} in three'.format(np.sum(three)))
            if np.any(three):
                field[three] = 'CDFS3'
            
    ## In the UDS (which current code)
#    three = (((ra >= xmm3ra) & (dec >= xmm3dec)) | ((ra >= xmm23ralow) & (dec < xmm3dec)))
#    if np.any(three):
#        field[three] = 'XMM3_DEEP'     
    
#    one = (((dec > decide_dec) & (ra < xmm12rahigh)) | ((dec <= decide_dec) & (ra < xmm12ra)))
#    if np.sum(one) > 0:
#        field[one] = 'XMM1_UDEEP'
        
    #print field
    #exit()   
    ## Finally, work out where we are in the UDEEP or DEEP part of HSC
    ## for the UDEEP HSC part
    ##circle(2:18:33.341,-4:51:28.132,2894.368")
    #circlera = 34.638917
    #circledec = -4.8578144
    #radius = 2894.368/3600.0 ## degrees

    if hsc:
        # also make another field column, with the circle
        #fieldHSC = np.chararray(ra.shape, itemsize = 10)
        fieldHSC = np.empty(ra.shape, dtype = 'object')
        #hscfile = '/users/bowlerr/vardygroupshare/data/codes/area/hsc_circles_deg_official.txt'
        hscfile = '/mnt/vardy/vardygroupshare/data/masks/hsc_circles_deg_mymask.txt'
        #hscFields = np.genfromtxt(hscfile, names = True, dtype = None, encoding = None)
        hscFields = Table.read(hscfile, format = 'ascii.commented_header')
        hscFields = hscFields[:-1]
        fieldHSC[:] = ''
        
    # use the official coordinates now
 #   circlera = 34.565000000000005
 #   circledec = -4.85
#    radius = 2700.0/3600.0
#    print "Using new circle now."
    
#    r = np.sqrt((ra - circlera)**2 + (dec - circledec)**2)
#    udeep = (field == 'XMM2_DEEP') & (r < radius)
#    if np.sum(udeep) >0:
#        field[udeep] = 'XMM2_UDEEP'

        # loop through labelling
        for l, lab in enumerate(hscFields['label']):
            circlera = hscFields['RA'][l]
            circledec = hscFields['DEC'][l]
            radius = hscFields['radius'][l]/3600.0

            # include the cos factor
            cosfactor = np.cos(np.pi*circledec/180.0)
            print('Applying the cos factor of ', cosfactor)
            radius = radius/cosfactor
            
            r = np.sqrt((ra - circlera)**2 + (dec - circledec)**2)
            here = (r < radius)
            if np.sum(here) > 0:
                
                fieldHSC[here] = fieldHSC[here] + lab#.replace('\x00','')#.decode('UTF-8')

    # Now for COSMOS
    cos= (ra > 140.0) & (ra < 160)
    field[cos] = 'COSMOS'
    
    #else:
    #    print "Field not recognised.  Code only works for XMM or COSMOS, sorry."
    #    return "none"    

    
    if hsc:
        return field, fieldHSC
    else:
        return field.astype(str)
    
def cut_out(imageName, ra, dec, idarray, filterName = 'blah', stampSizeAS = 10.0, header = False, overwrite = False, outDir = '', xyin = False, verbose = True, individual = False, astropycut = True):
    # ammended to match the once in image_codes.py
    # which uses astropy and fixes the astrometry!
    print("Running cut_out: ", '\n')
    from astropy.io import fits
    from astropy.wcs import WCS
    
    if filterName != 'blah':
        fitsOutput = outDir + filterName + '_tmp.fits'
    else:
        fitsOutput = 'stamps_output.fits'
 
    ## First check if the output exists!
    if overwrite == False:
        if os.path.isfile(fitsOutput):
            print("Fits cut-out file ", fitsOutput, " already exists!")
            print("cut_out finished", '\n')
            return fitsOutput


    ## Code to cut out stamps from a larger array
    ## open image
    hdu_image = fits.open(imageName, memmap=True)
    hdu = hdu_image[0]
    imageHeader = hdu.header
    imageData = hdu.data
    
    ## convert ra/dec to x,y pixels
    if xyin:
        x = ra
        y = dec
    else:
        w = WCS(imageHeader)
        x, y = w.all_world2pix(ra, dec, 1)
    
    ## get width, cut out
    pixScale = abs(imageHeader['CD1_1'])*3600.0
    widthPix = int(np.round(stampSizeAS/pixScale))
    if widthPix % 2 == 0:
        halfWidth = widthPix/2
    else:
        halfWidth = (widthPix-1)/2

    halfWidth = int(halfWidth)
    
    if verbose:
        print("The half width is ", halfWidth)
        print("Width = ", widthPix, pixScale, stampSizeAS)
    
    ###exit()
    ## create mef for all the objects here
    if individual == False:
        new_hdul = fits.HDUList()
        new_hdul.append(fits.PrimaryHDU())
    
    # loop through each object
    for oi, obj in enumerate(idarray):

        fitsout = outDir + 'ID{0}.fits'.format(obj)
        if overwrite == False:
            if os.path.isfile(fitsout):
                continue

        if astropycut:
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            
            from astropy.nddata import Cutout2D
            position = SkyCoord(ra[oi]*u.deg, dec[oi]*u.deg)
            sizeAS = u.Quantity((stampSizeAS, stampSizeAS), u.arcsec)
            stamp = Cutout2D(imageData, position, sizeAS, wcs=w)
            theheader = stamp.wcs.to_header()
            data = stamp.data
            
        else:
            
            xcentre = int(np.round(x[oi]))
            ycentre = int(np.round(y[oi]))
            
            xmin = xcentre - halfWidth
            xmax = xcentre + halfWidth+1
            ymin = ycentre - halfWidth
            ymax = ycentre + halfWidth+1
            
            #        print(xmin, xmax, ymin, ymax, xcentre, halfWidth)
            
            data = imageData[ymin:ymax, xmin:xmax]     

        ii = (data != 0)
        
        if individual:
            if np.any(ii):
                new_hdul = fits.HDUList()
                new_hdul.append(fits.PrimaryHDU())
                new_hdul.append(fits.ImageHDU(data))
                new_hdul.writeto(fitsout, overwrite = True)
                
        else:
            if astropycut:
                # add the hdr
                new_hdul.append(fits.ImageHDU(data, name = '{0}'.format(obj), header = theheader))    
            else:
                new_hdul.append(fits.ImageHDU(data, name = '{0}'.format(obj)))

    if individual == False:
        new_hdul.writeto(fitsOutput, overwrite = True)

    hdu_image.close()
   
    if verbose:
        print(fits.info(fitsOutput))

    print('Saved to ', fitsOutput)
    print("cut_out finished", '\n')
    return fitsOutput

def apply_mask(inputCatalogue, field, verbose = True, cutband = 'none', save = False, outputTable = 'None', hsc = True, imageDir = baseDir + 'data/', maskValueFlux = -99.0, maskValueMag = -99.0):
    """(Holly 03/03/25)
        Writes mask files by calling mask_col to get pixels to mask/keep. Sets masked magnitude pixel values to -99 (white) and if flux column exists, sets corresponding pixels to -99. Plots masked filter image and saves to pdf. 
    
    RB - cutband was used to make sure objects that were selected first in an optical catalogue also had the K-band coverage, so an extra masking step. For our high-z catalogues it's not needed. 
    
    """    

    print("Running apply_mask: ", '\n')

    from astropy.table import Table #, Column, hstack
    import astropy.io.fits as fits
    from astropy import units as u
    import numpy as np
    import os
    print(f"Looking for: {inputCatalogue}")
    print("Exists?", os.path.exists(inputCatalogue))
    print("Absolute path:", os.path.abspath(inputCatalogue))

    tb = Table.read(inputCatalogue)
    #tb = Table.read(inputCatalogue)
    raData = np.array(tb['RA'].quantity)
    decData = np.array(tb['DEC'].quantity)
    maskValue = np.zeros(raData.size)
    # If save, then I have provided the catalogue NAME not table
    #if save:
    # Read in

    # Loop through the filters in the catalogue
    dirHere = imageDir + field + '/'
    # read in filters
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()

    imagedata = Table.read(inputFile, format = 'ascii.commented_header')
#    imagedata = imagedata[2:3]
    
    availableFilters = np.array(imagedata['Name'])
    colnames = np.array(tb.colnames)
    #print(colnames)
    
    # where the reg files are
    regDir = '{0}masks/{1}/'.format(imageDir,field)

#     if decData[0] < 0:
#        field = 'XMM'
#    else:
#        field = 'COSMOS'
    #newFile = inputCatalogue  # Assign a default value in case it's never set
    
    if (save == True) and (cutband == 'none'):
        # define a new output file
        print("Need to define a new output file.")
        kk = inputCatalogue.find('UNMASKED')
        if kk > 0:
            newFile = inputCatalogue[:kk] + inputCatalogue[kk+2:]
            print("New file will be ", newFile)
        else:
            print("Please input file name, it's not obvious how to rename.")
            print(inputCatalogue)
            newFile = inputCatalogue
            print("new file assigned: ", newFile)
            
    elif (save == True) and (cutband != 'none'):
        #if cutband == 'Ks':
        if cutband == 'K':
            kk = inputCatalogue.find('MASKED')
            if kk > 0:
                newFile = inputCatalogue[:kk] + inputCatalogue[kk:kk+4] + 'VISTA' + inputCatalogue[kk+6:]
            else:
                print('Error, have not defined output file')
        else:
            kk = inputCatalogue.find('MASKED')
            jj = inputCatalogue.find('MASKVISTA')
            if kk > 0:
                newFile = inputCatalogue[:kk] + inputCatalogue[kk:kk+4] + 'DET' + inputCatalogue[kk+4:]
            elif jj > 0:
                newFile = inputCatalogue[:jj] + inputCatalogue[jj:jj+9] + 'DET' + inputCatalogue[jj+9:]
            else:
                print('Error, have not defined output file')
                exit()
                
    else:
        
        print("Please input output file name, it's not obvious how to rename this.")
        print(inputCatalogue)
        # exit code
        exit()

    if cutband != 'none':

        # cut according to a given filter!
        fi = (availableFilters == cutband)
        maskFile = imagedata['mask'][fi]
        #print(maskFile)

        # check xy!
        if len(maskFile) == 0:
            print(f"No mask file found for cutband {cutband}")
        else:
            kk = maskFile[0].find('xy')
            if kk != -1:
                print("Mask filename contains 'xy'")

        if kk > 0:
            tokeep = mask_column(tb['X_IMAGE'], tb['Y_IMAGE'], regDir + maskFile[0], tokeep = True, xy= True)
        else:
            tokeep = mask_column(tb['RA'], tb['DEC'], regDir + maskFile[0], tokeep = True)

        newtb = tb[tokeep]

        print("Only including data that has {0}-band coverage.".format(cutband))
        
        # save new file
        if save:
            newtb.write(newFile, overwrite = True)
            print("New masked file {0}/{1} objects ({2:.1f} percent) at : {3} ".format(len(newtb), len(tb), 100.0*float(len(newtb))/float(len(tb)), newFile))           
            return
        
        else:
            return newtb
        
    # make a nice diagnostic plot
    plotname = newFile.rstrip('.fits') + '.pdf'
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages(plotname)
    
    # loop through the filters
    for fi, filt in enumerate(availableFilters):

        # extract the appropriate column
        jj = (colnames == filt)
        
        if np.any(jj):
            
            # read in the mask file
            maskFile = imagedata['mask'][fi]

            # if we are doing hsc, then do that here
            if filt[0:3] == 'HSC':
  #              if ra[0] > 100:
 #                   hscfile = '../masks/HSC_circle_cosmos.reg'
   #             else:
    #                hscfile = '../masks/HSC_circle_xmm.reg'

                hscflag = True
            else:
#                hscfile = 'none'
                hscflag = False
                
            if (maskFile != 'none') or hscflag:

                fig = plt.figure()
                ax = plt.subplot()
                print(inputCatalogue)
                print("Masking filter {0} with file {1}".format(filt, maskFile))

                # check xy
                kk = maskFile.find('xy')
                if kk > 0:                   
                    # mask the mags
                    print('Masking with X/Y', maskFile)
                    if filt == 'HSC-R':
                        circlesFile = regDir + 'HSC_circle_cdfs_R_xy.reg'
                    else:
                        circlesFile = regDir + 'HSC_circle_cdfs_xy.reg'
                    
                    tomask = mask_column(tb['X_IMAGE'], tb['Y_IMAGE'], regDir + maskFile, plot = True, ax = ax, circlesFile = circlesFile, hsc = hscflag, xy= True)
                else:
                    tomask = mask_column(tb['RA'], tb['DEC'], regDir + maskFile, plot = True, ax = ax, hsc = hscflag)

                tb[filt][tomask] = maskValueMag

                # see if there is a flux column
                columns = np.array(tb.colnames)
                ff = (columns == 'flux_' + filt)
                if np.any(ff):
                    tb['flux_'+ filt][tomask] = maskValueFlux
                else:
                    print('No flux column for ', filt)

                ax.set_title('Masked for filter {0} in {1}'.format(filt, field))
                pdf.savefig(fig)
                plt.close()
                
            else:
                print("No mask available for filter {0}".format(filt))
                
        else:
            print("No filter {0} in file.".format(filt))

        
    pdf.close()
    
    # save new file
    if save:
        tb.write(newFile, overwrite = True)
        print("New masked file saved to ", newFile)
        print('Check fig saved to ', plotname)
        
        return

    else:
        return tb
    
    
    
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # read in the appropriate file
    ## Read in the reg file
    if os.path.isfile(regFile):
        lineNo = 0
        f = open(regFile, 'r')
        lines = f.readlines()[3:]
    else:
        lines = ''
        
    rasterize = True
    
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        ms = 0.05
        #fig = plt.figure()
        #ax = plt.subplot()
        ax.scatter(ra, dec, s = ms)
        
    # make a mask array!
    maskArray = np.ones(ra.size)
    
    if hsc:
        # apply this first
        if (ra[0] > 10) & (ra[0] < 50) & (dec[0] < 0.0) & (dec[0] > -15):
            hscfile = '../masks/HSC_circle_xmm.reg'
        elif ra[0] > 140:
            hscfile = '../masks/HSC_circle_cosmos.reg'
        else:
            hscfile = '../masks/HSC_circle_cdfs.reg'
            
        if circlesFile != 'NONE':
            hscfile = circlesFile
            
        hsc = open(hscfile, 'r')
        hsclines = hsc.readlines()[3:]
        keephsc = np.ones(ra.size, dtype = bool)
        
        for lli, lline in enumerate(hsclines):
            sline = lline.strip()
            shapeType = sline[:sline.find('(')]
         #   print('line = ', sline)
            
            # split the string
            split = sline[sline.find('(') + 1:-1].split(',')
            
            # the central coordinates
            if xy:
                racentre = float(split[0])
                deccentre = float(split[1])
            else:
                c = SkyCoord(split[0], split[1], unit = (u.hourangle, u.deg))
                racentre = c.ra.degree
                deccentre = c.dec.degree
                
            radiusAS = split[2]

            if xy:
                # in pixels
                radius = float(radiusAS[:-1])
                print('Radius in pixels is ', radius)
            else:
                # weight according to the declination
                radius = float(radiusAS[:-1])/3600.0
                cosfactor = np.cos(np.pi*deccentre/180.0) # accounts for DEC distortion when converting angular distances to physical dists. dec = 0 -> cosfactor = 1 means no adjustment required 
                print('Factor is ', cosfactor)
                radius = radius/cosfactor
            
            # cut!
            sep = np.sqrt((ra - racentre)**2 + (dec - deccentre)**2)

            if (ra[0] < 100.0) or xy:
                keep = sep < radius
                
                if lli < 1:
                    bb = np.copy(keep)
                else:
                    bb = keep | bb
            else:
                bb = sep > radius
                
            # plot the circle
            phis = np.arange(0,6.28,0.01)
            xcircle = racentre + radius*np.cos(phis)
            ycircle = deccentre + radius*np.sin(phis) 
            if plot:
                ax.plot(xcircle, ycircle, c='red',ls=':')
            
       # print(bb[0:5])
        if (ra[0] < 100.0) or xy:
            bb = np.invert(bb)
            #    print(bb[0:5])
            #print('Hello, inverting, removing {0}/{1} points'.format(np.sum(bb), ra.size))
            
        maskArray[bb] = 0.0
        if plot:
            ax.scatter(ra[bb], dec[bb], color = 'green', s = ms)
                    
    
    ##############################################
    for li, line in enumerate(lines):
        
        sline = line.strip()
        shapeType = sline[:sline.find('(')]
        
        ## split the string
        split = sline[sline.find('(') + 1:-1].split(',')
        
        if len(split) < 2:
            break

        if xy:
            racentre = float(split[0])
            deccentre = float(split[1])
        else:
            # the central coordinates
            c = SkyCoord(split[0], split[1], unit = (u.hourangle, u.deg))
            racentre = c.ra.degree
            deccentre = c.dec.degree
        if verbose:
            print(c, ra, dec)
        
        if shapeType == 'box':

            # mask in the box!
            raWidthArcsec = split[2]
            decWidthArcsec = split[3]
            if xy:
                raWidth = float(raWidthArcsec[:-1])/2.0
                decWidth = float(decWidthArcsec[:-1])/2.0
            else:
                raWidth = float(raWidthArcsec[:-1])/7200.0
                decWidth = float(decWidthArcsec[:-1])/7200.0 # convert to degrees

            # for the cosine
            #cosfactor = np.cos(np.pi*dec/180.0)
            #raWidth = raWidth/cosfactor
            
            # get the ra and dec dimentions
            decmin = deccentre - decWidth

            decmax = deccentre + decWidth
            ramin = racentre - raWidth
            ramax = racentre + raWidth

            #print(ramin, ramax, decmin, decmax, raWidth, decWidth)
            bb = (ra > ramin) & (ra < ramax) & (dec < decmax) & (dec > decmin)

            # plot the box!
            if plot:
                xbox = [ramin, ramin, ramax, ramax, ramin]
                ybox = [decmin, decmax, decmax, decmin, decmin]
                plt.plot(xbox, ybox, color = 'red')
            
        elif shapeType == 'ellipse':
            print("Ellipse.")
            
        elif shapeType == 'circle':
            if verbose:
                print("Circle.", split)
            
            radiusAS = split[2]
            if xy:
                radius = float(radiusAS[:-1])
            else:
                radius = float(radiusAS[:-1])/3600.0 # in degrees

            #cosfactor = np.cos(np.pi*deccentre/180.0)
            #radius = radius/cosfactor
            
            #print("The radius is ", radius, line)
            #exit()
            
            # cut!
            sep = np.sqrt((ra - racentre)**2 + (dec - deccentre)**2)

            bb = (sep < radius)
            
            # plot the circle
            if plot:
                phis=np.arange(0,6.28,0.01)
                xcircle = racentre + radius*np.cos(phis)
                ycircle = deccentre + radius*np.sin(phis)
                ax.plot(xcircle, ycircle, c='red',ls='-')
            
        elif shapeType == 'polygon':
            
            
            # For video (masked by Boris), and me...
            # First get all the corners into an array
            racorn = np.array(split[0::2])
            deccorn = np.array(split[1::2])

            if xy:
                x = racorn 
                y = deccorn
            else:
                
                c = SkyCoord(racorn, deccorn, unit = (u.hourangle, u.deg))
                x = c.ra.degree
                y = c.dec.degree

            #            print("Polygon", y)
            
            verbose = True
#            if (np.max(y) > -4.3) & (np.min(y) < -4.3):
#                verbose = True
#                print('Check, original ', x,y)

            
            # Make sure the order is right
            if np.argmax(y) > 0:
                value = np.argmax(y)
                x = np.roll(x, -value)
                y = np.roll(y, -value)

                if verbose:
                    print('new x, y = ', x, y)
                                
            if x[1] < x[-1]:
                # change direction to clockwise
                x = np.flip(x, axis = 0)
                y = np.flip(y, axis = 0)
                x = np.roll(x, 1)
                y = np.roll(y,1)
                if verbose:
                    print('new new x, y = ', x, y)
            
            # Section 1
            # define the y range for this section
            ymax = y[0]
            if y[1] > y[-1]:
                ymin = y[1]
            else:
                ymin = y[-1]
            
            # now get gradients of the lines
            mAB = (y[1]-y[0])/(x[1]-x[0])
            mAD = (y[0]-y[-1])/(x[0]-x[-1])
            mBC = (y[2]-y[1])/(x[2]-x[1])
            mCD = (y[3]-y[2])/(x[3]-x[2])
            
            cAD = y[0]-mAD*x[0]
            cAB = y[1]-mAB*x[1]
            cBC = y[2]-mBC*x[2]
            cCD = y[3]-mCD*x[3]

            lineAD = (dec-cAD)/mAD
            lineAB = (dec-cAB)/mAB
            lineCD = (dec-cCD)/mCD
            lineBC =  (dec-cBC)/mBC
            
            # first test
            condy = (dec > ymin) & (dec <= ymax)
            cond = condy & (ra < lineAB) & (ra > lineAD)
            
            # second section
            if y[1] > y[-1]:
                ymax = y[1]

                if y[2] > y[-1]:
                    ymin = y[2]
                else:
                    ymin = y[-1]
                
                # need BC gradient
                lineleft = lineAD
                lineright =lineBC
            else:
                ymax = y[-1]
                ymin = y[1]
                lineleft = lineCD
                lineright = lineAB

            condy = (dec > ymin) & (dec <= ymax)
            condtwo = condy & (ra < lineright) & (ra > lineleft)
            
            # final section
            # final section, three options?
            # C bottom:
            if (y[2] < y[-1]) & (y[2] < y[1]):
                ymin = y[2]
                lineleft = lineCD
                lineright = lineBC
                if y[1] > y[-1]:
                    ymax = y[-1]
                else:
                    ymax = y[1]
                    
            elif (y[1] < y[2]) & (y[1] < y[-1]):
                # B bottom
                ymin = y[1]
                lineleft = lineBC
                lineright = lineAB
                ymax = y[2]
            elif (y[-1] < y[2]) & (y[-1] < y[1]):
                # D bottom
                ymin = y[-1]
                lineleft = lineAD
                lineright = lineCD
                ymax = y[2]
                
            
            condy = (dec > ymin) & (dec <= ymax)
            condthree = condy & (ra < lineright) & (ra > lineleft)
            
            bb = cond | condtwo | condthree
            #if np.sum(bb) < 1:
            #    print('Nothing inside...')
                
            if verbose:
                print('Exciting.')
                exit()
            
            # this gives inside of the polygon
            if plot:
                xfull = np.append(x, x[0])
                yfull = np.append(y, y[0])
                ax.plot(xfull, yfull, color = 'red', zorder = 10)
            
            
        else:
            print("Shape not recognised: ", shapeType)
    
        # mask for this shape
        if plot:
            ax.scatter(ra[bb], dec[bb], color = 'green', s= ms)    
        maskArray[bb] = 0

    print('Created mask with {0} bad bits.'.format(maskArray.size -np.sum(maskArray)))
    if plot:
        ax.set_rasterized(rasterize)
        ax.set_ylim([np.min(dec), np.max(dec)])
        ax.set_xlim([np.min(ra), np.max(ra)])
        #fig.savefig(plotname)
        #print('check plot, ', plotname)
        #exit()
    # invert to mask in next code
    if tokeep:
        return maskArray.astype(bool)
    else:
        return np.invert(maskArray.astype(bool))

def combine_cats(field, detFilt, catalogueTypes, apDiametersAS = [1.8, 2.0, 3.0, 4.0, 5.0], imageDir = '/raid/scratch/data/', requiredSpit = ['2.8', '2.8', '2.8', '2.8', '2.8'], IRACapDiametersAS = [2.8, 3.8, 5.8, 9.8, 11.6]):
    print("Running combine_cats: ", '\n')

    # basically combine everything that is available in the directory of this
    dirHere = imageDir + field + '/'
    # read in filters
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
        
    imagedata = Table.read(inputFile, format = 'ascii.commented_header')
    availableFilters = np.array(imagedata['Name'])
    print("The available filters are ", availableFilters)
    
    ############################################################
    
    catDir = imageDir + 'catalogues/'
    indiDir = catDir + 'COSMOS/det_{0}/'.format(detFilt)
    
    ## Base table, make a basic table from the detection catalogue
    seCatalogue = indiDir + 'd' + detFilt + '_m' + detFilt + '.fits'
    print("seCatalogue: ", seCatalogue)
    tbdet = Table.read(seCatalogue, hdu = 2)
    columns = tbdet.colnames
    
    # loop through the required catalogue types
    for ci, catType in enumerate(catalogueTypes):
        
        if num_there(catType):
            typeString = catType + 'as'
            requiredMagKey = 'MAG_APER'
            
        else:
            #parts = catType.split('_')
            typeString = catType
            requiredMagKey = catType

        # get the det index
        detIndex = np.where((detFilt == availableFilters))
        detIndex = detIndex[0]


        # define output name
        # determine the flux key
        parts = requiredMagKey.split('_')
        if parts[1] == 'APER':
            requiredFluxKey = 'FLUX_' + parts[1]
            requiredSpitHere = requiredSpit[ci]
            if requiredSpitHere == 'NONE':
                outputName = indiDir + '{0}_DR2_UNMASKED_{1}_{2}.fits'.format(field, detFilt, typeString)
            else:
                outputName = indiDir + '{0}_DR2_UNMASKED_{1}_{2}_IRAC{3}as.fits'.format(field, detFilt, typeString, requiredSpitHere)

        else:
            requiredFluxKey = 'NONE'
            # e.g. for flux_radius, fwhm etc
            outputName = indiDir + '{0}_DR2_UNMASKED_{1}_{2}.fits'.format(field, detFilt, typeString)

            
        # sort out the base table for this one
        reqColumns = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', \
                          'DELTA_J2000', \
                          requiredMagKey, \
                          requiredFluxKey]
        
        # Create a base table, not fluxes/mags just coords
        baseTable = tbdet[reqColumns[0:5]]
        
        # Now update the column names
        baseTable['NUMBER'].name = 'ID'
        baseTable['ALPHA_J2000'].name = 'RA'
        baseTable['DELTA_J2000'].name = 'DEC'
        
        # Now extract the mags
        allMags = np.array(tbdet[requiredMagKey])
        magTable = Table()

        if requiredFluxKey != 'NONE':
            allFluxes = np.array(tbdet[requiredFluxKey])    
            fluxTable = Table()

        # split up the aper mags
        if parts[1] == 'APER':
            # match the aperture!
            ai = (np.array(apDiametersAS) == float(catType))
            magTable[detFilt] = allMags[:,ai]
            fluxTable['flux_' + detFilt] = allFluxes[:,ai]
        else:
            magTable[detFilt] = allMags

        ######################################################
        # loop through the available files!
        for fi, filterName in enumerate(availableFilters):

            seCatalogue = indiDir + 'd' + detFilt + '_m' + filterName + '.fits'
            exists = os.path.isfile(seCatalogue)

            # get the required aperture
            if parts[1] == 'APER':
                if filterName[0:2] == 'ch':
                    if requiredSpitHere == 'NONE':
                        ai = -1
                        print('No spitzer being included')
                    else:
                        ai = (np.array(IRACapDiametersAS) == float(requiredSpitHere))
                else:
                    ai = (np.array(apDiametersAS) == float(catType))
                    
            
            if filterName != detFilt and exists:
                
                print("Reading ", seCatalogue)
                tb = Table.read(seCatalogue, hdu = 2)
                columns = tb.colnames

                
                ## Do a check that the columns match up
                #offset = np.array(tb['X_IMAGE']) - np.array(baseTable['X_IMAGE'])
                check = np.array_equal(np.array(tb['X_IMAGE']), np.array(baseTable['X_IMAGE']))
                if check == False:
                    print("The catalogues are not the same length!!")
                    offset = np.abs(np.array(tb['X_IMAGE']) - np.array(baseTable['X_IMAGE']))
                    ind = np.where(offset > 0)
                    print(offset[ind])
                    print(tb['X_IMAGE'])
                    print(baseTable['X_IMAGE'])
                    print("EXITING")
                    exit()
                
                ## extract the required mag and flux
                # If auto doesn't exist, set as -99.0
                jj = (requiredMagKey == np.array(columns))
                if np.any(jj) == False:
                    allMeasuredMag = np.zeros(len(tbdet))
                    allMeasuredMag[:] = -99.0
                else:
                    allMeasuredMag = tb[requiredMagKey]
                    
                if requiredFluxKey != 'NONE':
                    measuredMag = Column(allMeasuredMag[:,ai], name = filterName)
                    # also the fluxes
                    allMeasuredFlux = tb[requiredFluxKey]
                    # HOLLY 11/03/2025: changed from measuredFlux = MaskedColumn(allMeasuredFlux[:,ai], name = 'flux_' +
                    measuredFlux = MaskedColumn(allMeasuredFlux[:,ai], name = 'flux_' + filterName)
                else:    
                    measuredMag = Column(allMeasuredMag, name = filterName)
                    
                ## add this to the new table, either before or after the mag,
                ## depending on reference to the other bands
                
                if fi < detIndex:
                    magTable.add_column(measuredMag, index = fi)
                    if parts[1] == 'APER':
                        fluxTable.add_column(measuredFlux, index = fi)
                else:
                    magTable.add_column(measuredMag)
                    if parts[1] == 'APER':
                        fluxTable.add_column(measuredFlux)

            elif filterName != detFilt:
                print("MISSING filter {0}".format(filterName))
                        
            # now combine the different parts of the table!
            if requiredFluxKey != 'NONE':
                # haven't saved mag auto flux so this code is taking aper flux
                # then just not saving it to file:
                finalTable = hstack([baseTable, magTable, fluxTable])
            else:
                finalTable = hstack([baseTable, magTable])


            
        finalTable.write(outputName, overwrite = True)
        print("The output catalogue is ", outputName)
        
    return

def num_there(s):
    return any(i.isdigit() for i in s)

def check_all_files(field, detFilt, imageDir = baseDir + 'data/'):
    """ (Holly 03/03/25)
    Checks .lis files exists and extracts data. Gets vailable filters from .lis header, gives this as user feedback. For each filter in availableFilters, checks if a .fits catalogue exists. If yes, gives user feedback by reporting the size of the catalogue file. 
    """
    print("Running check_all_files: ", '\n')
    catDir = imageDir + 'catalogues/{0}/det_{1}/'.format(field, detFilt)

    dirHere = imageDir + field + '/'
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
        
    imagedata = Table.read(inputFile, format = 'ascii.commented_header')
    availableFilters = np.array(imagedata['Name'])
    print("The available filters are ", availableFilters)

    for fi, filt in enumerate(availableFilters):

        expectedName = catDir + 'd{0}_m{1}.fits'.format(detFilt, filt)
        # get the length
        
        if os.path.isfile(expectedName):
            
            b = os.path.getsize(expectedName)
            if b > 10000:
                tt = Table.read(expectedName, hdu = 2)
                print('{0}\tEXISTS\t{1}'.format(filt, len(tt)))
            else:
                print('{0}\tTOO SMALL\t{1}'.format(filt, expectedName))
        else:
            print('{0}\tMISSING\t{1}'.format(filt, expectedName))
            
    print("check_all_files: finished running")        
    return

def run_source_extractor(field, detFilt, apDiameterAS, queue = 'none', reqFilters = ['all'], excludeFilters = ['ch1o', 'ch2o', 'ugriy', 'GRI', 'GR'], interactive = True, overwrite = False, imageDir = baseDir + 'data/', memory = 10.0, IRACapDiametersAS = [2.8, 3.8, 5.8, 9.8, 11.6]):
    print("Running run_source_extractor: ", '\n')
    # Note: removed 'YJ' from excludeFilters.
    # Note: removed 'JH', 'HKs' from excludeFilters.

    # standard inputs etc

    inputSex = baseDir + 'data/'+'bertin_config/video_mine.sex'
    #os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/2.25.0/share/sextractor'
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor-2.25.0/config'
    # define the directories
    ## Extract the information here
    dirHere = imageDir + field + '/'
    print("The directory is ", dirHere)
    
    catDir = imageDir + 'catalogues/{0}/'.format(field)
    indiDir = catDir + 'det_{0}/'.format(detFilt)
    
    if os.path.isdir(catDir) and os.path.isdir(indiDir):
        print("Directories created.")
    else:
        os.system('mkdir ' + catDir)
        os.system('mkdir ' + indiDir)

    
    ###################################################
    ################# Sort filters out ################
    ## Find the required images and their locations
    ## Read in the images.lis file here
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
        
    imagedata = Table.read(inputFile, format = 'ascii.commented_header')
    availableFilters = np.array(imagedata['Name'])
    print("The available filters are ", availableFilters)
    #for i, ii in enumerate(availableFilters):
    #    print(ii)
    # now remove the undesirable filters!
    
    excludeFilters = np.array(excludeFilters)
    
    if excludeFilters.size > 1:
        keep = np.ones(availableFilters.size, dtype=bool)
        for bf, badFilt in enumerate(excludeFilters):
            ii = np.where((availableFilters == badFilt) & (badFilt != detFilt))
            keep[ii] = False
            print("Excluding filter ", badFilt)
            
        imagedata = imagedata[keep]
        
    availableFilters = imagedata['Name']

    if reqFilters[0] == 'all':
        reqFilters = np.copy(availableFilters)

    else:
        print("Only running for a selection of filters ")
        availableFilters = np.array(reqFilters)

    # check!
    if interactive:
        print("I will now run SE for {0} filters in field {1}, on queue = {2}.".format(availableFilters.size, field, queue))
        print("The filters are ", availableFilters)
        #input("Press enter to continue.")
        
    # check the detection filter exists in catalogue
    detIndex = (availableFilters == detFilt)
    print('detIndex', detIndex)
    
    if np.sum(detIndex) < 1:
        print("The input detection filter does not exist.")
        print("Input = ", detFilt, " out of possible filters: ", availableFilters)
        
    else:
        detIndex = detIndex[0]
        
        
    # Change the deblending...
    
    for fi, filterName in enumerate(availableFilters):

        if filterName[0:2] == 'ch':
            apDiameterASuse = IRACapDiametersAS
            print('Changing apertures for filter {0}, to {1}'.format(filterName, apDiameterASuse))
            
        else:
            apDiameterASuse = apDiameterAS
            
        # Now finally run sextractor!

        outputCat = run_se(detFilt, filterName, imagedata, apDiameterASuse, inputSex, dirHere, queue = queue, field = field, outputDir = indiDir, overwrite = overwrite, memory = memory)

    return

def run_se(detFilt, filterName, imagedata, apDiametersAS, inputSex, dirHere, assoc_file = 'none', queue = 'none', field = 'NONE', outputDir = '', overwrite = False, memory = 10.0, verbose = True):
    print("Running run_se: ", '\n')
 
    from astropy.io import fits

    # extract the images etc for the detection filter
    availFilt = imagedata['Name']
    dd = (availFilt == detFilt)


    detImage = imagedata['Image'][dd][0]

    # and the measurement filter
    mm = (availFilt == filterName)
    if np.sum(mm) == 0:
        print(f"Warning: Filter '{filterName}' not found in available filters!")
        return None  # Exit func
    measImage = imagedata['Image'][mm][0]
    
    ## Extract directory for detection image
    if imagedata['directory'][dd] != 'here':
        detimageDir = imagedata['directory'][dd][0]
    else:
        detimageDir = dirHere
        
    # Extract the directory
    if imagedata['directory'][mm] != 'here':
        measimageDir = imagedata['directory'][mm][0]
    else:
        measimageDir = dirHere
        
    ## extract the pixel size
    hdulist = fits.open(detimageDir + detImage)
    #hdulist = fits.open(detimageDir + 'UVISTA_JH_DR6.fits') # holly
    imageHeader = hdulist[0].header
    cdone_o = -3600.0*imageHeader['CD1_1']
    pixScale = round(cdone_o, 5)
        
    ############## Aperture sizes ###################
    ## Do something for different for Spitzer, if required
    #if spitzerApertures[0] > 0 and (filterName[0:3] == 'ch1' or filterName[0:3] == 'ch2'#):
    #    print "I am using different apertures for Spitzer. "
    #    apDiametersAS = spitzerApertures
    #else:
    ## Use the default apertures
    #apDiametersAS = apDiametersAS_optical
    
    # and create my aperture string
    apStringPix = str(apDiametersAS[0]/pixScale) 
    for ai in range(1, len(apDiametersAS)):
        apStringPix = apStringPix + ',' + str(apDiametersAS[ai]/pixScale)        
    
    outputCat = outputDir + 'd' + detFilt + '_m' + filterName + '.fits'
    print("outputCat will be: ", str(outputCat), '\n')

    print('Here', filterName[0:2])
    if filterName[0:2] == 'ch':
        print('Running my own phot for filterName {0}'.format(filterName))
        # Read in the det filter catalogue
        inputCat = outputDir + 'd' + detFilt + '_m' + detFilt + '.fits'
        if field[0:2] == 'CD':
            # use the background subtracted images
            measimageDir = '../depths/{0}/images/'.format(field)
            measImage = '{0}_bgsub.fits'.format(filterName)
            print('Using instead the bg subtracted image ', measImage)
            
        command = 'python3 stupid_phot.py {0} {1} {2} {3} {4}'.format(measimageDir + measImage, inputCat, apStringPix, str(imagedata['zeropoint'][mm][0]), outputCat)
        
        if queue != 'none':
            # launch on glamdring queue
            # make a file
            shell = 'tmp_myfile_{0}_{1}_{2}.sh'.format(detFilt, filterName, field)

            f = open(shell, 'w')
            f.write('#!/bin/bash\n')
            f.write(command)
            f.close()
            # make executable
            os.system('chmod +x ' + shell)

            # launch sextractor!
            os.system("addqueue -c '{0}_{1}' -q {2} -m {4} -d ./{3}".format(detFilt, filterName, queue, shell, memory))
            
        else:
            # just run here
            do_phot_sex(measimageDir + measImage, inputCat, apStringPix, str(imagedata['zeropoint'][mm][0]), outputCat)
        print("run_se finished", '\n')

        return outputCat
    
    # Run Source Extractor
    keywords = ' -CATALOG_TYPE FITS_1.0 -CATALOG_NAME ' + outputCat + \
               ' -MAG_ZEROPOINT '+ str(imagedata['zeropoint'][mm][0]) + ' ' + \
               ' -PHOT_APERTURES ' + apStringPix

    # set the background parameters
   # if background == 'DEFAULT':
   #     bkgkeywords = ''
   # elif background == 'NO':
   #     bkgkeywords = ' -BACK_TYPE MANUAL -BACK_VALUE 0.0 '
   #     print('Setting a constant background.')
   # else:
   #     print('ERROR unknown background flag in run_se')
   #     exit()
    

    # Holly
    #bkgkeywords = ' -BACK_PEARSON "' + 3.5 \ #v2.28 or later unfortunately, we have v2.25.0 (2022-06-29)
    
    # if we have the detection image, change this slightly!
    if filterName == detFilt:

        imageString = measimageDir + measImage
        print("filterName == detFilt: ", filterName, '==', detFilt, '\n')
        keywordsWeight = ' -WEIGHT_TYPE ' + imagedata['Wht_type'][mm][0] + \
                         ' -WEIGHT_IMAGE ' +  measimageDir + imagedata['Weight'][mm][0]
    else:
        print("filterName != detFilt: ", filterName, '!=', detFilt, '\n')
        ## for dual mode do Detection image, then measurement image
        imageString = '"' + detimageDir + detImage + ',' + \
                      measimageDir + measImage + '"'
        keywordsWeight = ' -WEIGHT_TYPE "' + imagedata['Wht_type'][dd][0] \
                         + ',' + \
                         imagedata['Wht_type'][mm][0] +  '" -WEIGHT_IMAGE "' + \
                         detimageDir + imagedata['Weight'][dd][0] + ',' + \
                         measimageDir + imagedata['Weight'][mm][0] + '"'
            
    if field == 'COSMOS':
        ## change deblending...
        extraKeywords = ' -DEBLEND_MINCONT	0.00001'
    else:
        extraKeywords = ''

    extraKeywords = extraKeywords
    
    # Now actually run the catalogues through (if they don't exist)
    fileExists = os.path.isfile(outputCat)

    # check the size of the image
    if fileExists:
        b = os.path.getsize(outputCat)
        if b < 10000:
            fileExists = False
            print("The file size is ", b)
            
    if fileExists and overwrite == False:
        print("Catalogue for {0},{1} exists, not overwriting.".format(detFilt, filterName))
        
    else:
        # run sextractor!
        command = '/usr/local/sextractor/bin/sex ' + imageString +' -c ' + inputSex + keywords + keywordsWeight + extraKeywords
        #ßcommand = imageString +' -c ' + inputSex + keywords + keywordsWeight + extraKeywords
        if verbose:
            print("command: ", command == '/usr/local/sextractor/bin/sex ' + imageString +' -c ' + inputSex + keywords + keywordsWeight + extraKeywords)
        
        if queue != 'none':
            # launch on glamdring queue
            # make a file
            print("queue != 'none' = yes")
            shell = 'tmp_myfile_{0}_{1}_{2}.sh'.format(detFilt, filterName, field)
            f = open(shell, 'w')
            f.write('#!/bin/bash\n')
            f.write(command)
            f.close()
            # make executable
            os.system('chmod +x ' + shell)
            
            # launch sextractor!
            os.system("addqueue -c '{0}_{1}' -q {2} -m {4} -d ./{3}".format(detFilt, filterName, queue, shell, memory))
        else:
            # just run here
            os.system(command)
    print("run_se finished", '\n')
    return outputCat

def spawn_flux_errors(
    detectionFilters, fields, requiredCats, requiredSpit,
    catString='DR2_MASKVISTA_', queue='none', memory=6, n_jobs=12, short_test=False):

    n_jobs = len(detectionFilters)
    jobs = []  # list to hold flux_errors jobs for parallel execution

    for ff, fieldName in enumerate(fields):
        print('#############################################')
        print("Fluxes/errors for field ", fieldName)

        for df, detectionFilt in enumerate(detectionFilters):
            print('For det filter ', detectionFilt)

            for ai, apD in enumerate(requiredCats):
                spitstring = '_IRAC{0}as'.format(requiredSpit[ai])
                optstring = '{0}as'.format(apD)

                if requiredSpit[ai] == 'NONE':
                    spitstring = ''

                if apD != 'MAG_AUTO':
                    inputCat = (
                        '/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/catalogues/COSMOS/det_{0}/{1}_{4}{0}_{2}{3}.fits').format(detectionFilt, fieldName, optstring, spitstring, catString)

                    if os.path.isfile(inputCat):
                        if queue == 'none':
                            # Local parallel execution
                            if requiredSpit[ai] == 'NONE':
                                stringSpit = 'NONE'
                            else:
                                stringSpit = float(requiredSpit[ai])

                            jobs.append(
                                delayed(flux_errors)(
                                    inputCat, fieldName,
                                    IRACapertureAS=stringSpit,
                                    outputType='fits',
                                    verbose=False,
                                    inputDirProvided=True,
                                    removeCoadds=True))

                        else:
                            # Queue submission mode
                            print('Spawning in the queue ')
                            tmpName = 'tmp_myfile_{0}_{1}_{2}{3}'.format(
                                detectionFilt, fieldName, optstring, spitstring
                            )
                            with open(tmpName + '.sh', 'w') as f:
                                f.write('#!/bin/bash\n')
                                if requiredSpit[ai] == 'NONE':
                                    spitstringfloat = 'NONE'
                                else:
                                    spitstringfloat = float(requiredSpit[ai])
                                f.write(
                                    'python3 stupid_flux.py {0} {1} {2}\n'.format(
                                        inputCat, fieldName, spitstringfloat))

                            # make executable and submit
                            os.system('chmod +x ' + tmpName + '.sh')
                            command = "addqueue -c '{1}' -m {2} -q {0} -d ./{1}.sh".format(queue, tmpName, memory)
                            os.system(command)

                    else:
                        print('Input cat does not exist = ', inputCat)

    # Run local jobs in parallel if any
    if queue == 'none' and jobs:
        print(f"Running {len(jobs)} flux_errors jobs in parallel with n_jobs={n_jobs} at ", datetime.now().strftime("%d/%m %H:%M:%S"))
        results = Parallel(n_jobs=n_jobs, verbose=10)(jobs)
        return results

    return None
    
def flux_errors(inputTableFile, fieldName, dataDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/', overwrite = False, minError = 10.0, psfType = 'peak', minErrorIRAC = 10.0, IRACapertureAS = 2.8, depthType = 'grid', depths = [], provFilters = [], outputType = 'fits', verbose = False, inputDirProvided = False, removeCoadds = True, checkFile = 'none', microJyYes = True, xyName = ['X_IMAGE', 'Y_IMAGE'], badFilters = np.array(['NB118']), xmmflag = ['NONE'], countsKey = 'flux_', short_test=False):
    
    """ A function version of my flux errors code.
        
    ## Depth options are 'global', 'grid' or 'local' - where custom depths are calculated.

    ### NOTES FROM ROHAN
    ### 'rohangrid' does Rohan's custom grid depths.
    ### Also removed YJ from badFilters for my Y+J stack selection.
    ### Added psfType = 'grid' option to use PSFEx grid for the problematic HSC filters in CDFS.
    """

    ### NOTES FROM HOLLY
    ### Removed 'JH', 'HKs', 'HK' from badFilters for stacks
    
    psfType = 'peak'
    gridString = '_300_200'
    
    ####################### Import useful modules ###########################
    import astropy.io.fits as fits
    import os
    import time
    from astropy.table import Table, hstack, Column
    from hu_depth_codes import local_depths, grid_depths, grid_psf
    
    #################################################################
    # Read in the image info etc as before
    dirHere = dataDir + fieldName + '/'
    inputFile = dirHere + 'images.lis'
    
    # read in the input file
    # Check the input file exists
    if os.path.isfile(inputFile) == False:
        print("line 938")
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
        
    imagedata = Table.read(inputFile, format = 'ascii.commented_header')
    availableFilters = imagedata['Name']
    
    # Remove the coadds
    if removeCoadds:
        #print "Table before ", imagedata['Name']
        keep = np.ones(availableFilters.size, dtype=bool)
        #print keep
        for bf, badFilt in enumerate(badFilters):
            ii = np.where((availableFilters == badFilt))
            keep[ii] = False
        
        imagedata = imagedata[keep]
        print("After = ", imagedata['Name'])

    ###############################################################
    # Check that the required files exist
    # For both the depths and the psf correction
    
        
        
    ################################################################
    ## Extract the aperture size from the name
    #print "Hello", inputTableFile
    if depthType == 'input':
        apString = '1.8'
    else:
        beg = inputTableFile.find('as')
        if beg > 0:
            apString = inputTableFile[beg-3:beg]
        else:
            print('Cannot find aperture string!')
            print(inputTableFile)
            exit()
        
    apertureSize = float(apString)
    print("The aperture size is ", apertureSize)
    
#####################################################################
    
    if inputDirProvided:
        inputTableLocation = inputTableFile
        catDir = ''
    else:
        catDir = dirHere + 'catalogues/'
        inputTableLocation = catDir + inputTableFile
        
    #inputTableFile = '../catalogues/XMM1/det_HSC-G/XMM1_DR2_MASKED_HSC-G_1.8as_206testpsf.fits'
    print("Reading data from ", inputTableFile)
    inputTableLocation = inputTableFile
    
    tb = Table.read(inputTableLocation)

    if short_test:
        print('Testing on short cat!')
        tb = tb[np.random.randint(len(tb), size = (50))]
    
    #checkFile = '{0}_delete_check_random1000.fits'.format(fieldName)

    #checkFile = '22_08_HSC_masking_checkfile_HSCG.fits'
    
    print('Getting fluxes for {0} objects in {1}'.format(len(tb), inputTableLocation))
    
    if xmmflag[0] != 'NONE': # and np.any(h):
        print("I have an XMM flag, I want to duplicate the ch1/ch2 columns into ch1servs/ch2servs")
        
        if np.any(np.array(tb.colnames) == 'ch1servs') == False:
            servsOne = Column(tb['ch1'], name = 'ch1servs')
            tb.add_column(servsOne)
            
            servsTwo = Column(tb['ch2'], name = 'ch2servs')
            tb.add_column(servsTwo)
            
        # Also copy the fluxes
            servsOne = Column(tb['flux_ch1'], name = 'flux_ch1servs')
            tb.add_column(servsOne)
            
            servsTwo = Column(tb['flux_ch2'], name = 'flux_ch2servs')
            tb.add_column(servsTwo)
        
        #if np.any(np.array(tb.colnames) == 'ch1') == False:
        #    servsOne = Column(tb['ch1servs'], name = 'ch1')
        #    tb.add_column(servsOne)
        #    
        #    servsTwo = Column(tb['ch2servs'], name = 'ch2')
        #    tb.add_column(servsTwo)
        #    
        # Also copy the fluxes
         #   servsOne = Column(tb['flux_ch1servs'], name = 'flux_ch1')
         #   tb.add_column(servsOne)
         #   
         #   servsTwo = Column(tb['flux_ch2servs'], name = 'flux_ch2')
         #   tb.add_column(servsTwo)

    
    ## Do a check that we have all these filters...
    inputcolnames = np.array(tb.colnames)
    print("inputcolnames", inputcolnames)
    fullFilters = imagedata['Name']
    print("fullFilters", fullFilters)
    keep = np.ones(fullFilters.size, dtype=bool)
    keep[:] = False
    #inputcolnames = np.where(inputcolnames == "flux_Ks", "flux_K", inputcolnames)
    #inputcolnames = np.where(inputcolnames == "Ks", "K", inputcolnames)
    
    for fi, filt in enumerate(fullFilters):
        ii = (inputcolnames == filt)
        if np.any(ii):
            keep[fi] = True
            if keep[fi]:
                print("Keeping filter ", filt)
        else:
            print("Missing filter ", filt)
            
    imagedata = imagedata[keep]
    availableFilters = imagedata['Name']
#    print("the availableFilters are ", availableFilters)
    
    ## Define the catalogue names    
    #gg = inputTableFile.rfind('/')
    ee = inputTableFile.rfind('.fits')
    baseName = inputTableFile[:(ee)]
    
    
    # Create the output tables
    fluxmuJytb = Table()
    fluxCGStb = Table()
    errormuJytb = Table()
    errorCGStb = Table()
    depthsTable = Table()
    encFluxTable = Table()
    
    # check if ID other column exists...
    extra = 'IDother'
    if np.any(np.array(tb.colnames) == extra):
        speczzz = 'zspec'
        if np.any(np.array(tb.colnames) == speczzz):
            reqColumns = ['ID', extra, speczzz, xyName[0], xyName[1], 'RA', 'DEC']
        else:
            reqColumns = ['ID', extra, xyName[0], xyName[1], 'RA', 'DEC']
    elif depthType == 'input':
        reqColumns = ['ID']
    else:
        reqColumns = ['ID', xyName[0], xyName[1], 'RA', 'DEC']
        
    baseTable = tb[reqColumns]
    
    ##############################################################################
    # Define the aperture correction and local depth files
    if fieldName[0:3] == 'XMM':
        psfDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/' + 'psf/{0}/enclosedflux/'.format('XMM')
    elif fieldName[0:4] == 'CDFS':
        psfDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/' + 'psf/{0}/enclosedflux/'.format('CDFS')
    else:
        psfDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/' + 'psf/{0}/enclosedflux/'.format(fieldName)

    print('Enclosed flux is coming from directory:', psfDir)
    depthDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/depths/COSMOS/phot/' #TODO: check this dir is correct fo rfincal run
  
## read in enclosed flux stuff
#efFile = 'psf/final_enclosedFlux_Oct17.lis'
#enclosedFluxData = np.genfromtxt(efFile, names=True, dtype = None, deletechars='', skip_header = 1)
#print enclosedFluxData.dtype.names
#enclosedFluxAps = enclosedFluxData['aperture']
##print enclosedFluxAps
#api = np.where(float(apDiameterAS) == enclosedFluxAps)

###############################################################################
## Loop through the objects
## First convert to fluxes/local errors
## In both microJy and cgs
    if depthType != 'input':
    	x = tb[xyName[0]]
    	y = tb[xyName[1]]
    
# I need to loop through the filters
    for fi, filterName in enumerate(availableFilters):
        
        zeropoint = imagedata['zeropoint'][fi]
        
        if filterName[0:3] == 'ch1' or filterName[0:3] == 'ch2':
            apDiameterAS = IRACapertureAS
            if type(apDiameterAS) == float:
                apStringHere = '{0:.1f}'.format(apDiameterAS)
            else:
                apStringHere = apDiameterAS
                
            depthDiametersAS = np.array([2.8, 3.8, 5.8, 9.8, 11.6]) #2.8, 3.8, 8.2])
            minErrorHere = minErrorIRAC
            
        else:
            apDiameterAS = apertureSize
            apStringHere = apString
            depthDiametersAS = np.array([1.8, 2.0, 3.0, 4.0, 5.0])
            # this is only used for local depths
            minErrorHere = minError
            
        if depthType == 'local':
            
            print("Extracting local depths for filter ", filterName)
            # Calculate my own local depth from the aperture files
            tic = time.time()
            depthTableFile = depthDir + filterName + '_aperPhot.fits'
            inputTable = Table.read(depthTableFile)
            #print inputTable
            numApertures = 200
            toc = time.time()
            print("The time to read in the {0} file is {1:.4f}".format(depthTableFile, toc-tic))
            
        # cut this depending on the seg and wht values
            reqDepthInd = np.where((float(apDiameterAS) == depthDiametersAS))
            reqDepthInd = reqDepthInd[0]
            reqDepthInd = reqDepthInd[0]
            apStringLocal = '_' + str(reqDepthInd)
            print("The required index is ", reqDepthInd, " in ", depthDiametersAS)
            seg_sum = np.array(inputTable['SEG_flux' + apStringLocal])
            wht_sum = np.array(inputTable['WHT_flux' + apStringLocal])
            good_indicies = (seg_sum < 0.5)  & (wht_sum > 0.)
            print("There are ", sum(good_indicies), " good indicies.")
            good_indicies = np.where(good_indicies)
            
        ## This table has all the good apertures
            finalTable = inputTable[good_indicies]
            depthArray = local_depths(finalTable, apStringLocal, \
                                      x, y, numApertures,\
                                      zeropoint = zeropoint)
            
        if depthType[0:4] == 'grid':
            
            ## Read in a pre-prepared table of depth values
            ## find the closest
            #if "K" in filterName:
                #filterName = filterName.split('K')[0] + 'Ks'
            depthTableFile = depthDir + filterName +'_' + apStringHere + 'as_gridDepths{0}.fits'.format(gridString)
            if 'CFHT' in depthTableFile:
                continue
            print("####################### Filt = {0} GRID DEPTHS ################ ".format(filterName))
            tic = time.time()
            
            #if depthType == 'grid':
                
             #   depthTableFile = depthDir + filterName +'_' + apStringHere + 'as_gridDepths.fits'
            #elif depthType == 'grid19':
             #   print("Using new 2019 depths.")
            #depthTableFile = depthDir + filterName +'_' + apStringHere + 'as_gridDepths{0}.fits'.format(gridString)
            
            # Read in the depths table,  print("Input table is ", depthTableFile)
            print("depthTableFile: ", depthTableFile)
            inputTable = Table.read(depthTableFile)

            if np.isnan(np.sum(inputTable['depths'])):
                print('NAN in the input table of depths.')

            # Extract the depths at each input X/Y
            depthArray = grid_depths(inputTable, x, y)
            toc = time.time()
 #           print("The time to do the depths {:4f} seconds. ".format(toc-tic))
            
        if depthType == 'rohangrid':
            ## Read in Rohan's preprepared table of depth values.
            ## Find the closest coordinate and assign the corresponding depth.
            print("####################### Filt = {0} GRID DEPTHS ################ ".format(filterName))
            tic = time.time()
            
            ## _R = rohan's directory

            filterName_R = filterName
            
            # Deal with my grid naming convention
            if filterName_R.split('-')[0] == 'HSC':
                filterName_R = filterName.split('-')[1]
            # Check
            print(' MY FILTER NAME: ', filterName_R)

            # Depth directory
            depthDir_R = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/depths/{0}/grids/'.format(fieldName)

            # Depth file
            depthTableFile = depthDir_R + '{0}_{1}_2.0as_grid_depths.fits'.format(fieldName, filterName_R)

            # VIDEO depth directory is slightly different because my depth code is messy
            if (filterName == 'Y') | (filterName == 'J') | (filterName == 'H') | (filterName == 'Ks') | (filterName == 'YJ'):
                depthDir_R = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/depths/{0}/grids/VIDEO/'.format(fieldName)
                depthTableFile = depthDir_R + '{0}_{1}_2.0as_grid_depths.fits'.format(fieldName, filterName_R)

            # if Spitzer or other previous data, read Rebecca's grid.
            if (filterName[0:2] == 'ch') | (filterName[0:5]=='UKIRT') | (filterName=='SC-z') | (filterName[0:4]=='CFHT'):
                depthTableFile = depthDir + filterName +'_' + apStringHere + 'as_gridDepths{0}.fits'.format(gridString)

            # Check
            print('DEPTH GRID FILE IS: ', depthTableFile)

            # Read in the depths table,  print("Input table is ", depthTableFile)
            inputTable = Table.read(depthTableFile)
            print(inputTable)

            if np.isnan(np.sum(inputTable['depths'])):
                print('NAN in the input table of depths.')

            # Extract the depths at each input X/Y
            depthArray = grid_depths(inputTable, x, y)
            print(depthArray)
            toc = time.time()
            print("The time to do the depths {:4f} seconds. ".format(toc-tic))
            
            
        if depthType == 'global':
            
            ## use the global depth
            ## Read this in!
            depthFile = dirHere + 'depths.fits'
            t = Table.read(depthFile)
            names = np.array(t.colnames)
            apD = 1.8
            #if filterName[0:2] == 'ch':
            #    apD = 2.8
            
            ii = np.where(t['ApDiam'] == apD)
            #print t['ApDiam']
            #print "Extracting aperture ", apD#, " or ", t['ApDiam'][ii] #, ii
            ll = np.where(names == filterName)[0]
            if len(ll) > 0:
                
                depthHereCol = t[filterName][ii]
                depthHere =  depthHereCol.quantity[0]

                
                # deep strips
                if (filterName == 'Y') | (filterName == 'J') | (filterName == 'H') | (filterName == 'Ks'):
                    # deal with strips
                    minStrip = np.amin(depthHereCol[1:5])
                    maxStrip = np.amax(depthHereCol[1:5])
                    
                    # take the min
                    globalDepthHere = minStrip
                    
                else:
                    globalDepthHere = depthHere
                print("Filter {0} found!, depth = {1:.2f}".format(filterName, globalDepthHere))
                
            else:
                print("Error filter not found... ")
                print(filterName)
                exit()
                
            depthArray = np.zeros(x.size)
            depthArray[:] = globalDepthHere
             
        if depthType == 'input':
            
            # use the arrays provided!
            # find the filter name in the input
            #ii = (filterName == provFilters)
            #ii = [x is filterName.strip() for x in provFilters]
            #print ii
            
            ll = np.where(provFilters == filterName.strip())[0]
            #print "Hel ", ll, len(ll)
            if len(ll) < 1:   
            #if np.any(ii) == False:
                print("Error: using input depths, however depth for filter {0} does not exist!".format(filterName))
                
                print(filterName)
                print(provFilters)
                
                exit()
                
            depthArray = np.zeros(len(tb))
            depthArray[:] = depths[ll[0]]
            print("All fine adding depth for filter, {0}, depth = {1}".format(filterName, depths[ll[0]]))
            #exit()
            

        if verbose:
            ss = '{0:.2f}'.format(depthArray[0])
            for i in range(1,len(depthArray)):
                ss = ss +'\t'+ '{0:.2f}'.format(depthArray[i])
            print("The depths = ", ss)
            
    ################################################
    ## Also extract the enclosed flux, using the output of PSFEx
        apStringHere = apString
        
        # compare to IRAC mashby thing...
        if filterName[0:2] == 'ch':
     #       print(IRACapertureAS, type(IRACapertureAS))
            #apStringHere = '{0:.1f}'.format(IRACapertureAS)
            
            if type(IRACapertureAS) == float:
                apStringHere = '{0:.1f}'.format(IRACapertureAS)
            else:
                apStringHere = IRACapertureAS
            
        # Get the enclosed flux over the field
        if psfType == 'total':
            
            print("Not doing enclosed flux stuff")
            enclosedFlux = np.zeros(len(tb))
            enclosedFlux[:] = 100.0
            
        if psfType == 'peak':
            
            enclosedFlux = np.zeros(len(tb))
            
            if filterName[0:2] == 'ch':

                # Read in the servs values
                servsfile = dataDir + 'psf/SERVS/enclosedFlux.txt'
                servs = Table.read(servsfile, format = 'ascii.commented_header')
                rrrr = (servs['filter'] == filterName[0:3])
                value = servs[apStringHere][rrrr]
                enclosedFlux[:] = value*100.0
                print('Enclosed flux: Using SERVS ef of {0:.1f}, aperture {1}\n'.format(100.0*value.quantity[0], IRACapertureAS))
            else:
                # read in peak file
                psfTableFile = psfDir + filterName  + '_peak.txt'
                print("psfTableFile",psfTableFile)

                # Note from Rohan: reading in the median of the PSFs in CDFS1,2,3. So this will only work when running on CDFS.
                # I'll turn it off after using it!!
                #psfTableFile = psfDir + filterName  + '_peak_median.txt'
                if 'CFHT' in filterName:
                    continue
                encflux = Table.read(psfTableFile, format = 'ascii.commented_header')
                print("encflux: ",encflux)
                # extract the right apd
                ii = (encflux['apD'] > apertureSize-0.01) & (encflux['apD'] < apertureSize+0.01)

                if np.any(ii):

                    if (fieldName[0:3] == 'XMM') & (filterName[0:3] == 'HSC') & (filterName != 'HSC-NB0921'):

                        # I need to get a positional dependent value
                        blahh, hsc_circle = which_field(tb['RA'], tb['DEC'], hsc= True)
                        #print(hsc_circle)
                        print('Check1: ', blahh, hsc_circle)
                        # now pick the correct enclosed flux for each part of the circle!
                        one = (hsc_circle == '1')
                        if np.any(one):
                            enclosedFlux[one] = encflux['ef_XMM1'][ii][0]*100.0
                        print('Check2: ', enclosedFlux)
                        three = (hsc_circle == '3')
                        if np.any(three):
                            enclosedFlux[three] = encflux['ef_XMM3'][ii][0]*100.0
                        print('Check3: ', enclosedFlux)
                        twothree = (hsc_circle == '2U3') | (hsc_circle == '2L3') | (hsc_circle == '2U') | (hsc_circle == '2L') | (hsc_circle == '2L2U') | (hsc_circle == '2L2U3')
                        if np.any(twothree):
                            enclosedFlux[twothree] = encflux['ef_XMM3'][ii][0]*100.0
                        print('Check4: ', enclosedFlux)
                        # finally the ones in the overlap between 2 and 1
                        onetwo = (hsc_circle == '12U') | (hsc_circle == '12L') | (hsc_circle == '12L2U')
                        if np.any(onetwo):
                            enclosedFlux[onetwo] = (encflux['ef_XMM1'][ii][0]*100.0 + encflux['ef_XMM3'][ii][0]*100.0)/2.0
                        print('Check5: ', enclosedFlux)
                    else:

                        value = encflux['ef'][ii][0]*100.0
                        enclosedFlux[:] =  value

                        print('Enclosed flux: Using peak ef = {0:.1f}, aperture {1}\n'.format(value, apertureSize, encflux['apD'][ii][0]))
                    
                else:
                    print('Error, the aperture for this catalogue is not available.')
                    
        # Adding bit that uses PSFEx map
        elif psfType == 'grid':

            if filterName[0:2] == 'ch':

                # Read in the servs values
                servsfile = dataDir + 'psf/SERVS/enclosedFlux.txt'
                servs = Table.read(servsfile, format = 'ascii.commented_header')
                rrrr = (servs['filter'] == filterName[0:3])
                value = servs[apStringHere][rrrr]
                enclosedFlux[:] = value*100.0
                print('Enclosed flux: Using SERVS ef of {0:.1f}, aperture {1}\n'.format(100.0*value.quantity[0], IRACapertureAS))

            else:

                # Read in the psf table.
                psfTableFile = Table.read(dataDir + 'psf/{0}/enclosedflux/{1}.fits'.format(fieldName, filterName))

                # Run modified version of Rebecca's grid_depths on the psf, which finds the nearest measurement from PSF grid.
                enclosedFlux = grid_psf(psfTableFile, x, y)

        else:
            
            if filterName[0:2] == 'ch':

                enclosedFlux = np.zeros(len(tb))

                # Read in the servs values
                servsfile = dataDir + 'psf/SERVS/enclosedFlux.txt'
                servs = Table.read(servsfile, format = 'ascii.commented_header')
                rrrr = (servs['filter'] == filterName[0:3])
                value = servs['{0:.1f}'.format(IRACapertureAS)][rrrr]
                enclosedFlux[:] = value
                print('Using SERVS enclosed flux for Spitzer/IRAC', value)
                
            else:
                 
                psfTableFile = psfDir + filterName  + '.fits'
                psfTable = Table.read(psfTableFile)
                
                print("Extracting the enclosed flux... ")
                tic = time.time()
                enclosedFlux = aperture_correction(psfTable, apStringHere, x, y)
                toc = time.time()
                print("The time to do the enclosed flux is {:4f} seconds. ".format(toc-tic))
        print("The enclosed flux is ", enclosedFlux)
            
        #if filterName == 'Y':
        #    enclosedFlux[:] = 67.2
        #if filterName == 'J':
        #    enclosedFlux[:] = 71.2
        #if filterName == 'H':
        #    enclosedFlux[:] = 75.4
            
            
        #print "Filter = ", filterName, " depth = ", ss
        if verbose:
            ss = '{0:.2f}'.format(enclosedFlux[0])
            for i in range(1,len(enclosedFlux)):
                ss = ss +'\t'+ '{0:.2f}'.format(enclosedFlux[i])
        #print "Enclosed flux = ", ss
        #exit()
        
        
            #print "The enclosed flux = ", ss
            
    ################################################
        fluxCounts = tb[countsKey + filterName]
        
       
        ## Check that the object is covered by this data set, using the flag
    ## in the magnitudes
        valuesToMask = np.where(tb[filterName] < -90.0)
        # if magnitude is measured as <-90, mask over it 
        
    ## Convert to fluxes and errors in good units!
        if microJyYes:
            fluxMicroJy, errorMicroJy = flux_conv(fluxCounts, depthArray, \
                                                      zeropoint, microJy = True, \
                                                      enclosedFlux = enclosedFlux, \
                                                      minErrorPercent = minErrorHere)
            fluxMicroJy[valuesToMask[0]] = -99.0 
            # assign a mask value to measured magnitudes <-90
            errorMicroJy[valuesToMask[0]] = -99.0
            
    ## add to table
            fluxMicroJy_col = Column(fluxMicroJy, name = 'flux_' + filterName)
            errorMicroJy_col = Column(errorMicroJy, name = 'err_' + filterName)
            fluxmuJytb.add_column(fluxMicroJy_col)
            errormuJytb.add_column(errorMicroJy_col)
            
            finalTable = hstack([baseTable, fluxmuJytb, errormuJytb])

        #print "The final fluxes and errors in muJy are "
        #print fluxMicroJy_col, errorMicroJy_col
        
    ## and in CGS, this is for le Phare mainly
        ## now convert to cgs
            
        tic = time.time()
        fluxCGS, errorCGS = flux_conv(fluxCounts, depthArray, \
                                          zeropoint, cgs = True, \
                                          enclosedFlux = enclosedFlux, \
                                          minErrorPercent = minErrorHere)
        toc = time.time()
#        print("The time to do the flux conversion is {:4f} seconds. ".format(toc-tic))
        
        #if filterName == 'Y':
        #    
        #    print zeropoint
        #    print depthArray[0:5], enclosedFlux[0:5]
        #    print "Flux in ", fluxCounts[0:5]
        #    print "Flux out ", fluxCGS[0:5]
        #    print "Mag in ",  -2.5*np.log10(fluxCounts[0:5]) + zeropoint
        #    print "Mag out ", -2.5*np.log10(fluxCGS[0:5]) -48.6
         #   
        #    exit()
        
        fluxCGS[valuesToMask[0]] = -99.0
        errorCGS[valuesToMask[0]] = -99.0
        
        
        fluxCGS_col = Column(fluxCGS, name = 'flux_' + filterName)
        errorCGS_col = Column(errorCGS, name = 'err_' + filterName)
        fluxCGStb.add_column(fluxCGS_col)
        errorCGStb.add_column(errorCGS_col)
        
    ## If I am creating a checkfile, add columns to this with the 
        
        ## 1) local depths
        ## 2) the enclosed fluxes
        if checkFile != 'none':
            
            encFluxCol = Column(enclosedFlux, name = 'ef_' + filterName)
            depthCol = Column(depthArray, name = 'ld_' + filterName)
            depthsTable.add_column(depthCol)
            encFluxTable.add_column(encFluxCol)
            
## Finally combine into one big table
## and save!  Set the format here
        
    finalTableCGS = hstack([baseTable, fluxCGStb, errorCGStb])
    
    if xmmflag[0] != 'NONE':
        
        print("Editing the output file because of XMM flag.")
        servs = (xmmflag == "SERVS")
        
        if np.any(np.array(finalTableCGS.colnames) == 'flux_ch1') == False:
            # rename the servs columns
            finalTableCGS.rename_column('flux_ch1servs', 'flux_ch1')
            finalTableCGS.rename_column('flux_ch2servs', 'flux_ch2')
            finalTableCGS.rename_column('err_ch1servs', 'err_ch1')
            finalTableCGS.rename_column('err_ch2servs', 'err_ch2')
            print("Renaming columns here.")
            
        elif np.any(servs):
            finalTableCGS['flux_ch1'][servs] = finalTableCGS['flux_ch1servs'][servs]
            finalTableCGS['flux_ch2'][servs] = finalTableCGS['flux_ch2servs'][servs]
            finalTableCGS['err_ch1'][servs] = finalTableCGS['err_ch1servs'][servs]
            finalTableCGS['err_ch2'][servs] = finalTableCGS['err_ch2servs'][servs]
            print("I have edited the flux/error of {0} objects to SERVS.".format(np.sum(servs)))
            # also remove the servs columns
            finalTableCGS.remove_column('flux_ch1servs')
            finalTableCGS.remove_column('flux_ch2servs')
            finalTableCGS.remove_column('err_ch1servs')
            finalTableCGS.remove_column('err_ch2servs')
            
    if checkFile != 'none':
        checkTable = hstack([baseTable, depthsTable, encFluxTable])

        # also make a column of the which field
        
        newcol = Column(hsc_circle.astype(str), name = 'HSC_FIELD')
        checkTable.add_column(newcol)
        
        checkTable.write(checkFile, overwrite = True)
        print("The check table of depths and enclosed flux has been saved to ")
        print(checkFile)
        
    
    if outputType == 'fits':
        
        if microJyYes:
            fluxErrorsMicroJy = catDir + baseName + '_microJy' + '.fits'
            finalTable.write(fluxErrorsMicroJy, overwrite = True)
            print("The mJy catalogue is ", fluxErrorsMicroJy)
            
        if depthType == "local":
            fluxErrorsCGS = catDir + baseName + '_local_cgs' + '.fits'
        elif depthType == 'global':
            fluxErrorsCGS = catDir + baseName + '_global_cgs' + '.fits'
        else:
            fluxErrorsCGS = catDir + baseName + '_cgs' + '.fits'
            
        finalTableCGS.write(fluxErrorsCGS, overwrite = True)

    if outputType == 'ascii':
        
        if microJyYes:
            fluxErrorsMicroJy = catDir + baseName + '_microJy' + '.cat'
            finalTable.write(fluxErrorsMicroJy, format='ascii', overwrite = True)
            
        fluxErrorsCGS = catDir + baseName + '_cgs' + '.cat'
        finalTableCGS.write(fluxErrorsCGS, format='ascii', overwrite = True)
    
    print("The CGS catalogue is ", fluxErrorsCGS)
    
    if microJyYes:
        return fluxErrorsMicroJy
    else:
        return fluxErrorsCGS

    return 1

def return_instrips(xra, ydec, radec = False, pureGap = False, region = 'allstrips', dr5 = True, dr6 = True):
    
    ''' Code to return the True/False values if we are in the deep strips or not
    
    pureGap = True return coordinates really far from the strips, to get a clean measure of teh gaps.
    
    Created Mon 27th Nov 2017'''
    
    import numpy as np
    from astropy import wcs
    from astropy.io import fits
    
    ## Read in the header information for UltraVISTA
    refImage = baseDir + '/data/COSMOS/UVISTA_J_dr5_rc1.fits'

    if dr6:
        refImage = baseDir + '/data/COSMOS/UVISTA_J_DR6.fits'
        #xra = xra*2.0
        #ydec = ydec*2.0
        
    #refImage = '/users/bowlerr/vardygroupshare/hsc/final_mosaics/COSMOS_DR3/UVISTA_J_dr3_v5.fits'
    hdu = fits.open(refImage)
    w = wcs.WCS(hdu[0].header)
    
    #stripFile = '/users/bowlerr/vardygroupshare/hsc/final_mosaics/COSMOS/strip_definitions.txt'
    ## Load in the region file for UltraVISTA
    #with open(stripFile) as f:
    #    lines = f.readlines()
    #    for line in lines:
    #        myarray = np.fromstring(line, dtype = float, sep = ',')
    #        print myarray


    if dr5:
        # new strips defined by Bo in May 2022
        # dec range is the same for all, as before
        strip1_ra = [150.569300, 150.435925]
        strip1extra_ra = [150.533333, 150.435925] # this is the bad corner region of the detector
        strip1_dec = [2.7075, 1.89]# different due to the bad corner
        strip1extra_dec = [1.89, 1.7191] # the bottom bit
        strip2_ra = [150.202217, 150.067925]
        strip3_ra = [149.833884, 149.701300]
        strip4_ra = [149.467842, 149.332752]
        strips_dec = [2.7075, 1.7191]

        # for the gaps, deep part
        gap1_ra = [150.750724, 150.615634]
        gap2_ra = [150.386342, 150.249842]
        gap3_ra = [150.019592, 149.881842]
        gap4_ra = [149.651134, 149.515092]

        # these leave the region between empty, so the edges depend on what we are
        # using this for
        
    else:
        
        ## split into ra and dec
        ## Just write this out myself!!
        strip1_ra = [150.59, 150.41]
        strip1_dec = [2.72, 1.99]
        strip1_extra_ra = [150.52, 150.41]
        strip1_extra_dec = [1.99, 1.7] # changed 15/10/2018
        
        strip2_ra = [150.22, 150.05]
        strips_dec = [2.72, 1.7] ## same for all!
        
        strip3_ra = [149.85, 149.68]
        strip4_ra = [149.49, 149.32]
        
        gaps_y = [2500, 26500]
        gap1_x = [1000, 4200.0]
        gap1_y = [9500, 26500]
        gap2_x = [9800.0, 13100.0]
        gap3_x = [18600, 21900]
        gap4_x = [27450.0, 30700.0]
    
        # do this again with the new strips
        edges = np.genfromtxt('/mnt/vardy/vardygroupshare/data/codes/area/ud_d_cosmos.txt', dtype = None, names = True)
    
    
    #edgex, edgey = w.wcs_world2pix(strip1_ra, strip1_dec, 1)
    #print '{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}'.format(edgex[0], edgex[1], edgey[0], edgey[1])
    #edgex, edgey = w.wcs_world2pix(strip2_ra, strips_dec, 1)
    #print '{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}'.format(edgex[0], edgex[1], edgey[0], edgey[1])

    #edgex, edgey = w.wcs_world2pix(strip3_ra, strips_dec, 1)
    #print '{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}'.format(edgex[0], edgex[1], edgey[0], edgey[1])

    #edgex, edgey = w.wcs_world2pix(strip4_ra, strips_dec, 1)
    #print '{0:.0f}\t{1:.0f}\t{2:.0f}\t{3:.0f}'.format(edgex[0], edgex[1], edgey[0], edgey[1])

    if region == 'str1' or region == 's1':
        
        ## convert to x and y coords
        if dr5:
            # also take into account the dodgy bit in bottom corner
            edgex, edgey = w.wcs_world2pix(strip1_ra, strip1_dec, 1)
            edgex_extra, edgey_extra = w.wcs_world2pix(strip1extra_ra, strip1extra_dec, 1)
            indicies = ((xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])) | ((xra < edgex_extra[1]) & (xra > edgex_extra[0]) & (ydec > edgey_extra[1]) & (ydec < edgey_extra[0]))
            
        else:
            edgex = [edges['s1min'], edges['s1max']]
            indicies = (xra < edgex[1]) & (xra > edgex[0])# & (ydec > edgey[1]) & (ydec < edgey[0])) | ((xra < edgex_extra[1]) & (xra > edgex_extra[0]) & (ydec > edgey_extra[1]) & (ydec < edgey_extra[0]))
          
    elif region == 'str2' or region == 's2':
        ## convert to x and y coords
        if dr5:
            edgex, edgey = w.wcs_world2pix(strip2_ra, strips_dec, 1)
            indicies = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
        else:
            edgex = [edges['s2min'], edges['s2max']]    
            indicies = ((xra < edgex[1]) & (xra > edgex[0])) # & (ydec > edgey[1]) & (ydec < edgey[0]))
        
    elif region == 'str3' or region == 's3':
        ## convert to x and y coords
        if dr5:
            edgex, edgey = w.wcs_world2pix(strip3_ra, strips_dec, 1)
            indicies = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
        else:
            edgex = [edges['s3min'], edges['s3max']]    
            indicies = ((xra < edgex[1]) & (xra > edgex[0])) # & (ydec > edgey[1]) & (ydec < edgey[0]))

    elif region == 'str4' or region == 's4':
        ## convert to x and y coords
        if dr5:
            edgex, edgey = w.wcs_world2pix(strip4_ra, strips_dec, 1)
            indicies = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
        else:
            edgex = [edges['s4min'], edges['s4max']] 
            indicies = ((xra < edgex[1]) & (xra > edgex[0])) # & (ydec > edgey[1]) & (ydec < edgey[0]))
    
    elif region == 'gap1' or region == 'g1':

        if dr5:
            edgex, edgey = w.wcs_world2pix(gap1_ra, strips_dec, 1)
            indicies = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
        else:
            ## do something different here!
            indicies = (xra < gap1_x[1]) & (xra > gap1_x[0]) & (ydec < gap1_y[1]) & (ydec > gap1_y[0])

    elif region == 'gap2' or region == 'g2':

        if dr5:
            edgex, edgey = w.wcs_world2pix(gap2_ra, strips_dec, 1)
            indicies = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
        else:
            ## do something different here!
            indicies = (xra < gap2_x[1]) & (xra > gap2_x[0]) & (ydec < gaps_y[1]) & (ydec > gaps_y[0])
            
    elif region == 'gap3' or region == 'g3':

        if dr5:
            edgex, edgey = w.wcs_world2pix(gap3_ra, strips_dec, 1)
            indicies = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
        else:
            ## do something different here!
            indicies = (xra < gap3_x[1]) & (xra > gap3_x[0]) & (ydec < gaps_y[1]) & (ydec > gaps_y[0])
        
    elif region == 'gap4' or region == 'g4':
        if dr5:
            edgex, edgey = w.wcs_world2pix(gap4_ra, strips_dec, 1)
            indicies = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
        else:
            ## do something different here!
            indicies = (xra < gap4_x[1]) & (xra > gap4_x[0]) & (ydec < gaps_y[1]) & (ydec > gaps_y[0])
    elif region == 'allstripsold':
        
        ## Extract the values for the full strips!
        ## s1
        edgex, edgey = w.wcs_world2pix(strip1_ra, strip1_dec, 1)
        edgex_extra, edgey_extra = w.wcs_world2pix(strip1_extra_ra, strip1_extra_dec, 1)
        stripone = ((xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])) | ((xra < edgex_extra[1]) & (xra > edgex_extra[0]) & (ydec > edgey_extra[1]) & (ydec < edgey_extra[0]))
        
        ## s2
        edgex, edgey = w.wcs_world2pix(strip2_ra, strips_dec, 1)
        striptwo = ((xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0]))
        
        ## s3
        edgex, edgey = w.wcs_world2pix(strip3_ra, strips_dec, 1)
        stripthree = ((xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0]))
        
        ## s4
        edgex, edgey = w.wcs_world2pix(strip4_ra, strips_dec, 1)
        stripfour = ((xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0]))
        
        print("In s1, s2, s3, s4 = ", np.sum(stripone), np.sum(striptwo), np.sum(stripthree), np.sum(stripfour))
        
        indicies = stripone | striptwo | stripthree | stripfour
        
    elif region == 'puregap':

        if dr5:
            # redo this completely...
            edgex, edgey = w.wcs_world2pix(gap1_ra, strips_dec, 1)
            gapone = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])

            edgex, edgey = w.wcs_world2pix(gap2_ra, strips_dec, 1)
            gaptwo = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])

            edgex, edgey = w.wcs_world2pix(gap3_ra, strips_dec, 1)
            gapthree = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])

            edgex, edgey = w.wcs_world2pix(gap4_ra, strips_dec, 1)
            gapfour = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])
            
        else:
            
            
            buff = 1000
            ## only the gaps
            edgex1, edgey = w.wcs_world2pix(strip1_ra, strips_dec, 1)
            gapone = (xra < (edgex1[0]-buff))
            #print "Num in g1 = ", np.sum(gapone), edgex1[0]+buff
            
            ## gap two
            edgex2, edgey = w.wcs_world2pix(strip2_ra, strips_dec, 1)
            gaptwo = (xra < (edgex2[0]-buff)) & (xra > (edgex1[1]+buff))
            #print "Num in g2 = ", np.sum(gaptwo), edgex2[0]+buff, edgex1[1]-buff
            
            ## gap three
            edgex3, edgey = w.wcs_world2pix(strip3_ra, strips_dec, 1)
            gapthree = (xra < (edgex3[0]-buff)) & (xra > (edgex2[1]+buff))
            #print "Num in g3 = ", np.sum(gapthree), edgex3[0]+buff,edgex2[1]-buff
            
            ## gap four
            edgex4, edgey = w.wcs_world2pix(strip4_ra, strips_dec, 1)
            gapfour = (xra < (edgex4[0]-buff)) & (xra > (edgex3[1]+buff))
            #print "Num in g4 = ", np.sum(gapfour), edgex4[0]+buff, edgex3[1]-buff
            
        indicies = gapone | gaptwo | gapthree | gapfour
        
    elif region == 'allstrips' or region == 'notstrips':

        if dr5:
            # redo again!
            edgex, edgey = w.wcs_world2pix(strip1_ra, strip1_dec, 1)
            edgex_extra, edgey_extra = w.wcs_world2pix(strip1extra_ra, strip1extra_dec, 1)
            stripone = ((xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])) | ((xra < edgex_extra[1]) & (xra > edgex_extra[0]) & (ydec > edgey_extra[1]) & (ydec < edgey_extra[0]))
            
            
            edgex, edgey = w.wcs_world2pix(strip2_ra, strips_dec, 1)
            striptwo = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])

            edgex, edgey = w.wcs_world2pix(strip3_ra, strips_dec, 1)
            stripthree = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])

            edgex, edgey = w.wcs_world2pix(strip4_ra, strips_dec, 1)
            stripfour = (xra < edgex[1]) & (xra > edgex[0]) & (ydec > edgey[1]) & (ydec < edgey[0])

        else:
            
            
            edgex = [edges['s1min'], edges['s1max']]
            stripone = (xra < edgex[1]) & (xra > edgex[0])
            
            edgex = [edges['s2min'], edges['s2max']]
            striptwo = (xra < edgex[1]) & (xra > edgex[0])
            
            edgex = [edges['s3min'], edges['s3max']]
            stripthree = (xra < edgex[1]) & (xra > edgex[0])
            
            edgex = [edges['s4min'], edges['s4max']]
            stripfour = (xra < edgex[1]) & (xra > edgex[0])
            
        indicies = stripone | striptwo | stripthree | stripfour
        
    else:
        print("Region ", region, " not recognised... ")
        exit()
    
    if region == 'notstrips':
        ## reverse the boolean array!
        returnindicies = np.invert(indicies)
    else:
        returnindicies = indicies
        
    print("There are ", np.sum(returnindicies), " in ", region)
    
    return returnindicies

def label_ct(inputCat, fieldName, plot = False, tbprov = False, catDir = '/mnt/vardy/vardygroupshare/data/masks/crosstalk/'):
    
    from astropy.coordinates import SkyCoord, match_coordinates_sky
    from astropy import units as u
  
    if fieldName == 'COSMOS':
        fieldName = 'COSMOS_BLAH'
    
    if tbprov:
        tb = inputCat
    else:
        tb = Table.read(inputCat)
        
    raSample = tb['RA']
    decSample = tb['DEC']

    # add a class column
    colnames = tb.colnames
    jj = (np.array(colnames) == 'CLASS')
    if np.any(jj) == False:
        # make column
        from astropy.table import Column
        newcol = Column(np.zeros(len(tb)), name = 'CLASS')
        tb.add_column(newcol)
    
    
    ############################
    ## CT setup
    sep = 0.34*128.0 # arcseconds
    sep = 43.27 # reduced, gives a better fit
    numRep = 10
    radius = 2.5 # 5.5 # arcseconds, to search for CT from centroid
    
    ii = fieldName.find('_')
    field = fieldName[:ii]
    
    ## Read in the positions of bright stars
    #csvfile = '../mask/UVISTA_DR1cat_brightJmag.csv'
    if field == 'COSMOS':
        fitsFile = catDir +field+'_J_bright_5as_new.fits' #starsFull.fits'
    else:
        fitsFile =  catDir +field+'_J_bright_5as.fits'
        
    print("Input file is ", fitsFile)
    
    #sampData = np.genfromtxt(csvfile, delimiter = ',', names = True, dtype = None)
    sampData = Table.read(fitsFile)
    
    ## cut at bright stars
    Jlim = 15.0 ##14.0
    mag = sampData['J']
    bb = (mag < Jlim)
    
    redData = sampData[bb]
    raStar = redData['RA']
    decStar = redData['DEC']
    
    ## Create an array of CT artefacts
    print("There are ", raStar.size, " bright stars.")
    ctRA = np.zeros(raStar.size*(numRep*2))
    ctDec = np.zeros(raStar.size*(numRep*2))
    for si in range(raStar.size):
        ## populate this array
        indmin = si*numRep*2
        indmax = indmin + numRep*2
        sepArray =  np.arange(-numRep, numRep+1, 1)
        sepArray = np.delete(sepArray, numRep)
        
        if field == 'COSMOS':
            ctRA[indmin:indmax] = raStar[si]
            ctDec[indmin:indmax] = decStar[si] + (sepArray*sep/3600.0) ## deg
        else:
            ctRA[indmin:indmax] = raStar[si] + (sepArray*sep/3600.0)
            ctDec[indmin:indmax] = decStar[si]  ## deg
            
    ## Just do an RA/DEC match with sky coords
    ## Make this a coordinate object
    c = SkyCoord(ra = raSample, dec = decSample)
    catalog = SkyCoord(ra = ctRA*u.degree, dec = ctDec*u.degree)    
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
      
    ## I want to find all matches
    ##idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, *u.deg) 
    
    dmsSeps = d2d.degree # .dms
    separation = dmsSeps*3600.0
        
    ## look in a separation of a few arcseconds
    ii = (separation < radius)
    
    ## classify as CT
    #newClass[ii] = 5
    
    #print tb['ID'][ii]
    
    # find the missing objects
    #print "Missing are ", tb['ID'][np.invert(ii)]
    
    print("Labelling CT = ", np.sum(ii))
    
    #t = Table([ctRA, ctDec], names = ['RA', 'DEC'])
    #t.write('CT_COSMOS.fits')
    
    if plot:
        plt.scatter(ctRA, ctDec)

        plt.scatter(raSample, decSample, color = 'blue')
        
        plt.show()
        
        exit()
        
    else:
        if np.sum(ii) > 0:
            tb['CLASS'][ii] = 5
        
    if tbprov:
        return tb
    else:
        
        tb.write(inputCat, overwrite = True)
        return np.sum(ii)

def aperture_correction(inputTable, apString, x, y, faster = True):
    ''' Code to get the aperture correction for each point.'''
    
    xgrid = inputTable['x']
    ygrid = inputTable['y']
    keys = inputTable.colnames
    #print("The input keys are ", keys)
    
    # = inputTable['depths']
    
    if 'ef_' + apString in keys:
        enclosedFlux = inputTable['ef_' + apString]
    else:
        print("APERTURE_CORECTION FUNCTION - Error, no column called ",'ef_' + apString)
        exit()
        
   ## Make an output array
    efArray = np.zeros(x.size)
    efArray[:] = -99.0
    
    if faster:
        
        #print("Using faster method.")
        #print("Input array size is ", x.size)
        deltay = np.min(ygrid)
        deltax = np.min(xgrid)
                
        ## loop through the grid instead of each object
        for xi in range(xgrid.size):
            
            xmin = xgrid[xi] - deltax
            xmax = xgrid[xi] + deltax
            ymin = ygrid[xi] - deltay
            ymax = ygrid[xi] + deltay
            
            ii = (x > xmin) & (x <= xmax) & (y > ymin) & (y <= ymax)
            efArray[ii] = 100.0*enclosedFlux[xi]
            
    else:
        # Find the closest point to the objects x and y positions
        for xi in range(x.size):
            
            # make a radius array
            deltax = (xgrid - x[xi])
            deltay = (ygrid - y[xi])
            radius = np.sqrt(deltax*deltax + deltay*deltay)
            mini = np.argmin(radius)
            
            #print "The closest point is ", xgrid[mini], ygrid[mini], " to ", x[xi], y[xi], enclosedFlux[mini]
            efArray[xi] = 100.0*enclosedFlux[mini]
        
    return efArray

def flux_conv(fluxCounts, fiveSigmaArray, zpt, microJy = True, cgs = False, minErrorPercent = 10.0, enclosedFlux = 100.0, resetNegative = False):
    ''' Code to convert counts into flux etc '''
    
    import numpy as np

    if resetNegative:
        conv = fluxCounts < 0
        if np.any(conv):
            print("Negative fluxes have been set to zero: FLUX_CONV") #, fluxCounts)
            fluxCounts[conv] = 0.0
            #print fluxCounts
    
    # find places where the five sigma array has not converged
    # and is set to -99

    # check for NANs
    if np.isnan(np.sum(fiveSigmaArray)):
        print('NANS found here in fiveSigma array')
        exit()
        
    rem = fiveSigmaArray < 1.0 # if any of the elements in 5sigArr<1, True

    fiveSigmaArray_allPositive = fiveSigmaArray
    if np.any(rem):
        fiveSigmaArray_allPositive[rem] = 1.0 # if any element is rem, assign 1
        print("reassigned depth array element to 1 ", fiveSigmaArray_allPositive[rem])
    
    if microJy:
        # simply convert flux counts into microJy
        print("zpt, fluxCounts: ", zpt, fluxCounts)
        fluxFinal = 3631000000.0*(10**(-zpt/2.5))*fluxCounts.squeeze()
        print("fluxFinal ", fluxFinal)
        # now convert the five sigma array into fluxes
        errorArray_allPositive = 726200000.0*(10**(-fiveSigmaArray/2.5))
        
    if cgs:
        # convert to cgs!
        value = -(48.6 + zpt)/2.5
        fluxFinal = (10**value)*fluxCounts.squeeze()
        
        # and the five sigma array
        value = -(48.6 + fiveSigmaArray)/2.5
        errorArray_allPositive = 0.2*(10**value)

    ##print fluxCounts
    # Check that the minimum error is in place
    if minErrorPercent > 0.0:
        print("minErrorPercent > 0.0 ",  minErrorPercent)
        print("errorArray_allPositive: ", errorArray_allPositive)
        tooSmall = errorArray_allPositive < (minErrorPercent/100.0)*fluxFinal.squeeze()
        if np.any(tooSmall):
            print("too small: True ",errorArray_allPositive[tooSmall])
            errorArray_allPositive[tooSmall] = (minErrorPercent/100.0) * fluxFinal.squeeze()[tooSmall]
            print("too small: True ",errorArray_allPositive[tooSmall])

    errorFinal = errorArray_allPositive
    if np.any(rem):
        print("rem=yes", errorFinal[rem])
        errorFinal[rem] = -99.0

    ## correct for the enclosed flux!
    fluxFinal = fluxFinal*100.0/enclosedFlux
    errorFinal = errorFinal*100.0/enclosedFlux
    print("enclosedFlux!!!",enclosedFlux)
    print(fluxFinal)
    print(errorFinal)
    
    return fluxFinal, errorFinal

def extract_match(catOne, catTwo, outputName = 'match.fits', overwrite = False, tidy = True, oneTable = False, names = ['1', '2']):

    if os.path.isfile(outputName) and (overwrite == False):
        print('Catalogue exists', outputName, '\nNOT matching again, set overwrite = True in extract_match to change this.')
        return
        
    # Read cat 1

    if oneTable:
        one = catOne
    else:
        one = Table.read(catOne)
    two = Table.read(catTwo)
    
    print('Matching {0} objects in {1}'.format(len(one), catOne))
    print('Matching {0} objects in {1}'.format(len(two), catTwo))
    
    # get the ra and dec
    colsone = np.array(one.colnames)
    colstwo = np.array(two.colnames)
    if np.any(colsone == 'RA'):
        raone = one['RA']
        decone = one['DEC']
    elif np.any(colsone == 'RAJ2000'):
        raone = one['RAJ2000']
        decone = one['DEJ2000']
    elif np.any(colsone == 'ALPHA_J2000'):
        raone = one['ALPHA_J2000']
        decone = one['DELTA_J2000']
    elif np.any(colsone == 'ra'):
        raone = one['ra']
        decone = one['dec']
    
    if np.any(colstwo == 'RA'):
        ratwo = two['RA']
        dectwo = two['DEC']
    elif np.any(colstwo == 'RAJ2000'):
        ratwo = two['RAJ2000']
        dectwo = two['DEJ2000']
    elif np.any(colstwo == 'ALPHA_J2000'):
        ratwo = two['ALPHA_J2000']
        dectwo = two['DELTA_J2000']
    elif np.any(colstwo == 'ra'):
        ratwo = two['ra']
        dectwo = two['dec']
        
    if raone.unit == None:
        from astropy import units as u
        raone = raone*u.deg
        decone = decone*u.deg
        
    if ratwo.unit == None:
        from astropy import units as u
        ratwo = ratwo*u.deg
        dectwo = dectwo*u.deg

    nant = np.any(np.isnan(ratwo))
    if nant:
        nanii = np.isnan(ratwo) | np.isnan(dectwo)
        nani = np.invert(nanii)
        ratwo = ratwo[nani]
        dectwo = dectwo[nani]
        two = two[nani]
        print('Shortened the second array to {0}'.format(len(two)))
        
    
    keep, idx = match_cats(raone, decone, ratwo, dectwo, return_all = True, verbose = False)
    
    
    # now join the catalogues together
    goodone = one[idx[keep]]
    goodtwo = two[keep]

    print('The number matched is ', len(goodone))

    if names[0] != '1':

        for i, nn in enumerate(goodone.colnames):
            goodone.rename_column(nn, nn + '_' + names[0])

        for i, nn in enumerate(goodtwo.colnames):
            goodtwo.rename_column(nn, nn + '_' + names[1])

            
    # now join these together!
    from astropy.table import hstack
    fulltable = hstack([goodone, goodtwo])

    if tidy:
        print('Remove junk columns here!')
        
    # save the table
    fulltable.write(outputName, overwrite = overwrite)
    print('Joined table of {0} objects saved to: {1}'.format(len(fulltable), outputName))
    
    return

def match_cats(xin, yin, xnew, ynew, separation = 1.0, return_nonmatches = False, return_all = False, verbose = True):
    
    ''' Separation is in arcseconds... '''
    
    import numpy as np
    from astropy.coordinates import SkyCoord, match_coordinates_sky
    from astropy import units as u
    
    # main catalogue
    catalog = SkyCoord(ra = xin, dec = yin)
    
    # new catalogue
    newc = SkyCoord(ra = xnew, dec = ynew)
    idx, sep2d, d3d = match_coordinates_sky(newc, catalog)
    #idx, sep2d, d3d = catalog.match_to_catalog_3d(newc)
    arcsec = sep2d*3600.0/u.degree

    if verbose:
        print("There are ", xnew.size, " objects to match. ")
    
    if return_nonmatches:
        keep = (arcsec > separation)
        if verbose:
            print(np.sum(keep), " objects have no match")
        if return_all:
            return keep, idx
        else:
            return keep
    else:
        matches = (arcsec < separation)
        if verbose:
            print(np.sum(matches), " objects have a match.")
        if return_all:
            return matches, idx
        else:
            return matches         



def do_phot(inputCat, imageDir, reqFilters, join = True, imageNames = ['None'], zeropoint = 21.58, apDiametersAS = np.array([2.8]), overwrite = False, outputName = 'test.fits'):
    print("Running do_phot: ", '\n')
    import sep
    from astropy.io import fits
    
    # get the ra and dec coordinates
    cat = Table.read(inputCat)
    cols = np.array(cat.colnames)
    if np.any(cols == 'RA'):
        ra = cat['RA']
        dec = cat['DEC']
    elif np.any(cols == 'RAJ2000'):
        ra = cat['RAJ2000']
        dec = cat['DEJ2000']

    if ra.unit == None:
        from astropy import units as u
        ra = ra*u.deg
        dec = dec*u.deg

    # get the required images
    if imageNames[0] == 'None':

        imageNames = np.chararray(len(reqFilters), itemsize = 200)
        
        # read in images.lis file.
        imageFile = imageDir + 'images.lis'
        if os.path.isfile(imageFile):
            imagedata = Table.read(imageFile, format = 'ascii.commented_header')

            # find filter
            for fi, filt in enumerate(reqFilters):

                ff = (imagedata['Name'] == filt)
                if np.any(ff):
                    imageNames[fi] = imagedata['Image'][ff][0]

    # first get the header info
    imageHeader = fits.getheader(imageDir + imageNames[0].decode('utf-8'))
    cdone_o = -3600.0*imageHeader['CD1_1']
    pixScale = round(cdone_o, 5)
    rpix = apDiametersAS/(pixScale*2.0)

    # to convert ra and dec into x/y
    from astropy import wcs
    w = wcs.WCS(imageHeader)

    xArraySmall, yArraySmall = w.all_world2pix(ra, dec, 1)
    
    rpixHuge = np.broadcast_to(rpix, (xArraySmall.size, rpix.size))
    xArray = np.repeat(xArraySmall[:,np.newaxis], rpix.size, 1)
    yArray = np.repeat(yArraySmall[:,np.newaxis], rpix.size, 1)
    
    for fi, filt in enumerate(reqFilters):

        print('Performing photometry on image:', imageNames[fi].decode('utf-8'))
        hdulist = fits.open(imageDir + imageNames[fi].decode('utf-8'), memmap = True)
        data = hdulist[0].data
        hdulist.close()

        data = data.byteswap().newbyteorder()
        flux, fluxerr, flag = sep.sum_circle(data, xArray, yArray, rpixHuge, subpix = 5)

        badpix = np.where(flux < 0.0)
        flux[badpix] = 1.0
        mag = -2.5*np.log10(flux) + zeropoint
        
        mag[badpix] = 99.0
        flux[badpix] = 0.0
        
        ## Create columns for this new data
        fluxAper = Column(flux, name = 'sep_f_{0}'.format(filt))
        magAper = Column(mag, name = 'sep_m_{0}'.format(filt))
        
        cat.add_column(fluxAper)
        cat.add_column(magAper)
        print('Photometry for filter {0} added to table.'.format(filt))

    # Now save this new catalogue
    cat.write(outputName, overwrite = overwrite)
    print('Updated catalogue saved to ', outputName)
    print("do_phot line 3140 finished", '\n')
    return

def compare_plot_phot(inputCat, filterA = 'NONE', filterB = 'NONE', radiusCat = 'none', outputDir = '../checks/plots_phot/', psfCorr = False, autoCat = 'none'):
    print("Running compare_plot_phot: ", '\n')
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.table import join
    
    # read in the cat
    data = Table.read(inputCat)
    print(data.colnames)
    print(inputCat)
    print(radiusCat)
    
    if radiusCat != 'none':
        radius = Table.read(radiusCat)
       # print(radius.colnames)
      #  print(len(data), len(radius))
        
        fulltb = join(data, radius, keys = ['ID', 'X_IMAGE', 'Y_IMAGE', 'RA', 'DEC'], table_names = ['', '_RAD'], uniq_col_name = '{col_name}{table_name}')
     #   print(fulltb.colnames)


    if autoCat != 'none':
        auto = Table.read(autoCat)
    #    print(auto.colnames)
        
        fulltb = join(fulltb, auto, keys = ['ID', 'X_IMAGE', 'Y_IMAGE', 'RA', 'DEC'], table_names = ['_RAD', '_AUTO'], uniq_col_name = '{col_name}{table_name}')


    print('Length of full table is ', len(fulltb))
        
    # for each filter make a plot
    #    filterA = 'HSC-G'
    #    filterB = 'CFHT-g'
    if psfCorr > -99.0:
        # read in the ecnlosed flux values
        print('PSF corrrect')

    if filterB == 'AUTO':
        magB = fulltb[filterA+'_AUTO']

        neg = (fulltb['flux_' + filterA] < 0.0)
        fluxforlog = fulltb['flux_' + filterA]
        print(fluxforlog[0:10])
        
        fluxforlog[neg] = 1E-60
        
        magA = -2.5*np.log10(fluxforlog) -48.6
        magA[neg] = -99.0
        
    else:
        magB = fulltb[filterB]
        magA = fulltb[filterA]
    
    radA = fulltb[filterA+'_RAD']
    diff = magA - magB

    maxpix = 5.0
    minmag = 27.0
    maxmag = 13.0
    
    # clean the table
    good = (magA < 90.0) & (magB < 90.0) & (magA > 0.0) & (magB > 0.0) & (radA < maxpix)

    fig = plt.figure()
    ax = plt.subplot(221)
    ax.scatter(magA[good], diff[good], s= 0.5, color = 'k')
    ax.plot(ax.get_xlim(), [0.0, 0.0], color = 'grey')
    print(ax.get_xlim())
    ax.set_ylabel('{0} - {1}'.format(filterA, filterB))
    ax.set_ylabel('{0}'.format(filterA))

    ax = plt.subplot(222)
    ax.scatter(magA[good], diff[good], s= 0.5, color = 'k')
    ax.plot(ax.get_xlim(), [0.0, 0.0], color = 'grey')
    ax.set_ylim([-2.0, 2.0])
    ax.set_xlim([maxmag, minmag])
    ax.set_ylabel('{0} - {1}'.format(filterA, filterB))
    ax.set_xlabel('{0}'.format(filterA))

    ax = plt.subplot(223)
    ax.scatter(magA[good], diff[good], s= 0.5, color = 'k')
    ax.plot(ax.get_xlim(), [0.0, 0.0], color = 'grey')
    ax.set_ylim([-0.3, 0.3])
    ax.set_xlim([maxmag, minmag])
    ax.set_ylabel('{0} - {1}'.format(filterA, filterB))
    ax.set_xlabel('{0}'.format(filterA))

    verygood = good & (magA < 23.0)
    medoffset = np.median(diff[verygood])
    print('The median offset is {0:.2f}'.format(medoffset))
    ax.plot(ax.get_xlim(), [medoffset, medoffset], color = 'red', linestyle = ':')
    
    ax = plt.subplot(224)
    ax.scatter(magA, radA, s = 0.5, color = 'k')
    ax.plot(ax.get_xlim(), [maxpix, maxpix], linestyle = ':', color = 'red')
#    ax.scatter(magA[good], radA[good])
    #    ax.set_xlim([15, 25])
    ax.set_ylim([0.0, 15.0])
    ax.set_ylabel('flux radius')
    plotname = outputDir + '{0}_{1}.png'.format(filterA, filterB)
    fig.savefig(plotname)
    print(plotname)
    return
