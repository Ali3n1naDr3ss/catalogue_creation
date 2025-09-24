# make redshift master catalogues - 18/09/25
# For z = 8
    # Locate masked mag catalogues JH, J, H
    # find ra/dec of dets that are in J but not JH
    # vertically stack the detections in J that are not in JH
    # find ra/dec of dets that are in J but not JH
    # vstack detections that are in H but not in JH-J

    # ensure identical detections have consistent IDs

import os
import time
import numpy as np
import astropy.units as u
from datetime import datetime
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord, search_around_sky
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

baseDir = '/raid/scratch/hullyott/cataloguing/current/'
globalCatDir = os.path.join(baseDir, 'data/catalogues/COSMOS_experimental/')

z8_dfs = ['JH', 'J', 'H'] # why not K for z=8 detections?
# ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y'] 


testing = True
verbose = True

def find_dets_not_in_other():
    """
    Finds the RA/Dec of detections that are unique to J, when compared between J and JH. i.e. finds detections that are in J but not in JH
    """
    # locate files
    catPaths = []
    for df in z8_dfs:
        fDir = f'det_{df}/COSMOS_DR2_MASKVISTADET_{df}_1.8as_IRAC2.8as.fits'
        catpath = os.path.join(globalCatDir, fDir)
        catPaths.append(catpath)

    # get each table
    JHdir = 'det_JH/COSMOS_DR2_MASKVISTADET_JH_1.8as_IRAC2.8as.fits'
    JHpath = os.path.join(globalCatDir, JHdir)

    Jdir = 'det_J/COSMOS_DR2_MASKVISTADET_J_1.8as_IRAC2.8as.fits'
    Jpath = os.path.join(globalCatDir, Jdir)

    JHtable = Table.read(JHpath)
    Jtable = Table.read(Jpath)

    # get coords for filt
    JHra = np.array(JHtable['RA'], dtype=float)
    JHdec = np.array(JHtable['DEC'], dtype=float)

    Jra = np.array(Jtable['RA'], dtype=float)
    Jdec = np.array(Jtable['DEC'], dtype=float)

    # Handle maksed values
    mask = ~np.isnan(JHra) & ~np.isnan(JHdec) # ~: if not nan...
    JHra, JHdec = JHra[mask], JHdec[mask]  ## Keep only valid entries

    mask = ~np.isnan(Jra) & ~np.isnan(Jdec)
    Jra, Jdec = Jra[mask], Jdec[mask]  ## Keep only valid entries

    # make into a 1D coord arrays
    coordsJH = SkyCoord(ra=JHra * u.deg, dec=JHdec * u.deg) # all coords in MASKED cat
    coordsJ = SkyCoord(ra=Jra * u.deg, dec=Jdec * u.deg)    # all coords in MASKED cat

    print(datetime.now().strftime("%H:%M:%S"), "Finding matched detections...")

    maxsep_deg = 1 / 3600 # 1as in deg
    maxsep_deg = maxsep_deg * u.deg # astropy quantity for search_around_sky

    ## Find unique J sources
    idx, angsep, _ = coordsJ.match_to_catalog_sky(coordsJH) # nearest neighbours in J and JH

    unique_to_J_mask = angsep > maxsep_deg # exclude if NN is separated by more than 1as
    unique_to_J = coordsJ[unique_to_J_mask] # apply mask

    topcat_unique_to_J = 123540
    print(len(unique_to_J))
    print("diff: ", topcat_unique_to_J - len(unique_to_J))
    print((topcat_unique_to_J - len(unique_to_J))/topcat_unique_to_J *100)

find_dets_not_in_other()

"""
13:14:58 Finding matched detections...
116193
diff:  7347
5.947061680427392

"""




