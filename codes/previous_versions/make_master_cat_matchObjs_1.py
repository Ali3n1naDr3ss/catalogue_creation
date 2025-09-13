# make_master_cat.py, started 04/09/25

# modify col names so that each col is traceable to the detection catalogue i.e. flux_Y_Y is from data/catalogues/dY_mY.fits
# join all tables together so that a single table contains all detections 
# match by RA/Dec with maximum 1as seperation 
# give the matched objects a unique ID
# this constitutes the Master Catalogue

import os
import time
import numpy as np
from datetime import datetime
from astropy.table import Table, hstack, join
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

baseDir = '/raid/scratch/hullyott/cataloguing/current/'
globalCatDir = os.path.join(baseDir,  'data/catalogues/COSMOS_make_master_test/')

detectionFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']
#####################################################################################################################################
print('\n')

def run_parallel(func, items, max_workers=24, use_processes=False,):
    """
    Run `func(item)` in parallel over a list of items.

    INPUTS
    func(callable)                  The function to run, must accept a single argument.
    items(list)                     The list of items to pass to `func`.
    max_workers(int, optional)      Maximum number of workers. Defaults to number of items or cores.
    use_processes(bool, optional)   If True, use processes (good for CPU-heavy tasks).
                                    If False, use threads (good for Input/Output*-heavy tasks).
                                    * I/O Tasks: 
                                    - Reading/writing large FITS files from disk.
	                                - Copying files, renaming columns, saving catalogues.
	                                - Downloading/uploading data over the network.
	                                - Waiting for an external program (e.g. running SExtractor via a subprocess).
	                                - Anything that spends most of its time waiting, not computing.

    OUTPUTS
    results(dict)                   Mapping of item -> result (or exception if failed).
    """

    print(f"Beginning parallel process for: {func.__name__}",  datetime.now().strftime("%d-%m-%Y %H:%M:%S"))
    start_time = time.time()

    Executor = ProcessPoolExecutor if use_processes else ThreadPoolExecutor

    results = {}
    with Executor(max_workers=max_workers) as executor:
        futures = {executor.submit(func, item): item for item in items}

        for fut in as_completed(futures):
            item = futures[fut]
            try:
                results[item] = fut.result()
                print(datetime.now().strftime("%H:%M:%S"), f"Finished: {item}")
            except Exception as e:
                results[item] = e
                print(f">>>>>>>>>> ERROR with {item}: {e}")

    elapsed = time.time() - start_time
    print(f"Finished {func.__name__} in {elapsed:.2f} seconds.")

    return results

### modify col names so that each col is traceable to the detection catalogue i.e. flux_Y_Y is from data/catalogues/dY_mY.fits
def modify_colnames(detFilt):
    """
    Writes new catalogue tables with updated column names so that the detection filter is 
    preserved when the filter-specific catalogues are joined together to create the master 
    catalogue. e.g. flux_Y is a column name in all catalogues up to this point. By renaming
    to flux_Y_Y, flux_Y_J, etc, the detection filter is preserved as the first filter name.
    """
    import shutil

    catDir = os.path.join(globalCatDir, 'det_'+detFilt+'/')
    catPrefix = 'COSMOS_DR2_MASKVISTADET_'
    catSuffix = '_1.8as_IRAC2.8as_cgs.fits'
    catName = catPrefix + detFilt+ catSuffix
    catPath = os.path.join(catDir, catName)

    # make a compy of original so you don't have to re-run everything if u fuck up ^^
    copycatPath = os.path.join(catDir, catPrefix + detFilt+ '_1.8as_IRAC2.8as_cgs_duplicate.fits')
    shutil.copy2(catPath, copycatPath)
    catTable = Table.read(copycatPath)

    renames = {}

    for colname in catTable.colnames[:]: # [:] - make a copy of colnames to avoid inadvertant overwrite
        if not colname.endswith(tuple(detectionFilters)):
            # new col name: ID_{detFilt}
            renames[colname] = colname +'_'+ detFilt

        elif colname.endswith(tuple(detectionFilters)):
            # new col name: param{Detection_measurement}
            param, measFilt = colname.split('_', 1)
            renames[colname] = param + detFilt + '_' + measFilt

    # Apply all renamings in one go
    catTable.rename_columns(list(renames.keys()), list(renames.values()))
    catTable.write(copycatPath, overwrite=True)

def join_cats():
    """ Locates catalogues with the modifed colnames and joins them into a single, large table. Writes table to  globalCatDir. """

    ### find the duplicated cats ###
    copycatPaths = []

    for detFilt in detectionFilters:
        catDir = os.path.join(globalCatDir, 'det_'+ detFilt +'/')
        catPrefix = 'COSMOS_DR2_MASKVISTADET_'
        copycatPath = os.path.join(catDir, catPrefix + detFilt+ '_1.8as_IRAC2.8as_cgs_duplicate.fits')
        if os.path.isfile(copycatPath):
            copycatPaths.append(copycatPath)

    combinedCatPath = copycatPaths[0]
    combinedCatTable = Table.read(combinedCatPath)

    ### join them and write file ###
    print(datetime.now().strftime("%H:%M:%S"), "Joining cats...")
    start_time = time.time()

    allTables = [Table.read(catPath, memmap=True) for catPath in copycatPaths]
    combinedCatTable = hstack(allTables)

    elapsed = time.time() - start_time
    print("Combined cats: ", combinedCatPath)
    print(f"Finished in {elapsed:.2f} seconds.")

    print(datetime.now().strftime("%H:%M:%S"), "writing combined cat...")
    start_time = time.time()

    combinedCatPath = os.path.join(globalCatDir, "combined_cat.fits")
    combinedCatTable.write(combinedCatPath, overwrite=True)

    elapsed = time.time() - start_time
    print("Written combined cat: ", combinedCatPath)
    print(f"Finished in {elapsed:.2f} seconds.")

    ### manual check ###
    command = f'topcat {combinedCatPath} &'
    #os.system(command)

# match by RA/Dec with maximum 1as seperation 
def match_objs():

    import astropy.units as u
    from astropy.coordinates import SkyCoord, search_around_sky
    
    combinedCatPath = os.path.join(globalCatDir, "combined_cat.fits")
    print("Reading large table...")
    table = Table.read(combinedCatPath) # the huge table with cols from every detFilt
    seplimit = 1 * u.arcsec 
    
    # make coordinate-like objects to use with search_around_sky
    coordsDict = {} # pre-load all coords b4 matching
    
    for filt in detectionFilters:
        # get coords for filt
        ra = table[f'RA_{filt}']
        dec = table[f'DEC_{filt}']
        # replace maksed values with nan
        ra = np.ma.MaskedArray(ra).filled(0.0)
        dec = np.ma.MaskedArray(dec).filled(0.0)
        #make into arrays of floats
        ra = np.array(ra, dtype=float)
        dec = np.array(dec, dtype=float)

        # make into a 1D coord array
        coordsDict[filt] = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    for i, df in enumerate(detectionFilters[:-1]):
        coords1 = coordsDict[df]
        if detectionFilters[i+1] in coordsDict:
            coords2 = coordsDict[detectionFilters[i+1]]

        if coords1 is None or coords2 is None:
            print(f"Skipping {df} vs {detectionFilters[i+1]} (no valid coords)")
            continue

        print("Matching:", df, "vs", detectionFilters[i+1])

        print("lllll", df, detectionFilters[i+1], coords1, coords2)
        idx, sep2d, sep3d = match_coordinates_sky(coords1, coords2)


    # neither search_around_sky, nor match_coordinates_sky work
    # alternate approach
    from scipy.spatial import cKDTree
    import numpy as np

    # Convert to unit vectors
    x1 = np.cos(coords1.ra.rad) * np.cos(coords1.dec.rad)
    y1 = np.sin(coords1.ra.rad) * np.cos(coords1.dec.rad)
    z1 = np.sin(coords1.dec.rad)
    xyz1 = np.column_stack([x1, y1, z1])

    x2 = np.cos(coords2.ra.rad) * np.cos(coords2.dec.rad)
    y2 = np.sin(coords2.ra.rad) * np.cos(coords2.dec.rad)
    z2 = np.sin(coords2.dec.rad)
    xyz2 = np.column_stack([x2, y2, z2])

# cKDTree - class provides an index into a set of k-dimensional points which can be used to quickly look-up the nearest neighbors of any point
# https://youtu.be/TLxWtXEbtFE
# is only an approximate technique - some nearest neightbours may be missed
tree1 = cKDTree(xyz1)
tree2 = cKDTree(xyz2)

# Find all pairs within a given angular separation (in radians)
r = 2 * np.sin(0.5 * seplimit.to_value(u.rad))
pairs = tree1.query_ball_tree(tree2, r)


   
    #idx1, idx2, sep2d, _ = search_around_sky(coords1, coords2, seplimit)

    print(f"Matched {df} with {detectionFilters[i+1]}: {len(idx1)} matches")
    # idx1 -> indices in cat
    # idx2 -> indices in cat2
    # sep2d -> angular separation



#####################################################################################################################################
# call functions - move this to an execution script when sub-pipeline finished

#modify colnames to preserve detFilt_measFilt
#run_parallel(modify_colnames, detectionFilters, max_workers=len(detectionFilters), use_processes=False)
#join_cats()
match_objs()
#run_parallel(func, input, max_workers=len(detectionFilters), use_processes=False)




