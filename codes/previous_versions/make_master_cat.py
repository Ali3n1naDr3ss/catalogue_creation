# make_master_cat.py, started 04/09/25

# modify col names so that each col is traceable to the detection catalogue i.e. flux_Y_Y is from data/catalogues/dY_mY.fits
# join all tables together so that a single table contains all detections 
# match by RA/Dec with maximum 1as seperation 
# give the matched objects a unique ID
# this constitutes the Master Catalogue

import os
import time
import h5py
import numpy as np
import astropy.units as u
from random import randrange
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree # works in Cartesian, https://youtu.be/TLxWtXEbtFE
from itertools import combinations
from matplotlib.patches import Circle
from astropy.table import Table, hstack, join
from matplotlib.collections import PatchCollection
from astropy.coordinates import SkyCoord, search_around_sky #TODO: look at source code to see if/how Astropy make their search robust against missed nearest neighbours
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

baseDir = '/raid/scratch/hullyott/cataloguing/current/'
globalCatDir = os.path.join(baseDir, 'data/catalogues/COSMOS_make_master_test/')
masterCatPath = os.path.join(globalCatDir, 'master.fits')
matchedDetsPath = os.path.join(globalCatDir, 'matchedDets.h5')

detectionFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y'] 

start_fresh = False
size = 5000 # number of rows to use in test cats
maxsep_as = 1 
testing = True
verbose = True
#####################################################################################################################################
print('\n')

def run_parallel(func, items, max_workers=len(detectionFilters), use_processes=False):
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

    # make a copy of original so you don't have to re-run everything if u fuck up ^^
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

def make_small_test_cat(filters, size=50):
    """
    Makes a new catalogue with randomly sampled rows.
    """    

    baseDir = '/raid/scratch/hullyott/cataloguing/current/'
    globalCatDir = os.path.join(baseDir, 'data/catalogues/COSMOS_make_master_test/')

    combinedCatPath = os.path.join(globalCatDir, "combined_cat.fits")

    print("Reading large table to make short cat...")
    table = Table.read(combinedCatPath)

    # make list of random indices
    random_indis = []
    for _ in range(size):
        random_idx = randrange(len(table)//2) #TODO: there's probs gunna be an issue when there are none detections in some filters - do masking first
        random_indis.append(random_idx)
    #print("random indices: ", random_indis)

    # build the short catalogue
    shortCat = table[random_indis]

    # write test cat
    shortCatPath = combinedCatPath.replace(".fits", "_short.fits")
    shortCat.write(shortCatPath, overwrite=True)

    # open
    command = f'topcat {shortCatPath} &'
    #os.system(command)

def find_matched_dets(maxsep_as, matchedDetsPath, testing=True, verbose=False):
    """ 
    Writes a new table where identical detections that have been identified in multiples filters maintain the same ID number 

    maxsep(int/float) Maximum allowed seperation across different band images (in arcseconds) to be considered the same object 
    """
    
    if testing:
        combinedCatPath = os.path.join(globalCatDir, "combined_cat_short.fits")
    else:
        combinedCatPath = os.path.join(globalCatDir, "combined_cat.fits")
    
    print("Reading table ", os.path.basename(combinedCatPath), "...")

    table = Table.read(combinedCatPath) # the huge table with cols from every detFilt

    # find largest ID
    idCols = []
    for col in table.colnames:
        if 'ID_' in col:
            print(col)
            idCols.append(table[col][-1:])
    largestID = max(idCols).item()
    newIDBase = int(round(largestID/1e6)*1e6)

    # make coordinate-like objects to use with search_around_sky
    allCoordsDict = {} # pre-load all coords b4 matching for speed
    idDict = {} # the row indices of duplicate detections in the 'combined' catalogue
    newIDs = [] # new IDs assigned to matched dets

    for filt in detectionFilters:

        # get coords for filt
        ra = np.array(table[f'RA_{filt}'], dtype=float)
        dec = np.array(table[f'DEC_{filt}'], dtype=float)

        ### Handle maksed values
        mask = ~np.isnan(ra) & ~np.isnan(dec)
        ra, dec = ra[mask], dec[mask]  ## Keep only valid entries

        # make into a 1D coord array
        allCoordsDict[filt] = SkyCoord(ra=ra * u.deg, dec=dec * u.deg) # all coords in MASKED combined cat
        idDict[filt] = np.where(mask)[0]  # indices into the combined table (not the same as allCoords due to mask handling)
        
        # add new col so old ID can be traced
        new_col_name = f'New_ID_{filt}'
        if new_col_name not in table.colnames:
            index = table.colnames.index(f'ID_{filt}')
            table.add_column(table[f'ID_{filt}'].copy(), name=f'New_ID_{filt}', index=index) # Insert at pos of old ID, initialise with copy of old IDs

    print(datetime.now().strftime("%H:%M:%S"), "Finding matched detections...")

    for df1, df2 in combinations(detectionFilters, 2):  # compare all unique pairs
        coords1 = allCoordsDict[df1]
        coords2 = allCoordsDict[df2]

        maxsep_deg = maxsep_as / 3600
        seplimit = maxsep_deg * u.deg # astropy quantity for search_around_sky
        # Y detections dont have matched pairs in figure - is it bc S_a_s allows Y-Y matching?
        idx1, idx2, sep2d, _ = search_around_sky(coords1, coords2, seplimit)
        # idx1 the indics of the matching coords in coords1 
        # coords1[idx1[0]] matches with coords2[idx2[0]]
        sep2d = sep2d.value * 3600 * u.arcsec

        matchedCoords_df1 = coords1[idx1]
        matchedCoords_df2 = coords2[idx2]

        orig_idx1 = idDict[df1][idx1] # idx in original combined cat
        orig_idx2 = idDict[df2][idx2]

        # re-assign unique ID numbers to the matched dets
        orig_df1_ID = table[f'ID_{detectionFilters[0]}'][orig_idx1] # allows extra traceability
        newID  = newIDBase + orig_df1_ID.data

        table[f'New_ID_{df1}'][orig_idx1] = newID # write over the old ID with unique ID
        table[f'New_ID_{df2}'][orig_idx2] = newID # NB: matched dets will not nece be at the same row
        newID = newID.data.tolist()
        newIDs.append(newID)

    # write matched detections info to file for later use
    print(datetime.now().strftime("%H:%M:%S"), "Writing matched detections file...")

    newIDs = sum(newIDs, [])
    newIDs = np.array(newIDs, dtype=np.int64) # convert to arr for h5 functionality and memory efficiency
    with h5py.File(matchedDetsPath, "w") as f:
        f.create_dataset('newIDs', data=newIDs)

    # write changes
    if testing:
        matchedCatPath = os.path.join(globalCatDir, "master_short.fits",)
    else:
        matchedCatPath = masterCatPath

    print(datetime.now().strftime("%H:%M:%S"), f"Writing {masterCatPath}")
    table.write(matchedCatPath, overwrite=True)

    # open
    command = f'topcat {matchedCatPath} &'
    os.system(command)


def matched_plotter(maxsep_as, matchedDetsPath, testing=True):
    """
    AIM: 
    plots the mathced objects and a validation circle of radius = 1as around each 
    (points and circles have colours that relate to their filters)

    maxsep_as(int/float)   Maximum allowed seperation across different band images (in arcseconds) to be considered the same object 
    
    CURRENT FUNCTIONALITY: 
    plots the mathced objects and Y OBJECTS  without a matching counterpart as well as validation circle of radius = 1as around each 
    """
    print(datetime.now().strftime("%H:%M:%S"), f"Plotting matched detections from {masterCatPath}")

    # assuming the IDs are correct - check them - get the ra dec of the objects at that id
    maxsep_deg = maxsep_as / 3600
    a = len(detectionFilters)

    fig, ax = plt.subplots(figsize=(6, 6))

    # Set limits and aspect
    ax.set_xlim(148.5, 151)
    ax.set_ylim(2, 2.3)
    ax.set_aspect('equal', adjustable='datalim')

    for df in detectionFilters:
        with h5py.File(matchedDetsPath, "r") as f:     
            newIDs = f["newIDs"][:] # [:] read entire file
        
        table = Table.read(masterCatPath) # the huge table with cols from every detFilt
            
        # Build mask: True if table ID is in newIDs
        mask = np.isin(table[f"New_ID_{df}"], newIDs)

        # Apply mask
        matchedRows = table[mask]
        excludedRows = table[~mask]

        ra = matchedRows[f'RA_{df}']
        dec = matchedRows[f'DEC_{df}']
        
        if len(ra) > 5000:
            print(datetime.now().strftime("%H:%M:%S"), f"Getting random sample for {df}...")
            random_indis = np.random.choice(len(matchedRows), size=500, replace=False)

            # build the short catalogue
            ra = ra[random_indis]
            dec = dec[random_indis]

        # colours for plot
        cm = plt.colormaps.get_cmap("rainbow")
        i = detectionFilters.index(df)
        colour = cm(i/a)

        ## plot the scatter points
        ax.scatter(ra, dec, marker='.', color=colour, label=df)

        ## plot validation circle
        for ra1, dec1 in zip(ra, dec):
            validation_circle = Circle((ra1, dec1),
                    maxsep_deg, fill=False, edgecolor=colour, linewidth=1)
            ax.add_patch(validation_circle)

    plt.title(f"Matched multi-band detections Maxsep: {maxsep_as}as")
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    plt.legend()
    plt.show()



#####################################################################################################################################
# call functions - move this to an execution script when sub-pipeline finished
if testing:
    masterCatPath = os.path.join(globalCatDir, 'master_short.fits')
    matchedDetsPath = matchedDetsPath.replace('.h5', '_short.h5') 

if start_fresh==True or os.path.isfile(masterCatPath)==False: 
    # modify colnames to preserve detFilt_measFilt
    run_parallel(modify_colnames, detectionFilters)
    join_cats()

    if testing:
        make_small_test_cat(detectionFilters, size)

    find_matched_dets(maxsep_as, matchedDetsPath, testing, verbose)
    matched_plotter(maxsep_as, matchedDetsPath)

elif start_fresh==False and os.path.isfile(masterCatPath)==True:
    find_matched_dets(maxsep_as, matchedDetsPath, testing, verbose)
    matched_plotter(maxsep_as, matchedDetsPath)





