# make_master_cat.py, started 04/09/25

# modify col names so that each col is traceable to the detection catalogue i.e. flux_Y_Y is from data/catalogues/dY_mY.fits
# join all tables together so that a single table contains all detections 
# match by RA/Dec with maximum 1as seperation 
# give the matched objects a unique ID
# this constitutes the Master Catalogue

import os
import time
import numpy as np
from random import randrange
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree # works in Cartesian, https://youtu.be/TLxWtXEbtFE
from itertools import combinations
from matplotlib.patches import Circle
from astropy.table import Table, hstack, join
from matplotlib.collections import PatchCollection
from astropy.coordinates import SkyCoord, search_around_sky
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

baseDir = '/raid/scratch/hullyott/cataloguing/current/'
globalCatDir = os.path.join(baseDir,  'data/catalogues/COSMOS_make_master_test/')

detectionFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']

size = 5 # number of rows to use in test cats
maxsep_as = 1 
testing = False
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
    os.system(command)

def make_small_test_cat(filters, size=50):
    """
    Makes a new catalogue with randomly sampled rows.
    """    
    
    import os
    import numpy as np  
    import astropy.units as u
    from astropy.table import Table
    from astropy.coordinates import SkyCoord

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
    os.system(command)

def find_matched_objects(maxsep_as, testing=True, verbose=False):
    """ 
    Writes a new table where identical detections that have been identified in multiples filters maintain the same ID number 

    maxsep(int/float) Maximum allowed seperation across different band images (in arcseconds) to be considered the same object 
    """

    import astropy.units as u
    from scipy.spatial import cKDTree # works in Cartesian, https://youtu.be/TLxWtXEbtFE
    from astropy.coordinates import SkyCoord, search_around_sky, match_coordinates_sky #TODO: look at source code to see if/how Astropy make their search robust against missed nearest neighbours
    
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

    for filt in detectionFilters:

        # get coords for filt
        ra = np.array(table[f'RA_{filt}'], dtype=float)
        dec = np.array(table[f'DEC_{filt}'], dtype=float)

        ### Handle maksed values
        mask = ~np.isnan(ra) & ~np.isnan(dec)
        ra, dec = ra[mask], dec[mask]  ## Keep only valid entries

        # make into a 1D coord array
        allCoordsDict[filt] = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        idDict[filt] = np.where(mask)[0]  # indices into the combined table (not the same as allCoords due to mask handling)
        
        # add new col so old ID can be traced
        new_col_name = f'UniqueID_{filt}'
        if new_col_name not in table.colnames:
            index = table.colnames.index(f'ID_{filt}')
            table.add_column(table[f'ID_{filt}'].copy(), name=f'New_ID_{filt}', index=index) # Insert at pos of old ID, initialise with copy of old IDs

    print("Finding matched objects...", datetime.now().strftime("%H:%M:%S"))

    for df1, df2 in combinations(detectionFilters, 2):  # compare all unique pairs
        coords1 = allCoordsDict[df1]
        coords2 = allCoordsDict[df2]

        maxsep_deg = maxsep_as / 3600
        seplimit = maxsep_deg * u.deg # astropy quantity for search_around_sky
        idx1, idx2, sep2d, _ = search_around_sky(coords1, coords2, seplimit)
        # idx1 the indics of the matching coords in coords1 
        # coords1[idx1[0]] matches with coords2[idx2[0]]
        sep2d = sep2d.value * 3600 * u.arcsec

        orig_idx1 = idDict[df1][idx1] # idx in combined cat
        orig_idx2 = idDict[df2][idx2]

        # re-assign unique ID numbers to the matched dets
        orig_df1_ID = table[f'ID_{detectionFilters[0]}'][orig_idx1] # allows extra traceability
        newUniqueID  = newIDBase + orig_df1_ID.data
        
        table[f'New_ID_{df1}'][orig_idx1] = newUniqueID # write over the old ID with unique ID
        table[f'New_ID_{df2}'][orig_idx2] = newUniqueID # NB: matched dets will not nece be at the same row

        #breakpoint()
    # write changes
    if testing:
        matchedCatPath = combinedCatPath.replace("_short.fits", "_matched_short.fits",)
    else:
        matchedCatPath = combinedCatPath.replace(os.path.basename(combinedCatPath), "master.fits",)
    table.write(matchedCatPath, overwrite=True)

    # open
    command = f'topcat {matchedCatPath} &'
    os.system(command)


def nn_plotter_save(coordsDict, nnIDdict, maxsep, testing=True):
    """
    Plots a circle of radius maxsep around each of the points provided. 
    Provided points should be the nearest neighbours found by KDTree.

    coordsDict{dict}    All coordinates in the combined (or test) catalogue for each of the detection filters {filt: Ra, Dec}
    nnIDdict{dict}      Indices of the nearest neighbours found by KDTree {filt: [idx]}
    maxsep(int/float)   Maximum allowed seperation across different band images (in arcseconds) to be considered the same object 
    """
    maxsep = maxsep / 3600
    print("This bit takes a lil while..........")

    fig, ax = plt.subplots(figsize=(6, 6))

    # Set limits and aspect
    ax.set_xlim(min(coordsDict[df1][0]) - 0.1, max(coordsDict[df1][0]) + 0.1)
    ax.set_ylim(min(coordsDict[df1][1]) - 0.1, max(coordsDict[df1][1]) + 0.1)
    ax.set_aspect('equal', adjustable='datalim')

    # RA/Dec for filter 1
    df1_ra  = coordsDict[df1][0][nnIDdict[df1]]
    df1_dec = coordsDict[df1][1][nnIDdict[df1]]

    for NNra, NNdec in zip(df1_ra, df1_dec):
        for indiRA, indiDec in zip(NNra, NNdec):
            validation_circle = Circle((indiRA, indiDec),
                maxsep, fill=False, edgecolor='magenta', linewidth=1)
            ax.add_patch(validation_circle)

    # RA/Dec for filter 2
    df2_ra  = coordsDict[df2][0][nnIDdict[df2]]
    df2_dec = coordsDict[df2][1][nnIDdict[df2]]

    for NNra, NNdec in zip(df2_ra, df2_dec):
        for indiRA, indiDec in zip(NNra, NNdec):
            validation_circle = Circle((indiRA, indiDec),
                maxsep, fill=False, edgecolor='magenta', linewidth=1)
            ax.add_patch(validation_circle)

    for df in detectionFilters:
        # colours for plot
        cm = plt.colormaps.get_cmap("rainbow")
        a =  len(detectionFilters)
        i = detectionFilters.index(df)
        colour = cm(i/a)

        ## plot the scatter points
        df_ra = coordsDict[df][0][nnIDdict[df]]
        df_dec = coordsDict[df][1][nnIDdict[df]]
        ax.scatter(df_ra, df_dec, marker='.', color=colour, label=df)

    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    if testing:
        plt.title(f"Nearest Neighbours to each det for middle {size} rows")
        print(">>>>>> NOTE: for small test catalogues, detections will appear as strips since all detections are takens from rows in the main catalogue.")
    else:
        plt.title(f"Nearest Neighbours to each det. Maxsep: {maxsep}as")
    plt.legend()
    plt.show()


#####################################################################################################################################
# call functions - move this to an execution script when sub-pipeline finished

# modify colnames to preserve detFilt_measFilt
run_parallel(modify_colnames, detectionFilters)
join_cats()

if testing:
    make_small_test_cat(detectionFilters, size)

find_matched_objects(maxsep_as, testing, verbose)


