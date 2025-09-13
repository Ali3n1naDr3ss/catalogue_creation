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
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.table import Table, hstack, join
from astropy.coordinates import SkyCoord, search_around_sky
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

def make_small_test_cat(filters, size=50):
    import os
    import numpy as np  
    import astropy.units as u
    from astropy.table import Table
    from astropy.coordinates import SkyCoord

    baseDir = '/raid/scratch/hullyott/cataloguing/current/'
    globalCatDir = os.path.join(baseDir,  'data/catalogues/COSMOS_make_master_test/')

    combinedCatPath = os.path.join(globalCatDir, "combined_cat.fits")
    print("Reading large table to make short cat...")
    table = Table.read(combinedCatPath) # the huge table with cols from every detFilt

    coordsDict = {} # pre-load all coords b4 matching
    
    for filt in filters:

        ra = table[f'RA_{filt}']
        ra = np.ma.MaskedArray(ra).filled(np.nan)
        ra = np.array(ra, dtype=float)
        idx = len(ra) // 2

        mid_ra = table[idx][f'RA_{filt}']
        mid_dec = table[idx][f'DEC_{filt}']

        start = int(idx - size/2)
        stop  = int(idx + size/2)

        shortCat = table[start:stop]
        shortCat.write(combinedCatPath.replace(".fits", "_short.fits"), overwrite=True)

# match by RA/Dec with maximum 1as seperation 
def match_objs(testing=True, verbose=False):

    import astropy.units as u
    from scipy.spatial import cKDTree # works in Cartesian, https://youtu.be/TLxWtXEbtFE
    from astropy.coordinates import SkyCoord, search_around_sky, match_coordinates_sky # couldn't get s_a_s nor m_c_s to work # TODO: look at source code to see if/how Astropy make their search robust against missed nearest neighbours
    
    if testing :
        combinedCatPath = os.path.join(globalCatDir, "combined_cat_short.fits") #TODO: "combined_cat.fits" full version
    else:
        combinedCatPath = os.path.join(globalCatDir, "combined_cat.fits") #TODO: "combined_cat.fits" full version
    
    print("Reading table ", os.path.basename(combinedCatPath), "...")
    table = Table.read(combinedCatPath) # the huge table with cols from every detFilt

    # make coordinate-like objects to use with search_around_sky
    coordsDict = {} # pre-load all coords b4 matching for speed
    
    for filt in detectionFilters:
        # get coords for filt
        ra = table[f'RA_{filt}']
        dec = table[f'DEC_{filt}']

        # replace maksed values with nan
        ra = np.ma.MaskedArray(ra).filled(np.nan)
        dec = np.ma.MaskedArray(dec).filled(np.nan)


        # make into arrays of floats
        ra = np.array(ra, dtype=float)
        dec = np.array(dec, dtype=float)

        # make into a 1D coord array
        #coordsDict[filt] = SkyCoord(ra=ra * u.deg, dec=dec * u.deg) # for search_around_sky attempt
        coordsDict[filt] = np.array([ra, dec])                       # for manual attempt

    for i, df in enumerate(detectionFilters[:-1]):
        coords1 = coordsDict[df]
        if detectionFilters[i+1] in coordsDict:
            coords2 = coordsDict[detectionFilters[i+1]]

        if coords1 is None or coords2 is None:
            print(f"Skipping {df} vs {detectionFilters[i+1]} (no valid coords)")
            continue

        """
        # search_around_sky attempt

        seplimit = 1 * u.arcsec # maximum allowed seperation to be considered the same object
        idx1, idx2, sep2d, _ = coords1.search_around_sky(coords2, seplimit)
        """
 
        # manual attempt
        # Build KDTree for catalog 1
        print(coords1)
        tree1 = cKDTree(coords1)

        # Query nearest neighbor in cat1 for each point in cat2
        distances, indices = tree1.query(coords2, k=1e6) #TODO: when using full-size cat, need to increase k (no of NN to return)
        idx1, idx2 = indices # idx1 indeces of NN in cat1
        # distances to the nearest neighbors, indecies of the NNs in the cats
        # Missing neighbors are indicated with infinite distances. Hits are sorted by distance (nearest first).


        ### Do manual and visual checks of nearest neighbours
        print("Found nearest neighbours... ")
        #for j in range(len(idx1)):
            #print("coords: ", coords1[0][idx1[j]], coords2[0][idx2[j]], "seperation: ", abs(coords2[0][idx2[j]] - coords1[0][idx1[j]]))
        ## visual checks
        central_points1 =  coords1[0][idx1], coords1[1][idx1] # RA/Dec of NN from cat 1
        central_points2 =  coords2[0][idx2], coords2[1][idx2] # RA/Dec of NN from cat 2
        
        # draw circle of radius 1as for validation
        maxsep_dimless = 1/3600 # 1 arcsec
        maxsep = maxsep_dimless * u.arcsec # maximum allowed seperation to be considered the same object
        
        central_x1 = coords1[1][idx1] # dec
        central_y1 = coords1[0][idx1] # RA

        fig, ax = plt.subplots(figsize=(6, 6))

        # Set limits and equal aspect before adding circles
        ax.set_xlim(min(central_points1[1].min(), central_points2[1].min()) - maxsep_dimless,
                    max(central_points1[1].max(), central_points2[1].max()) + maxsep_dimless)
        ax.set_ylim(min(central_points1[0].min(), central_points2[0].min()) - maxsep_dimless,
                    max(central_points1[0].max(), central_points2[0].max()) + maxsep_dimless)
        ax.set_aspect('equal', adjustable='datalim')  # ensures 1:1 scaling

        for l in range(len(central_x1)):
            validation_circle = Circle((central_x1[l], central_y1[l]), maxsep_dimless, fill=False, edgecolor='magenta', linewidth=1)
            ax.add_patch(validation_circle)

        ax.scatter(central_points1[1], central_points1[0], marker='.', color='red')
        ax.scatter(central_points2[1], central_points2[0], marker='.', color='green')
        ax.set_xlabel("Dec (deg)")
        ax.set_ylabel("RA (deg)")
        plt.title(f"Nearest Neighbours in Band {df}")

        plt.show()
        breakpoint()
       

        '''
        # cKDTree - class provides an index into a set of k-dimensional points which can be used to quickly look-up the nearest neighbors of any point
        # https://youtu.be/TLxWtXEbtFE
        # is only an approximate technique - some nearest neightbours may be missed
        tree1 = cKDTree(xyz1) # ra/dec space of points split by the median along each axis into quadrants
        tree2 = cKDTree(xyz2) # th nearest neighbour (NN) in each quadrat is assigned as NN - not perfect, but fast

        # Find all pairs within a given angular separation
        maxsep = 1 * u.arcsec # maximum allowed seperation to be considered the same object
        #r = 2 * np.sin(0.5 * maxsep.to_value(u.rad)) # assuming the sky is a sphere (in radians)
        r = maxsep # assuming the sky is flat
        pairs = tree1.query_ball_tree(tree2, r)'''


  


#####################################################################################################################################
# call functions - move this to an execution script when sub-pipeline finished

#modify colnames to preserve detFilt_measFilt
#run_parallel(modify_colnames, detectionFilters, max_workers=len(detectionFilters), use_processes=False)
#join_cats()
#make_small_test_cat(detectionFilters, size=100)
match_objs()
#run_parallel(func, input, max_workers=len(detectionFilters), use_processes=False)




