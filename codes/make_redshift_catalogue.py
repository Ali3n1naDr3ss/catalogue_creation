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

### Use mainVenv environment ###
# astropy v5.2
# numpy v1.24


baseDir = '/raid/scratch/hullyott/cataloguing/current/'
globalCatDir = os.path.join(baseDir, 'data/catalogues/COSMOS_experimental/')

z8_dfs = ['JH', 'J', 'H'] # why not K for z=8 detections?
# ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y'] 

testing = True
verbose = True

def new_ID_col(path):
    """ Inserts a new column which acts to preserve the original ID"""
    table = Table.read(path)

    # add new col so original ID can be traced
    new_col_name = f'original_ID'
    if new_col_name not in table.colnames:
        table.add_column(table[f'ID'].copy(), name=new_col_name, index=0)
    outputPath = path.replace('.fits', '_ID_col.fits')
    table.write(outputPath, overwrite=True)
    print("Written new table to: ", outputPath)

    return outputPath

def prep_coords(paths):
    """
    Prepares coordinates from FITS files as astropy coordinate on=bjects, ready for further processing
    """

    raw_coords = {}
    prepared_coords = {} # returned: coords as SkyCoord objects
    df2 = None

    for path in paths:

        table = Table.read(path)

        ### get coords for filt
        # make into a 1D coord arrays
        ID = np.array(table['ID'], dtype=int)
        ra = np.array(table['RA'], dtype=float)
        dec = np.array(table['DEC'], dtype=float)

        # Handle maksed values
        mask = ~np.isnan(ra) & ~np.isnan(dec)
        ra, dec, ID = ra[mask], dec[mask], ID[mask] ## Keep only valid entries

        if 'det_' in path:
            df = path.split('det_')[1].split('/')[0]
        #elif '_distinct_from_' and '_and_' in path:
            #df = path.split('/')[-1].split('_and_')[0].split('_distinct_from_')[1]
            #df2 = path.split('/')[-1].split('_and_')[1][0]ß

        elif '_distinct_from_' in path:
            df = path.split('/')[-1].split('_')[0]

        else:
            raise RuntimeError(f"\n\n>>>>>>>>> Could not get filter name from path: {path} \n")

        raw_coords[df] = {"ID": ID, "ra": ra, "dec": dec}

        if df2 is not None:
            raw_coords[df2] = {"ID": ID, "ra": ra, "dec": dec}

    ### prepare coords as SkyCoord object
    for df in raw_coords:
        ID = raw_coords[df]["ID"]
        ra = raw_coords[df]["ra"]
        dec = raw_coords[df]["dec"]

        coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg) # ALL coords in MASKED cat

        prepared_coords[df] = {"ID": ID, "coords": coords}

    return prepared_coords

def find_unique(targetUniquesPath, comparisonCatPaths):
    """
    Finds the RA/DEC of detections that are unique to the first catalogue when compared to RA/DEC in the second catalogue.
    i.e. finds detections that are in targetUniques but not in comparisonCat

    targetUniquesPath(str)  Path to the catalogue you want to find unique detections in.

    comparisonCatPath(str)  Path to the catalogue you want to compare targetUniques to.
    """
    
    print(datetime.now().strftime("%H:%M:%S"), "Finding matched detections...")

    ### get coords for comparison

    ## define the unique detection filter: Udf

    Udf = targetUniquesPath.split('det_')[1].split('/')[0]
    
    comparisonFilters = []
    for path in comparisonCatPaths:
        df = path.split('det_')[1].split('/')[0]
        comparisonFilters.append(df)

    ## get coords in correct format
    targetCoords = prep_coords([targetUniquesPath])
    comparisonCoords = prep_coords(comparisonCatPaths)

    ### Prepare search radius as Astropy quantity
    maxsep_as  = 1.0 # arcsec
    maxsep_deg = maxsep_as / 3600

    maxsep_deg = maxsep_deg * u.deg  # astropy quantity
    maxsep_as = maxsep_as * u.arcsec # astropy quantity

    ### Find unique sources
    ## Find number of neighbours each source has within radius
    distinctFromPaths = []

    for df, coords in zip(comparisonFilters, comparisonCoords):
        idx_df, idx_U, sep2d, _ = targetCoords[Udf]["coords"].search_around_sky(comparisonCoords[df]["coords"], maxsep_as)
        # coords[df]["coords"][idx_df] matches with taregtCoords[Udf]["coords"][idx_U]

        ## Remove all detections returned by search_around_sky

        print(f"Number of first pass {Udf}-{df} matches: ", len(idx_U))
        table = Table.read(targetUniquesPath) # table with unique and non-unique detections

        # Remove matched rows with neighbours from old table and add the remaining rows to base table
        table.remove_rows(idx_U)
        distinctTableName = f'{Udf}_distinct_from_{df}.fits'
        distinctFromPath = os.path.join(globalCatDir, distinctTableName)
        distinctFromPaths.append(distinctFromPath)
        table.write(distinctFromPath, overwrite=True)
        print("Distinct detections table written to: ", distinctFromPath)

        ### Check if any still have neighbours within 1 as
        ## Repeat search around each det in unique-only table
        # get unique-only coords
        
        distinctFromCoords = prep_coords([distinctFromPath])

        idx_df, idx_U, sep2d, _ = distinctFromCoords[Udf]["coords"].search_around_sky(comparisonCoords[df]["coords"], maxsep_as)
        neighbour_counts = np.bincount(idx_U, minlength=len(distinctFromCoords)) 
        # uniqueOnlyCoords[i] has neighbour_counts[i] neighbours within 1"

        # Find Udf sources with more than zero neighbours: these are not unique detections
        has_neighbours = np.where(neighbour_counts > 0)[0]
        max_neighbours = neighbour_counts.max()
        print("Max neighbours:", max_neighbours)
        print("Sources with >0 neighbours:", len(has_neighbours))

    ### Do distinct from df1 and distinct from df2

    if len(distinctFromPaths)>1:
        # get H not in distinct from JH_J
        # join distinct from JH and distinct from J tables and do s_a_s

        # stack distinct_from_JH and distinct_from_J
        distinctFrom2Path = vstack_uniques(globalCatDir, distinctFromPaths[0], distinctFromPaths[1]) 
        distinctFromFilts = os.path.basename(distinctFrom2Path).split('_distinct_from_')[1].replace('.fits','').split('_')[0:3:1]
        distinctFromFilts = '_'.join(distinctFromFilts)

        comparisonCoords = prep_coords([distinctFrom2Path])

        # H not in (HnotinJH and HnotinJ)
        df = Udf # compare H's from original H-cat to H's in stacked cats
        idx_df, idx_U, sep2d, _ = targetCoords[Udf]["coords"].search_around_sky(comparisonCoords[df]["coords"], maxsep_as)

        ## Remove all detections returned by search_around_sky

        print(f"Number of {Udf}-{distinctFromFilts} matches: ", len(idx_U))
        table = Table.read(targetUniquesPath) # table with unique and non-unique detections

        # Remove matched rows with neighbours from old table and add the remaining rows to base table
        table.remove_rows(idx_U)
        uniqueTableName = f'{Udf}_distinct_from_{distinctFromFilts}.fits'
        uniqueOnlyPath = os.path.join(globalCatDir, uniqueTableName)
        table.write(uniqueOnlyPath, overwrite=True)
        print("Unique detections table written to: ", uniqueOnlyPath)

        ### TESTING/DEBUGGING ###

        topcat_H = 658829
        numUnique = len(table)
        print(f"Number of unique {Udf}-detections: " , numUnique)
        print("difference: ", abs(topcat_H - numUnique))
        print(abs((topcat_H - numUnique))/topcat_H *100, '%')

        ### TESTING/DEBUGGING ###

        return uniqueTableName


    ### TESTING/DEBUGGING ###

    '''### Remove all sources with more than zero neighbours and write this as the unique table 
    table = Table.read(targetUniquesPath) # table with unique and non-unique detections

    # Remove matched rows with neighbours from old table and add the remaining rows to base table
    table.remove_rows(has_neighbours)
    uniqueTableName = f'only_unique_{Udf}.fits'
    table.write(os.path.join(globalCatDir, uniqueTableName), overwrite=True)'''

    topcat_unique_J = 123540
    numUnique = len(table)
    print(f"Number of unique {Udf}-detections: " ,numUnique)
    #print("difference: ", abs(topcat_unique - numUnique))
    #print(abs((topcat_unique - numUnique))/topcat_unique *100, '%')

    ### TESTING/DEBUGGING ###

    return distinctFromPaths


def vstack_uniques(globalCatDir, baseCatName, uniqueTableName):
    """
    
    baseCatPath(str)
    uniqueDetsToStack
    """
    # locate files
    baseCatPath = os.path.join(globalCatDir, baseCatName)
    uniquePath = os.path.join(globalCatDir, uniqueTableName)

    # read tables
    tablesToStack = [] # [base, unique detections]

    baseTable = Table.read(baseCatPath)
    uniqueTable = Table.read(uniquePath)
    print(len(baseTable))
    print(len(uniqueTable))
    targetLength = len(baseTable) + len(uniqueTable)
    print("target length: ",targetLength)
    # check col names
    for i in baseTable.colnames:
        if i not in uniqueTable.colnames:
            print("NO MATCH!: ", i, f"is not a column name in {uniquePath}")

    tablesToStack.append(baseTable) 
    tablesToStack.append(uniqueTable)

    # stack tables
    newTable = vstack(tablesToStack)
    print("Final length: ", len(newTable))
    print("Target length achieved?: ", len(newTable) == targetLength)

    if len(newTable) == targetLength:
        ## Write new fits

        # get filter names for new file name
        if 'det_' in baseCatPath:
            baseFiltName = baseCatPath.split('det_')[1].split('/')[0]
            uniqueFiltName = os.path.basename(uniqueTableName).split('_')[0]

            newPath = os.path.join(globalCatDir, f'{uniqueFiltName}_distinct_from_{baseFiltName}_vstacked.fits')
            newTable.write(newPath, overwrite=True)
        
            print("Written vertically stacked file to ", newPath)

        elif '_distinct_from_' in baseCatPath:
            path1 = baseCatPath
            path2 = uniqueTableName

            uniqueFiltName = os.path.basename(path1).split('_distinct_from_')[0]

            filtName1 = os.path.basename(path1).split('_distinct_from_')[1].replace('.fits','')
            filtName2 = os.path.basename(path2).split('_distinct_from_')[1].replace('.fits','')

            newPath = os.path.join(globalCatDir, f'{uniqueFiltName}_distinct_from_{filtName1}_and_{filtName2}_vstacked.fits')
            newTable.write(newPath, overwrite=True)
        
            print("Written vertically stacked file to ", newPath)

    return newPath


### FUNCTION CALL ################################################################################################
baseCatName = 'det_JH/COSMOS_DR2_MASKVISTADET_JH_1.8as_IRAC2.8as.fits'
JH_path = os.path.join(globalCatDir, baseCatName)
J_path = os.path.join(globalCatDir, 'det_J/COSMOS_DR2_MASKVISTADET_J_1.8as_IRAC2.8as.fits')
H_path = os.path.join(globalCatDir, 'det_H/COSMOS_DR2_MASKVISTADET_H_1.8as_IRAC2.8as.fits')

J_path = new_ID_col(J_path)
JH_path = new_ID_col(JH_path)
H_path = new_ID_col(H_path)

uniqueJTable = find_unique(J_path, [JH_path])[0]

# continuing as if J-uniques are found correctly. 
# Next, vertically stakc J-uniques to JH cat -> J-JH cat
J_JH_path = vstack_uniques(globalCatDir, baseCatName, uniqueJTable)

# Find detections distinct from JH and J: unique to H
# Vertically stakc H-uniques to J-JH cat 
H_distinct_JH_and_J = find_unique(H_path, [JH_path, J_path])






### SAVEPOINTS  ################################################################################################
def find_unique_single_comparison_backup(targetUniquesPath, comparisonCatPath):
    """
    Finds the RA/DEC of detections that are unique to the first catalogue when compared to RA/DEC in the second catalogue.
    i.e. finds detections that are in targetUniques but not in comparisonCat

    targetUniquesPath(str)  Path to the catalogue you want to find unique detections in.

    comparisonCatPath(str)  Path to the catalogue you want to compare targetUniques to.
    """
    
    print(datetime.now().strftime("%H:%M:%S"), "Finding matched detections...")

    ### get coords for comparison

    ## define the unique detection filter: Udf

    Udf = targetUniquesPath.split('det_')[1].split('/')[0]
    df = comparisonCatPath.split('det_')[1].split('/')[0]

    ## get coords in correct format
    coords = prep_coords([targetUniquesPath, comparisonCatPath])

    ### Prepare search radius as Astropy quantity
    maxsep_as  = 1.0 # arcsec
    maxsep_deg = maxsep_as / 3600

    maxsep_deg = maxsep_deg * u.deg  # astropy quantity
    maxsep_as = maxsep_as * u.arcsec # astropy quantity

    ### Find unique sources

    ## Find number of neighbours each source has within radius
    idx_df, idx_U, sep2d, _ = coords[Udf]["coords"].search_around_sky(coords[df]["coords"], maxsep_as)
    # coords[df]["coords"][idx_df] matches with coords[Udf]["coords"][idx_Udf]

    # remove all of these: search_around_sky

    print(f"Number of first pass {Udf}-{df} matches: ", len(idx_U))
    table = Table.read(targetUniquesPath) # table with unique and non-unique detections

    # Remove matched rows with neighbours from old table and add the remaining rows to base table
    table.remove_rows(idx_U)
    uniqueTableName = f'{Udf}_distinct_from_{df}.fits'
    uniqueOnlyPath = os.path.join(globalCatDir, uniqueTableName)
    table.write(uniqueOnlyPath, overwrite=True)
    print("Unique detections table written to: ", uniqueOnlyPath)

    ### Check if any still have neighbours within 1 as
    ## Repeat search around each det in unique-only table
    # get unique-only coords
    
    uniqueOnlyCoords = prep_coords([uniqueOnlyPath])

    idx_df, idx_U, sep2d, _ = uniqueOnlyCoords[Udf]["coords"].search_around_sky(coords[df]["coords"], maxsep_as)
    neighbour_counts = np.bincount(idx_U, minlength=len(uniqueOnlyCoords)) 
    # uniqueOnlyCoords[i] has neighbour_counts[i] neighbours within 1"

    # Find Udf sources with more than zero neighbours: these are not unique detections
    has_neighbours = np.where(neighbour_counts > 0)[0]
    max_neighbours = neighbour_counts.max()
    print("Max neighbours:", max_neighbours)
    print("Sources with >0 neighbours:", len(has_neighbours))

    ### TESTING/DEBUGGING ###

    '''### Remove all sources with more than zero neighbours and write this as the unique table 
    table = Table.read(targetUniquesPath) # table with unique and non-unique detections

    # Remove matched rows with neighbours from old table and add the remaining rows to base table
    table.remove_rows(has_neighbours)
    uniqueTableName = f'only_unique_{Udf}.fits'
    table.write(os.path.join(globalCatDir, uniqueTableName), overwrite=True)'''

    topcat_unique = 123540
    numUnique = len(table)
    print(f"Number of unique {Udf}-detections: " ,numUnique)
    print("difference: ", abs(topcat_unique - numUnique))
    print(abs((topcat_unique - numUnique))/topcat_unique *100, '%')

    ### TESTING/DEBUGGING ###

    return uniqueTableName



def find_unique_testing(coords, Udf, df):
    """
    Finds the RA/Dec of detections that are unique to the first catalogue referenced in 'coords', 
    when compared between the first two catalogues in 'coords'. 
    i.e. finds detections that are in coords[df1] but not in coords[df2]

    coords{dict}    Coordinates of detections in detection filter, 'df', 
                    in a prepared format for use in Astropy's match_to_catalog. 
                    coords[df] = {"ID": ID, "coords": coords} where coords = SkyCoord[(ra,dec)]

    Udf(str)        The detection filter you want to find Unique detections in. 
                    i.e. for z = 8, detections unique to J- and H-bands (compared to JH) are needed 

    df(str)         The base detection filter from which you want to compare i.e. for z = 8, JH
    """
    
    print(datetime.now().strftime("%H:%M:%S"), "Finding matched detections...")

    ### Prepare search radius as Astropy quantity
    maxsep_as  = 1.0 # arcsec
    maxsep_deg = maxsep_as / 3600

    maxsep_deg = maxsep_deg * u.deg  # astropy quantity
    maxsep_as = maxsep_as * u.arcsec # astropy quantity

    ### Find unique coords[Udf] sources

    ## Find how many neighbours each source has within radius
    idx_df, idx_U, sep2d, _ = coords[Udf]["coords"].search_around_sky(coords[df]["coords"], maxsep_as)
    # coords[df]["coords"][idx_df] matches with coords[Udf]["coords"][idx_Udf]
    print("Number of first pass J-JH matches: ", len(idx_U))
    neighbour_counts = np.bincount(idx_U, minlength=len(coords[Udf]["coords"])) 
    # coords[Udf]["coords"][i] has neighbour_counts[i] neighbours within 1"

    # Find Udf sources with more than zero neighbours
    has_neighbours = np.where(neighbour_counts > 0)[0]
    Udf_with_neighbours = coords[Udf]["ID"][has_neighbours]

    max_neighbours = neighbour_counts.max()
    more_than_one = (neighbour_counts > 1).sum()
    print("max neighbours:", max_neighbours)
    print("sources with >1 neighbours:", more_than_one)
    print("IDs with neighbours: ", Udf_with_neighbours)

    ### all detections within 1as are matches. remove tham all
    ##remove all dets form coords[Udf] that have > 0 neighbours

    path = os.path.join(globalCatDir, 'det_J/COSMOS_DR2_MASKVISTADET_J_1.8as_IRAC2.8as.fits')
    J_table = Table.read(path)
    # Remove matched rows from J_table and add the remaining rows to JH_table
    J_table.remove_rows(has_neighbours)
    uniqueTableName = f'unique_{Udf}.fits'
    J_table.write(os.path.join(globalCatDir, uniqueTableName), overwrite=True)
    unique_J = len(J_table)

    ### Assume this is ok for now - debuggind/testing stuff below this ###

    # how many neighbours do each of the neighbours_removed detections have? expect zero
    new_path = os.path.join(globalCatDir, 'neighbours_removed.fits')
    J_table = Table.read(new_path)
    df = 'J'
    ### get coords for filt
    # make into a 1D coord arrays
    ID = np.array(J_table['ID'], dtype=int)
    ra = np.array(J_table['RA'], dtype=float)
    dec = np.array(J_table['DEC'], dtype=float)

    # Handle maksed values
    mask = ~np.isnan(ra) & ~np.isnan(dec)
    ra, dec, ID = ra[mask], dec[mask], ID[mask] ## Keep only valid entries
    raw_coords = {}
    prepared_coords = {}
    raw_coords[df] = {"ID": ID, "ra": ra, "dec": dec}

    ID = raw_coords[df]["ID"]
    ra = raw_coords[df]["ra"]
    dec = raw_coords[df]["dec"]

    coordsJ = SkyCoord(ra=ra * u.deg, dec=dec * u.deg) # ALL coords in MASKED cat

    prepared_coords['J'] = {"ID": ID, "coords": coordsJ}

    idx_df, idx_U, sep2d, _ = prepared_coords['J']["coords"].search_around_sky(coords['JH']["coords"], maxsep_as)
    print("Number of second pass J-JH matches: ", len(idx_U))
    neighbour_counts = np.bincount(idx_U, minlength=len(prepared_coords['J']["coords"])) 
    # coords[Udf]["coords"][i] has neighbour_counts[i] neighbours within 1"

    # Find Udf sources with more than zero neighbours
    has_neighbours = np.where(neighbour_counts > 0)[0]

    max_neighbours = neighbour_counts.max()
    more_than_one = (neighbour_counts > 1).sum()
    print("NEW max neighbours:", max_neighbours)
    print("NEW sources with >1 neighbours:", more_than_one)
    print("NEW indices with neighbours: ", has_neighbours)

    # find the distance to the nearest neighbour

    idx, angsep, _ = prepared_coords['J']["coords"].match_to_catalog_sky(coords['JH']["coords"], nthneighbor=1)
    angsep = angsep.to(u.arcsec) # astropy quantity
    sorted_angsep = np.sort(angsep)
    close_dist = 1.23059
    mask = sorted_angsep < close_dist * u.arcsec
    print(f"Angles < {close_dist}arcsec:", sorted_angsep[mask])
    print("number of close neighbours: ", len(sorted_angsep[mask]))

    topcat_unique = 123540

    print("Number of unique J-detections: " ,unique_J)
    print("difference: ", abs(topcat_unique - unique_J))
    print(abs((topcat_unique - unique_J))/topcat_unique *100, '%')


    ### Assume this is ok for now - debuggind/testing stuff above this ###

    return unique_J, uniqueTableName



