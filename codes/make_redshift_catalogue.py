# make redshift master catalogues - 18/09/25
# For z = 8
    # Locate masked mag catalogues JH, J, H                     - DONE
    # find ra/dec of dets that are in J but not JH              - DONE
    # vertically stack the detections in J that are not in JH   - DONE
    # find ra/dec of dets that are in H but not JH              - DONE
    # find ra/dec of dets that are in H but not J               - DONE
    # vstack detections that are in H but not in JH or in J     - DONE
    # ensure identical detections have consistent IDs           - DONE

# This produces c master cat only 667849 detections long

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

def new_ID_col(path, df):
    """ Inserts a new column which acts to preserve the original ID"""
    table = Table.read(path)

    # add new col so original ID can be traced
    new_col_name = f'New_ID'
    if new_col_name not in table.colnames:
        if df == 'JH':
            new_IDs = table[f'ID'].copy()
        elif df == 'J':
            new_IDs = table[f'ID'].copy() + int(3e6)
        elif df == 'H':
            new_IDs = table[f'ID'].copy() + int(6e6)

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
        #new_ID_col(path)

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

        # Remove rows with neighbours from Udf table and add the remaining rows to base table
        table.remove_rows(idx_U)
        distinctTableName = f'{Udf}_distinct_from_{df}.fits'
        distinctFromPath = os.path.join(globalCatDir, distinctTableName)
        distinctFromPaths.append(distinctFromPath) # len of list determines programme flow
        table.write(distinctFromPath, overwrite=True)
        print("Distinct detections table written to: ", distinctFromPath)

        ### Check if any still have neighbours within 1 as
        
        distinctFromCoords = prep_coords([distinctFromPath])

        idx_df, idx_U, sep2d, _ = distinctFromCoords[Udf]["coords"].search_around_sky(comparisonCoords[df]["coords"], maxsep_as)

        neighbour_counts = np.bincount(idx_U, minlength=len(distinctFromCoords)) 
        # uniqueOnlyCoords[i] has neighbour_counts[i] neighbours within 1"

        # Find Udf sources with more than zero neighbours: these are not unique detections
        has_neighbours = np.where(neighbour_counts > 0)[0]
        max_neighbours = neighbour_counts.max()
        #print("Max neighbours:", max_neighbours)
        #print("Sources with >0 neighbours:", len(has_neighbours))

    ### Do distinct from df1 and distinct from df2

    if len(distinctFromPaths)>1:
        ### get unique H

        # (H not in J) - in distinctFromPaths
        # (H not in JH) - in distinctFromPaths
        # ((H not in J) concat (H not in JH))
        # (H not in ((H not in J) concat (H not in JH))) = unique H

        # stack distinct_from_JH and distinct_from_J -> ((H not in J) concat (H not in JH))
        distinctFrom2Path = vstack_uniques(globalCatDir, distinctFromPaths[0], distinctFromPaths[1])


#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO#TODO
#TODO: unique_H_path is currently this: '/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/z8_master.fits' and has a length of only len(Table.read(unique_H_path))=232000
# it should be more like 600,000
#TODO: sort out the path for unique_H detections and join this with JH and unique J properly



        if 'master' in distinctFrom2Path:
            return distinctFrom2Path
        else:
            distinctFromFilts = os.path.basename(distinctFrom2Path).split('_distinct_from_')[1].replace('.fits','').split('_')[0:3:1]

            distinctFromFilts = '_'.join(distinctFromFilts)

            comparisonCoords = prep_coords([distinctFrom2Path])

            # (H not in ((H not in J) concat (H not in JH))) = unique H
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

            '''### TESTING/DEBUGGING ###

            topcat_H_not_JH = 64912
            topcat_H_not_J = 181166
            table_JH = Table.read(os.path.join(globalCatDir, 'H_distinct_from_JH.fits'))
            table_J = Table.read(os.path.join(globalCatDir, 'H_distinct_from_J.fits'))

            numUnique_JH = len(table_JH)
            numUnique_J = len(table_J)
            print(f"Number H not JH: " , numUnique_JH)
            print("difference: ", abs(topcat_H_not_JH - numUnique_JH))
            print(abs((topcat_H_not_JH - numUnique_JH))/topcat_H_not_JH *100, '%')

            print(f"Number H not J: " , numUnique_J)
            print("difference: ", abs(topcat_H_not_J - numUnique_J))
            print(abs((topcat_H_not_J - numUnique_J))/topcat_H_not_J *100, '%')

            # ((H not in J) concat (H not in JH))
            topcat_concat = 246078
            table_vstack = Table.read(distinctFrom2Path)
            numUnique_vstack= len(table_vstack)

            print(f"Number vstack {os.path.basename(distinctFrom2Path)}: " , numUnique_vstack)
            print("difference: ", abs(topcat_concat - numUnique_vstack))
            print(abs((topcat_concat - numUnique_vstack))/topcat_concat *100, '%')


            ### TESTING/DEBUGGING ###'''

            return uniqueTableName

    return distinctFromPaths


def vstack_uniques(globalCatDir, baseCatName, uniqueTableName):
    """
    Vertically stacks the tables provided 
    
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
    targetLength = len(baseTable) + len(uniqueTable)

    ### check col names and add new ID col
    for i in baseTable.colnames:
        if i not in uniqueTable.colnames:
            print("NO MATCH!: ", i, f"is not a column name in {uniquePath}")

    tablesToStack.append(baseTable) 
    tablesToStack.append(uniqueTable)

    # stack tables
    newTable = vstack(tablesToStack)

    if len(newTable) == targetLength:
        ## Write new fits
        # get filter names for new file name
        if 'det_' in baseCatPath:
            baseFiltName = baseCatPath.split('det_')[1].split('/')[0]
            uniqueFiltName = os.path.basename(uniqueTableName).split('_')[0]

            newPath = os.path.join(globalCatDir, f'{uniqueFiltName}_distinct_from_{baseFiltName}_vstacked.fits')
            newTable.write(newPath, overwrite=True)
        
            print("Written vertically stacked file to ", newPath)

        elif 'vstacked' in baseCatPath and '_distinct_from_' in baseCatPath:
            # J_distinct_from_JH_vstacked.fits 
            # H_distinct_from_JH_and_J.fits

            newPath = os.path.join(globalCatDir, 'z8_master.fits')
            newTable.write(newPath, overwrite=True)

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


J_path = new_ID_col(J_path, 'J')
JH_path = new_ID_col(JH_path, 'JH')
H_path = new_ID_col(H_path, 'H')

unique_J_table = find_unique(J_path, [JH_path])[0]

# Next, vertically stakc J-uniques to JH cat -> J-JH cat
JH_and_unique_J = vstack_uniques(globalCatDir, baseCatName, unique_J_table)

# Find detections distinct from JH and J: unique to H
# Vertically stakc H-uniques to J-JH cat 

unique_H = find_unique(H_path, [JH_path, J_path])

# stack JH and unique J, stack this with unique H

breakpoint()
master_z8_path = vstack_uniques(globalCatDir, JH_and_unique_J, unique_H)
print(master_z8_path)

### TESTING/DEBUGGING ###

topcat_z8 = 1637293
unique = len(Table.read(master_z8_path))

print(f"total detections in master: " , unique)
print("difference: ", abs(topcat_z8 - unique))
print(abs((topcat_z8 - unique))/topcat_z8 *100, '%')

### TESTING/DEBUGGING ###

























### SAVEPOINTS  ##########################################################################################
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



