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

def prep_coords():
    """
    Prepares coordinates from FITS files as astropy coordinate on=bjects, ready for further processing
    """

    raw_coords = {}
    prepared_coords = {} # returned: coords as SkyCoord objects

    for df in z8_dfs:
    ### locate files
        Dir = f'det_{df}/COSMOS_DR2_MASKVISTADET_{df}_1.8as_IRAC2.8as.fits'
        path = os.path.join(globalCatDir, Dir)
        
        table = Table.read(path)

        ### get coords for filt
        # make into a 1D coord arrays
        ID = np.array(table['ID'], dtype=int)
        ra = np.array(table['RA'], dtype=float)
        dec = np.array(table['DEC'], dtype=float)

        # Handle maksed values
        mask = ~np.isnan(ra) & ~np.isnan(dec)
        ra, dec, ID = ra[mask], dec[mask], ID[mask] ## Keep only valid entries

        raw_coords[df] = {"ID": ID, "ra": ra, "dec": dec}

    ### prepare coords as SkyCoord object
    for df in raw_coords:
        ID = raw_coords[df]["ID"]
        ra = raw_coords[df]["ra"]
        dec = raw_coords[df]["dec"]

        coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg) # ALL coords in MASKED cat

        prepared_coords[df] = {"ID": ID, "coords": coords}

    return prepared_coords

def find_unique(coords, Udf, df):
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
    maxsep_as  = 1 # arcsec
    maxsep_deg = maxsep_as / 3600

    maxsep_deg = maxsep_deg * u.deg  # astropy quantity
    maxsep_as = maxsep_as * u.arcsec # astropy quantity

    ### Find unique coords[Udf] sources

    ## Find how many neighbours each source has within radius
    idx_U, idx_df, sep2d, _ = coords[Udf]["coords"].search_around_sky(coords[df]["coords"], maxsep_as)
    # coords[df]["coords"][idx_U] matxhes with coords[Udf]["coords"][idx_df]

    neighbour_counts = np.bincount(idx_U, minlength=len(coords[Udf]["coords"])) 
    # number of mathces (neighbours closer than 1as) each detection in Udf has 
    # i.e. coords[Udf][0] has neighbour_counts[0] neighbours closer than 1as
    max_neighbours = neighbour_counts.max()
    more_than_one = (neighbour_counts > 1).sum()
    print("max neighbours:", max_neighbours)
    print("sources with >1 neighbours:", more_than_one)


    unique_IDs_dict = {}
    unique_coords_dict = {}
    for i in range(max_neighbours):
        i = i+1

        ## check if ith NN is within 1 as
        idx, angsep, _ = coords[Udf]["coords"].match_to_catalog_sky(coords[df]["coords"], nthneighbor=i)
        # nearest neighbours in Udf and df
        # https://docs.astropy.org/en/latest/_modules/astropy/coordinates/sky_coordinate.html#SkyCoord.match_to_catalog_sky
        #print(coords["coords"][Udf], coords[df]["coords"][idx]) # Each source in coords[Udf], is ra/dec matched to coords[df][idx]

        angsep = angsep.to(u.arcsec) # astropy quantity

        mask = angsep > maxsep_as # include as unique if NN is separated by more than 1as

        unique_coords = coords[Udf]["coords"][mask] # apply mask resulting in detections unique to Udf
        unique_IDs = coords[Udf]["ID"][mask]        # ID numbers of unique detections as listed in MASKVISTADET table

        unique_IDs_dict[f"{i}th_neighbour"] = unique_IDs
        unique_coords_dict[f"{i}th_neighbour"] = unique_coords
        print(f"Number of {i}th order unique detections in {Udf}: {len(unique_coords)}")
        print("Removed this round: ", len(coords[Udf]["ID"]) - len(unique_IDs))


    ### Write over the Udf entry in coords dict, retaining the IDs from the MASKVISTADET cat, 
    ### and return coords dict so find_unique can be repeated
    coords[Udf]["coords"] = unique_coords
    coords[Udf]["ID"] = unique_IDs

    return coords

### FUNCTION CALL ################################################################################################

coords = prep_coords()
unique_J_dict = find_unique(coords, 'J', 'JH')

topcat_unique = 123540
unique_J = unique_J_dict['J']["coords"]
print("Number of unique J-detections: " ,len(unique_J))
print("difference: ", abs(topcat_unique - len(unique_J)))
print(abs((topcat_unique - len(unique_J)))/topcat_unique *100, '%')

# re-run find_unique so secondary neighbours within 1as are removed




