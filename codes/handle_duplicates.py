from astropy.table import Table, unique
from astropy.io import fits
import numpy as np
import os

def remove_duplicates():
    # Input and output file paths
    infile = "/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/topcat_J_not_JH.fits"
    outfile = "topcat_J_not_JH_full_no_duplicates.fits"

    # Read the FITS table
    table = Table.read(infile, format="fits")

    # Remove duplicate IDs (keeping the first occurrence)
    table_unique = unique(table, keys="ID")

    # Write to a new FITS file
    table_unique.write(outfile, overwrite=True)

    print(f"Removed duplicates. New table written to {outfile}")
    print(len(outfile))



def find_duplicated_ids(fits_file, id_column="ID"):
    # Open the FITS file and read the table
    with fits.open(fits_file) as hdul:
        data = hdul[1].data  # usually the table is in extension 1
        ids = data[id_column]

    # Find duplicated IDs
    unique_ids, counts = np.unique(ids, return_counts=True)
    duplicated_ids = unique_ids[counts > 1]

    return duplicated_ids


paths = ["/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/det_J/COSMOS_DR2_MASKVISTADET_J_1.8as_IRAC2.8as.fits", 
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/det_JH/COSMOS_DR2_MASKVISTADET_JH_1.8as_IRAC2.8as.fits",
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/det_H/COSMOS_DR2_MASKVISTADET_H_1.8as_IRAC2.8as.fits",
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/det_J/COSMOS_DR2_MASKVISTADET_J_1.8as_IRAC2.8as_short.fits", 
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/det_JH/COSMOS_DR2_MASKVISTADET_JH_1.8as_IRAC2.8as_short.fits",
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/det_H/COSMOS_DR2_MASKVISTADET_H_1.8as_IRAC2.8as_short.fits",
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/found_by_topcat.fits", 
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/J_distinct_from_JH.fits", 
"/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/short_test/J_distinct_from_JH_short.fits",
 "/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/short_test/not_found_by_python_short.fits"]

for path in paths:
    dupes = find_duplicated_ids(path, id_column="ID")
    print(f"Duplicated IDs in {os.path.basename(path)}:", dupes)








