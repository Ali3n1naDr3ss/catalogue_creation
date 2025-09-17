# make smaller test cats to manually check whether matched objects are found

def make_small_test_cat(filters, size=50):
    import os
    import numpy as np  
    import astropy.units as u
    from astropy.table import Table
    from astropy.coordinates import SkyCoord

    baseDir = '/raid/scratch/hullyott/cataloguing/current/'
    globalCatDir = os.path.join(baseDir,  'data/catalogues/COSMOS_make_master_test/')

    combinedCatPath = os.path.join(globalCatDir, "combined_cat.fits")
    print("Reading large table...")
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
        shortCat.write(combinedCatPath.replace(".fits", "_short.fits"))

filters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']
make_small_test_cat(filters)

