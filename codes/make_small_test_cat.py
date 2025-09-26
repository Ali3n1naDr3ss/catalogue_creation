
def make_small_test_cat(catName, size=50):
    """
    Makes a new catalogue with randomly sampled rows.
    """    
    import os
    import time
    from random import randrange
    from astropy.table import Table

    baseDir = '/raid/scratch/hullyott/cataloguing/current/'
    globalCatDir = os.path.join(baseDir, 'data/catalogues/COSMOS_experimental/')

    combinedCatPath = os.path.join(globalCatDir, catName)

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

for filt in ['J', 'JH', 'H']:   
    make_small_test_cat(f'det_{filt}/COSMOS_DR2_MASKVISTADET_{filt}_1.8as_IRAC2.8as.fits', 2000)





