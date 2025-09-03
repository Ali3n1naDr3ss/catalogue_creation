import os

def open_cats(cats, open_subset=True, overwrite=False):
    """
        Opens catalogues for filters of interest in TOPCAT.
        Also creates the subsets and returns them for future use.
        NB: has only been tested for master catalogues. Opening depth cats requires debugging.

    INPUT(s)
        cats[list(str)]     list of paths to the cats you want to open. 
                            Typically, the ones you just made by looping through image_depths.
    OUTPUT(s)
        subsets(list)       List of subset paths  """
    import numpy as np
    from astropy.table import Table
    
    # discern which type of catalogue i.e. is it a depth cat with many columns, or a master cat with fewer column names?

    command = f"topcat "
    catsToOpen = []

    for cat in cats:
        catDir = os.path.dirname(cat)
        catName = os.path.basename(cat)[:-len('.fits')]
        subsetBase = os.path.join(catDir, catName)
        if '_subset' in subsetBase: # if a subset already exists, skip
            continue
        else:
            catsToOpen.append(cat)
    files_str = " ".join(catsToOpen)
    files_str = " ".join(catsToOpen)
    command += f"{files_str} "
   
    if open_subset == False:
        command += " &"
        print(command)
        os.system(command)

    elif open_subset == True:
        # assuming subset = MAG_APER[1]>50
        subsets = []
        badratiosDict = {}
        for cat in cats:
            catDir = os.path.dirname(cat)
            catName = os.path.basename(cat)[:-len('.fits')]
            subsetPath = os.path.join(catDir, catName + '_subset.fits')

            # skip if subset already exists
            if os.path.isfile(subsetPath):
                subsets.append(subsetPath)
                subset = Table.read(subsetPath)
                table = Table.read(cat)
                badrat = len(subset)/len(table)
                if 'det_' in subsetPath: # this only works for /data/COSMOS/cataloges, not depths
                    filt = subsetPath.split('det_')[1].split('/')[0]
                badratiosDict[filt] = badrat
                continue

            # create subset
            table = Table.read(cat)
            if "MAG_APER" in table.colnames:
                subset = table[table["MAG_APER"][:, 1] > 50]  # criterion
                subset.write(subsetPath, overwrite=overwrite)
                subsets.append(subsetPath)
                badrat = len(subset)/len(table)
                badratiosDict[filt] = badrat
            else:
                # get filter name from filt path
                filt = subsetPath.split('det_')[1].split('/')[0]
                subset = table[table[f"flux_{filt}"][:] < 0]
                subset.write(subsetPath, overwrite=overwrite)
                subsets.append(subsetPath)
                badrat = len(subset)/len(table)
                badratiosDict[filt] = badrat

        files_str = " ".join(subsets)
        command += f"{files_str} &"
        #print(command)
        #os.system(command) #TODO: turn on

        for key in badratiosDict.keys():
            # prints out nicely
            if len(key)%len('HSC_X') == 0:
                filt = key
            elif len(key)%len('HSC_X') == 1:
                filt = key + '    '
            elif len(key)%len('HSC_X') == 2:
                filt = key + '   '
            print("Filter: ", filt, "bad-ratio: ", round(badratiosDict[key], 3))
            #filt = subsetPath.split('det_')[1].split('/')[0]

        return subsets



filters = reqFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']
catDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/catalogues/COSMOS/det_'

cats = []

for filt in filters:
    catDir = catDir + filt + '/'
    targetFile = f'COSMOS_DR2_MASKVISTADET_{filt}_1.8as_IRAC2.8as_cgs.fits'
    catPath = os.path.join(catDir, targetFile)
    if os.path.exists(catPath):
        #print("Found: ", catPath)
        cats.append(catPath)
    catDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/catalogues/COSMOS/det_'
if len(cats) == len(filters):
    open_cats(cats, open_subset=True)
else:
    print(">>>>>>> WARNING: did not find corresponding catalogue for each of the filters. Found: \n", cats) 



