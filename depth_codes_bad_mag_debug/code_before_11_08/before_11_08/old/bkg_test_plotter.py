import numpy as np
from astropy.table import Table
import os

baseDir = '/raid/scratch/'
depthDir = baseDir + 'depths/COSMOS_test/'
dataDir = baseDir + '/data/COSMOS_test/'
verbose = True


def bkg_test_plotter(param_combos, depthDir='/raid/scratch/depths/COSMOS_test/', dataDir='/raid/scratch/data/COSMOS_test', verbose=True):
    """
    what affect do each set of SE bkg modelling parameters have?
    Do visual check with this side-by-side plotter

    param_combos(list): SE parameters used in bkg modelling test: [BACK_SIZE, BACK_FILTERSIZE]
    depthDir(str):      directory of SE'd find_depth catalogues
    dataDir(str):       directory of science images (full-size and cutouts)
    """

    import re
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from collections import defaultdict
    from astropy.visualization import ZScaleInterval

    ImgDir = depthDir + 'images/' # bkg map and bkg-sub'd img

    ## get all filenames in ImgDir and verify given test parameters are present
    if os.listdir(testImgDir)==False:
        print(f">>>>>>>>>>>>>>>>>>>>>>>> WARNING: ImageDir NOT found: {testImgDir}")
    else:
        all_files = [f for f in os.listdir(testImgDir) if os.path.isfile(os.path.join(testImgDir, f))]
        comparison_imgs = []    
        for f in all_files:
            for param_combo in param_combos:
                if str(param_combo) in f:
                    test_imgs.append(testImgDir+f)
                    if verbose:
                        print(f"Bkg model comparison image found: {f}")
                else:
                     print(f">>>>>>>>>>>>>>>>>>>>>>>> WARNING: Bkg model comparison image NOT found: {f}")

    ## Specify which bkg params you want to compare  
    pattern = re.compile(r"^(?P<band>.+?)bk_sz(?P<bksz>.+?)bkfisz(?P<bkfsz>.+?)_(?P<type>.+?)\.fits")

    # Group files by (bksz, bkfsz)
    grouped = {}
    bands = []
    for fname in comparison_imgs:
        match = pattern.match(os.path.basename(fname)) # extracts info from filename if name matches supplied pattern
        if not match:
            continue  # Skip files that don't match the pattern
            if verbose:
                print(f"Skipping unmatched filename: {fname}")

        band = match.group("band") # filter/band
        bands.append(band)
        bksz = match.group("bksz")
        bkfsz = match.group("bkfsz")
        ftype = match.group("type")

        key = (bksz, bkfsz) # combination of test params

        if key not in grouped:
            grouped[key] = {"J": {}, "JH": {}} # add keys for subdicts 
            # TODO: make flexible in terms of filter choice
            print(f"Adding: grouped[{key}][{band}][{ftype}]")
        grouped[key][band][ftype] = os.path.join(ImgDir, fname)

    # get all filenames in science image dir and verify test bands are present
    all_files = [f for f in os.listdir(dataDir) if os.path.isfile(os.path.join(dataDir, f))]
    sci_imgs = []   
    orig_pattern = re.compile(r"^UVISTA_(?P<band>.+?)\_DR6.fits") # TODO: make flexible

    for f in all_files:
        basename = os.path.basename(f)
        match = orig_pattern.match(basename)
        if match:
            matched_band = match.group("band")
            if verbose:
                print(f"Science image found: {f}")
            if matched_band in ['J', 'JH']:
                sci_imgs.append(os.path.join(origImgDir, f))
                if verbose:
                    print(f"Visualising background model for: {f}")

    # display bgmap, bgsub, original for both filters

    # Find group with matching files for J and JH
    selected_group = None
    for (bksz, bkfsz), bands in grouped.items():
        if all(ftype in bands["J"] for ftype in ["bgmap", "bgsub"]) and \
           all(ftype in bands["JH"] for ftype in ["bgmap", "bgsub"]):
            selected_group = (bksz, bkfsz)
            break

    if not selected_group:
        raise ValueError("No group found with complete J and JH bgmap, and bgsub images.")

        bksz, bkfsz = selected_group

    ### Load bgmap and bgsub images TODO: hard-coded to J and JH only
    j_images = []
    jh_images = []
    ftypes = ["bgmap", "bgsub"]

    for ftype in ftypes:
        j_path = grouped[(bksz, bkfsz)]["J"][ftype]   # TODO: not flexible
        jh_path = grouped[(bksz, bkfsz)]["JH"][ftype] # TODO: not flexible

        j_data = fits.getdata(j_path)
        jh_data = fits.getdata(jh_path)

        j_images.append((j_data, ftype))
        jh_images.append((jh_data, ftype))

    ### Load science images
    sciImg_pattern = re.compile(r"^UVISTA_(?P<band>.+?)_DR6\.fits")

    for img_path in sci_imgs:
        basename = os.path.basename(img_path)
        match = orig_pattern.match(basename)
        if not match:
            continue

        band = match.group("band")
        data = fits.getdata(img_path)

        if band == "J":
            j_images.append((data, "original"))
        elif band == "JH":
            jh_images.append((data, "original"))

    # TODO: add check that science image corresponds to the bkg e.g. cutout coords
    
    ### Plot
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    zscale = ZScaleInterval()

    ## Top row: J-band
    for i, (data, ftype) in enumerate(j_images):
        vmin, vmax = zscale.get_limits(data)
        ax = axes[0, i]
        ax.imshow(data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"J bk_sz{bksz} bkfisz{bkfsz} {ftype}", fontsize=10)
        ax.axis('off')

    # Bottom row: JH-band
    for i, (data, ftype) in enumerate(jh_images):
        vmin, vmax = zscale.get_limits(data)
        ax = axes[1, i]
        ax.imshow(data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(f"JH bk_sz{bksz} bkfisz{bkfsz} {ftype}", fontsize=10)
        ax.axis('off')

    plt.tight_layout()
    plt.show()
                  

