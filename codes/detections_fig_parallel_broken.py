import os

def get_fits_data(fitsImg, verbose=True):
    """ Returns data, header, and wcs of a fits image. """

    from astropy.io import fits
    from astropy.wcs import WCS
    from datetime import datetime

    if verbose:
        print(f'\n>>>>>>> Getting data for {fitsImg} ', datetime.now().strftime("%d-%m-%Y %H:%M:%S"), "<<<<<<<\n")
    with fits.open(fitsImg) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

    return data, header, wcs

def detections_fig(imagePaths, subsetPaths, overwrite=True, verbose=True):
    """Plots subsets in RA/Dec plane with multiple images per page, saves to PDF."""

    import os
    import math
    import matplotlib
    import numpy as np
    matplotlib.use("TkAgg")  # or "QtAgg"
    from astropy.io import fits
    from datetime import datetime
    from astropy.table import Table
    import matplotlib.pyplot as plt
    import matplotlib.backends.backend_pdf
    from matplotlib.ticker import ScalarFormatter
    from astropy.visualization import ZScaleInterval

    from open_cats import call_open_cats

    figName = 'location_bad_fluxes.pdf'
    figDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/catalogues/COSMOS_full/'
    figPath = os.path.join(figDir, figName)

    if overwrite and os.path.exists(figPath):
        os.remove(figPath)

    pdf = matplotlib.backends.backend_pdf.PdfPages(figPath)
    zscale = ZScaleInterval()
    formatter = ScalarFormatter(useMathTeimport os

def get_fits_data(fitsImg, verbose=True):
    """ Returns data, header, and wcs of a fits image. """

    from astropy.io import fits
    from astropy.wcs import WCS
    from datetime import datetime

    if verbose:
        print(f'\n>>>>>>> Getting data for {fitsImg} ', datetime.now().strftime("%d-%m-%Y %H:%M:%S"), "<<<<<<<\n")
    with fits.open(fitsImg) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

    return data, header, wcs         # TODO: parallelise?

def prepare_plot_data(args):
    from astropy.table import Table

    imagePath, subsetPath, verbose = args

    imageData, _, wcs = get_fits_data(imagePath, verbose=verbose)

    # read subset

    subsetTable = Table.read(subsetPath)

    subRA = subsetTable['RA']
    subDec = subsetTable['DEC']

    # convert RA/Dec -> pixels
    print("world to pixel: ", subsetPath)
    x_pix, y_pix = wcs.world_to_pixel_values(subRA, subDec)
    print("finished world to pixel: ", subsetPath)
    return imagePath, subsetPath, imageData, x_pix, y_pix

def detections_fig(imagePaths, subsetPaths, overwrite=True, verbose=True):
    """Plots subsets in RA/Dec plane with multiple images per page, saves to PDF."""

    import os
    import math
    import matplotlib
    import numpy as np
    matplotlib.use("TkAgg")  # or "QtAgg"
    from astropy.io import fits
    from datetime import datetime
    from astropy.table import Table
    import matplotlib.pyplot as plt
    import matplotlib.backends.backend_pdf
    from matplotlib.ticker import ScalarFormatter
    from astropy.visualization import ZScaleInterval
    from concurrent.futures import ProcessPoolExecutor

    from open_cats import call_open_cats

    figName = 'location_bad_fluxes.pdf'
    figDir = '/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/catalogues/COSMOS_full/'
    figPath = os.path.join(figDir, figName)

    if overwrite and os.path.exists(figPath):
        os.remove(figPath)

    pdf = matplotlib.backends.backend_pdf.PdfPages(figPath)
    zscale = ZScaleInterval()
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    _, badratiosDict = call_open_cats(open_subset=True)

    # prepare args for pool
    args_list = [(ip, sp, verbose) for ip, sp in zip(imagePaths, subsetPaths)]

    print(">>>>>>> Loading data in parallel...")
    with ProcessPoolExecutor(max_workers=20) as pool:
        results = list(pool.map(prepare_plot_data, args_list))
    print(">>>>>>> Finished loading data.")

    print(">>>>>>> Making plots...")

    n_per_page = 4
    n_pages = math.ceil(len(results) / n_per_page)
    nrows = int(n_per_page/2)
    ncols = int(n_per_page/2)

    for page in range(n_pages):
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 8))
        axes = axes.flatten()

    for i in range(n_per_page):
        idx = page * n_per_page + i
        if idx >= len(results):
            axes[i].axis('off')
            continue

        imagePath, subsetPath, imageData, x_pix, y_pix = results[idx]

        plotTitle = os.path.basename(imagePath).replace('.fits', '')
        if 'HSC_' in plotTitle:
            plotTitle = plotTitle[:len('HSC_X')]
            badrat = str(round(round(badratiosDict[plotTitle], 3) * 100, 3)) + '%'
        else:
            plotTitle = plotTitle.split('MASKVISTADET_')[0].split('_')[1]
            badrat = str(round(round(badratiosDict[plotTitle], 3) * 100, 3)) + '%'

        ax = axes[i]
        vmin, vmax = zscale.get_limits(imageData)
        ax.imshow(imageData, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.scatter(x_pix, y_pix, color='green', alpha=0.5, s=1)
        ax.set_title(f"Filter: {plotTitle}", fontsize=10)
        ax.text(0.05, 0.95, badrat, transform=ax.transAxes,
                fontsize=18, color="red", va="top")
        ax.grid(True)
        ax.xaxis.set_major_formatter(formatter)
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("Dec (deg)")
        ax.axis('off')

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    print(">>>>>> Finished making detection figures.: ", figPath)

def call_detections_fig(
        filters =['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y'], 
        catBaseDir='/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/catalogues/COSMOS_full/det_', 
        imageDir='/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/COSMOS/'):

    cats = []
    imagePaths = []
    subsetPaths = []

    for filt in filters:
        catDir = catBaseDir + filt + '/'
        targetFile = f'COSMOS_DR2_MASKVISTADET_{filt}_1.8as_IRAC2.8as_cgs.fits'
        subsetFile = f'COSMOS_DR2_MASKVISTADET_{filt}_1.8as_IRAC2.8as_cgs_subset.fits'
        catPath = os.path.join(catDir, targetFile)
        subsetPath = os.path.join(catDir, subsetFile)

        if os.path.exists(catPath):
            cats.append(catPath)
        if os.path.exists(subsetPath):
            subsetPaths.append(subsetPath)

        if 'HSC_' in filt:
            imagePath = imageDir + filt + '_DR3.fits'
            if os.path.exists(imagePath):
                imagePaths.append(imagePath)
        else:
            imagePath = imageDir + 'UVISTA_' + filt + '_DR6.fits'
            if os.path.exists(imagePath):
                imagePaths.append(imagePath)

    detections_fig(imagePaths, subsetPaths, overwrite=True)


call_detections_fig()


t=False)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    _, badratiosDict = call_open_cats(open_subset=True)

    print("\n>>>>>>> Making detections figures. Started at: ", datetime.now().strftime("%d-%m-%Y %H:%M:%S"),    "<<<<<<<")

    n_per_page = 4
    n_pages = math.ceil(len(imagePaths) / n_per_page)
    nrows = int(n_per_page/2) 
    ncols = int(n_per_page/2)
    for page in range(n_pages):
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 8))
        axes = axes.flatten()

        for i in range(n_per_page):
            idx = page * n_per_page + i
            if idx >= len(imagePaths):
                axes[i].axis('off')  # turn off empty subplots
                continue

            imagePath = imagePaths[idx]
            subsetPath = subsetPaths[idx]

            plotTitle = os.path.basename(imagePath).replace('.fits', '')
            if 'HSC_' in plotTitle:
                plotTitle = plotTitle[:len('HSC_X')]
                badrat = round(round(badratiosDict[plotTitle], 3)*100, 3)
                badrat = str(badrat) + '%'
            else:
                plotTitle = plotTitle.split('MASKVISTADET_')[0].split('_')[1]
                badrat = round(round(badratiosDict[plotTitle], 3)*100, 3)
                badrat = str(badrat) + '%'

            step = 1
            imageData, _, wcs = get_fits_data(imagePath, verbose=verbose)

            subsetTable = Table.read(subsetPath)
            subRA = subsetTable['RA']
            subDec = subsetTable['DEC']

            # Convert RA/Dec -> pixel coordinates
            x_pix, y_pix = wcs.world_to_pixel_values(subRA, subDec)

            ax = axes[i]
            vmin, vmax = zscale.get_limits(imageData)
            ax.scatter(x_pix[::step], y_pix[::step], color='green', alpha=0.5, s=1)
            ax.set_title(f"Filter: {plotTitle}", fontsize=10)
            ax.text(150, 2.2, badrat, fontsize=18, color="red")
            ax.grid(True)
            ax.xaxis.set_major_formatter(formatter)
            ax.set_xlabel("RA (deg)")
            ax.set_ylabel("Dec (deg)")
            ax.axis('off')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    pdf.close()
    print(">>>>>> Finished making detection figures.: ", figPath)

def call_detections_fig(
        filters =['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y'], 
        catBaseDir='/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/catalogues/COSMOS_full/det_', 
        imageDir='/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/COSMOS/'):

    cats = []
    imagePaths = []
    subsetPaths = []

    for filt in filters:
        catDir = catBaseDir + filt + '/'
        targetFile = f'COSMOS_DR2_MASKVISTADET_{filt}_1.8as_IRAC2.8as_cgs.fits'
        subsetFile = f'COSMOS_DR2_MASKVISTADET_{filt}_1.8as_IRAC2.8as_cgs_subset.fits'
        catPath = os.path.join(catDir, targetFile)
        subsetPath = os.path.join(catDir, subsetFile)

        if os.path.exists(catPath):
            cats.append(catPath)
        if os.path.exists(subsetPath):
            subsetPaths.append(subsetPath)

        if 'HSC_' in filt:
            imagePath = imageDir + filt + '_DR3.fits'
            if os.path.exists(imagePath):
                imagePaths.append(imagePath)
        else:
            imagePath = imageDir + 'UVISTA_' + filt + '_DR6.fits'
            if os.path.exists(imagePath):
                imagePaths.append(imagePath)

    detections_fig(imagePaths, subsetPaths, overwrite=True)


call_detections_fig()



