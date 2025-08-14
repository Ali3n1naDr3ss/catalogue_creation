

"""
cut large catalogue into more manageable sizes according to cutout-size
"""
import os
import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord


baseDir = '/raid/scratch/'
############ Cut out science image ############ 
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import os

def make_cutout(image_path, center_ra, center_dec, size_arcsec, output_path=None):
    with fits.open(image_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

        center = SkyCoord(center_ra, center_dec, unit='deg')

        try:
            scale = abs(header['CD1_1']) * 3600  # deg to arcsec
        except KeyError:
            scale = abs(header['CDELT1']) * 3600

        size_pix = int(size_arcsec / scale)

        cutout = Cutout2D(data, center, size_pix, wcs=wcs)

        # Copy the original header and update it with the cutout WCS
        new_header = header.copy()
        new_header.update(cutout.wcs.to_header())

        hdu = fits.PrimaryHDU(data=cutout.data, header=new_header)
        hdul_out = fits.HDUList([hdu])

        if output_path is None:  
            out = os.path.splitext(image_path)[0][:-len('_backup')]
            output_path = out + '.fits'

        hdul_out.writeto(output_path, overwrite=True)
        print(f"Saved cutout to {output_path}")
    return output_path

############ chop cat accordingly ############ 
def get_fits_data(fits_img):
    with fits.open(fits_img) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
    return data, header, wcs

def chop_cat(sciCutout, cutout_size, filters=['J', 'JH'], catDir=baseDir+'depths/COSMOS_original_run/catalogues/', outputDir=baseDir+'depths/COSMOS_test/catalogues/'):
    """ Writes a cutout of the original, full-size depth catalogue(s) based on the size and coordinates of the science image.
    INPUTs
        sciCutout(str)      Full path to the science image cutout you are testing on
        cutout_size(int)    Length in arcsec of the square science cutout 
        filters(list[str])  Filter names NB: only tested on UVISTA
        catDir(str)         Dir to the original, full-size depth catalogues
        outputDir(str)      Dir to save cutout catalogues to i.e. the test folder

    OUTPUT
        outputCat(str)      Full path to cutout catalogue
    """

    # get bounds of science image cutout so the table doesn't go beyond this 
    data, header, wcs = get_fits_data(sciCutout)
    n_dec, n_ra = data.shape # shape = (rows, cols)

    # get pix scale
    try:
        scale = abs(header['CD1_1']) * 3600  # deg to arcsec
    except KeyError:
        scale = abs(header['CDELT1']) * 3600

    size_pix = int(cutout_size / scale)

    ## get limits from original image
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    # Center of cutout from WCS
    center_pixel = [data.shape[1] / 2, data.shape[0] / 2]
    center_sky = wcs.pixel_to_world(*center_pixel)

    center = SkyCoord(ra=center_sky.ra, dec=center_sky.dec)

    for filt in filters:
        inputCat = Table.read(catDir + 'd' + filt + '.fits')
        cat_coords = SkyCoord(ra=inputCat['ALPHA_J2000'], dec=inputCat['DELTA_J2000'], unit='deg')
        
        # Keep only sources within cutout radius
        sep = cat_coords.separation(center)
        radius = (cutout_size / 2.0) * u.arcsec
        mask = sep <= radius

        choppedCat = inputCat[mask]

        outputCat = outputDir + 'd' + filt + "_cutout.fits"
        choppedCat.write(outputCat, overwrite=True)
        print(f"Saved: {outputCat}")
    return outputCat

def plot(sciCutout, outputCat):
    """ broken  """
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm

    data, header, wcs = get_fits_data(sciCutout) # needed to get corresponding ra/dec in img/depth cat
    cat = Table.read(outputCat) #TODO JH...

    fig = plt.figure()
    ax = fig.add_subplot() 
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.grid(color='white', ls='dotted')
    norm = simple_norm(data, 'sqrt', percent=99)

    ax.imshow(data, origin='lower', norm=norm, cmap='gray')
    ax.scatter(cat['ALPHA_J2000'], cat['DELTA_J2000'],
               c='green',
               s=10, edgecolor='red', 
               facecolor='none')

    plt.show()

    



