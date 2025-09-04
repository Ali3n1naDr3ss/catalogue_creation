from astropy.io import fits
import sep
import numpy as np
image = '/raid/scratch/data/COSMOS/UVISTA_J_DR6.fits'


'''print('Performing photometry on image:',image)
with fits.open(image) as hdulist:
    data = hdulist[0].data  
    #data = data.byteswap().newbyteorder()
print(data)   
unbad_fluxes = data[data > 0]
print(len(unbad_fluxes),"/",data.size,"=", len(unbad_fluxes)/data.size)
#np.savetxt("JH_pixel_values.txt", unbad_fluxes, fmt='%.6f') '''

table = '/raid/scratch/data/catalogues/COSMOS/det_JH/COSMOS_DR2_MASKED_JH_1.8as_IRAC2.8as.fits'
from astropy.io import fits

table = '/raid/scratch/data/catalogues/COSMOS/det_JH/COSMOS_DR2_MASKED_JH_1.8as_IRAC2.8as.fits'

with fits.open(table) as hdul:
    data = hdul[1].data
mask = data['flux_JH'] > -99
print(data['flux_JH'][mask], len(data['flux_JH'][mask]),'/',data.size,'=', len(data['flux_JH'][mask])/data.size)
