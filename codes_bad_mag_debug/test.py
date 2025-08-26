from astropy.io import fits
import numpy as np

# Open FITS file
filename = ("/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/COSMOS/cutouts/UVISTA_JH_DR6_wht_1495_22_size500.fits")
fixed_filename = "/raid/scratch/hullyott/cataloguing/DepthsTestDir/data/COSMOS/cutouts/fixed.fits"

# Open file
hdul = fits.open(filename)

# Get header
hdr = hdul[0].header

# Remove all WCS-related keywords that might be invalid
for key in list(hdr.keys()):
    if key.startswith(("CRPIX", "CRVAL", "CDELT", "CTYPE", "CUNIT", "CD", "PC", "EQUINOX", "RADESYS", "LONPOLE", "LATPOLE")):
        del hdr[key]

# Save repaired file
hdul.writeto(fixed_filename, overwrite=True)
hdul.close()

print(f"Saved repaired FITS as {fixed_filename}")
