
import os
import time
import numpy as np
import astropy.units as u
from datetime import datetime
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord, search_around_sky
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

baseDir = '/raid/scratch/hullyott/cataloguing/current/'
globalCatDir = os.path.join(baseDir, 'data/catalogues/COSMOS_experimental/')

# locate files
catPaths = []
for df in ['JH', 'J']:
    fDir = f'det_{df}/COSMOS_DR2_MASKVISTADET_{df}_1.8as_IRAC2.8as.fits'
    catpath = os.path.join(globalCatDir, fDir)
    catPaths.append(catpath)

# get each table
JHdir = 'det_JH/COSMOS_DR2_MASKVISTADET_JH_1.8as_IRAC2.8as.fits'
JHpath = os.path.join(globalCatDir, JHdir)

Jdir = 'det_J/COSMOS_DR2_MASKVISTADET_J_1.8as_IRAC2.8as.fits'
Jpath = os.path.join(globalCatDir, Jdir)

JHtable = Table.read(JHpath)
Jtable = Table.read(Jpath)

# get coords for filt
JHra = np.array(JHtable['RA'], dtype=float)
JHdec = np.array(JHtable['DEC'], dtype=float)

Jra = np.array(Jtable['RA'], dtype=float)
Jdec = np.array(Jtable['DEC'], dtype=float)

# Handle maksed values
mask = ~np.isnan(JHra) & ~np.isnan(JHdec)
JHra, JHdec = JHra[mask], JHdec[mask]  ## Keep only valid entries

mask = ~np.isnan(Jra) & ~np.isnan(Jdec)
Jra, Jdec = Jra[mask], Jdec[mask]  ## Keep only valid entries

# make into a 1D coord arrays
cat1 = SkyCoord(ra=JHra * u.deg, dec=JHdec * u.deg) # all coords in MASKED cat
cat2 = SkyCoord(ra=Jra * u.deg, dec=Jdec * u.deg)    # all coords in MASKED cat

# --- KD-tree search (via Astropy) ---
idx, sep2d, _ = cat1.match_to_catalog_sky(cat2)

# --- Brute force (all pairwise separations) ---
# For each source in cat1, compute angular distances to all in cat2
brute_idx = []
brute_sep = []
for c in cat1:
    sep = c.separation(cat2)
    i = np.argmin(sep)
    brute_idx.append(i)
    brute_sep.append(sep[i])

brute_idx = np.array(brute_idx)
brute_sep = u.Quantity(brute_sep)

# --- Compare results ---
exact_match = np.sum(idx == brute_idx)
accuracy = exact_match / n_cat1

# Distances can differ by tiny floating-point epsilons, so check closeness
sep_diff = (sep2d - brute_sep).to(u.arcsec)  # convert to arcseconds

print(f"Accuracy (same neighbor found): {accuracy:.2%}")
print(f"Mean |sep difference|: {np.mean(np.abs(sep_diff)):.3e}")
print(f"Max |sep difference|: {np.max(np.abs(sep_diff)):.3e}")
