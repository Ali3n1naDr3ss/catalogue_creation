# find seperation between objs det in 2 filters
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky #TODO: look at source code to see if/how Astropy make their search robust against missed nearest neighbours

path = '/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/found_by_python.fits'

table = Table.read(path)

# get coords for filt
JHra = np.array(table['RA_JH'], dtype=float)
JHdec = np.array(table['Dec_JH'], dtype=float)

Jra = np.array(table['RA_J'], dtype=float)
Jdec = np.array(table['Dec_J'], dtype=float)


# make into a 1D coord arrays
coordsJH = SkyCoord(ra=JHra * u.deg, dec=JHdec * u.deg) # all coords in MASKED cat
coordsJ = SkyCoord(ra=Jra * u.deg, dec=Jdec * u.deg)    # all coords in MASKED cat


maxsep_as  = 0.5 #as
maxsep_deg = maxsep_as / 3600

maxsep_deg = maxsep_deg * u.deg # astropy quantity
maxsep_as = maxsep_as * u.arcsec # astropy quantity

idxJ, idxJH, sep2d, _ = search_around_sky(coordsJ, coordsJH, maxsep_deg)
# idx1 the indics of the matching coords in coords1 
# coords1[idx1[0]] matches with coords2[idx2[0]]
print("number of matches: ", len(coordsJ[idxJ]), len(coordsJH[idxJH])) 
print("separations: ", sep2d.arcsec)
print("len tbale", len(JHra)) # the search_around_sky function is reporting more matches the the length of the table
print(len(sep2d.arcsec))      # there are multiple objects within a 1 as radius??


# --- Count number of matches per J source ---
# bincount will give counts of how many times each J index appears
J_match_counts = np.bincount(idxJ, minlength=len(coordsJ))

# attach the count for each row of the matches table
match_counts = J_match_counts[idxJ]
# build output table
matches_table = Table(
    [idxJ, coordsJ[idxJ].ra.deg, coordsJ[idxJ].dec.deg,
     idxJH, coordsJH[idxJH].ra.deg, coordsJH[idxJH].dec.deg,
     sep2d.arcsec, match_counts],
    names=('idxJ', 'RA_J', 'Dec_J',
           'idxJH', 'RA_JH', 'Dec_JH',
           'Sep_as', 'N_matches_within_radius'))

# --- Identify ambiguous cases ---
ambiguous = matches_table[matches_table['N_matches_within_radius'] > 1]
ambiguous_J_sources = len(np.unique(ambiguous['idxJ']))
print(f"Number of ambiguous J sources: {ambiguous_J_sources}")
#print(ambiguous[:10])

num_py_matches = 712203
num_topcta_matches = 691559
missing_from_topcat = num_py_matches - num_topcta_matches
print("missing_from_topcat: ", missing_from_topcat)
print('min sep', min(matches_table['Sep_as']))
print('max sep', max(matches_table['Sep_as']))





