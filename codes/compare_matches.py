import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

baseDir = '/raid/scratch/hullyott/cataloguing/current/data/catalogues/COSMOS_experimental/short_test/'

# get ra dec for both 
topcat_match = os.path.join(baseDir, 'topcat_J_not_JH_short.fits')
python_match = os.path.join(baseDir, 'J_distinct_from_JH_short.fits')
python_not_match = os.path.join(baseDir, 'not_found_by_python_short.fits')

topcat_match = Table.read(topcat_match)
python_match = Table.read(python_match)
python_not_match = Table.read(python_not_match)

# make into a 1D coord arrays
ID = np.array(topcat_match['ID'], dtype=int)
ra = np.array(topcat_match['RA'], dtype=float)
dec = np.array(topcat_match['DEC'], dtype=float)

IDp = np.array(python_match['ID'], dtype=int) #p = python
rap = np.array(python_match['RA'], dtype=float)
decp = np.array(python_match['DEC'], dtype=float)

IDnp = np.array(python_not_match['ID'], dtype=int) #np = not found by python
ranp = np.array(python_not_match['RA'], dtype=float)
decnp = np.array(python_not_match['DEC'], dtype=float)

# plot

plt.scatter(ra, dec, s=1, label='topcat matches')
plt.scatter(ranp, decnp, s=3, label='not found by python')
#plt.scatter(rap, decp, s=3, label='python matches')
plt.title('J distinct from JH Comparison')
plt.legend()
plt.show()

topcat_match = os.path.join(baseDir, 'topcat_J_not_JH_full.fits')
python_match = os.path.join(baseDir, 'J_distinct_from_JH_full.fits')
python_not_match = os.path.join(baseDir, 'not_found_by_python_full.fits')
