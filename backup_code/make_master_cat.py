# Process catalogues to make final Master/prov0 cat
from check_for_duplicates import check_for_duplicate_IDs
from update_catalogues import update_catalogues
import numpy as np

# global definitions
baseDir = '/raid/scratch/'
catDir = 'data/catalogues/COSMOS_OLD_COPY/'

original_cat = 'cat_for_testing_on.fits'
cats_to_check = [baseDir+catDir+original_cat]

det_filters = ['JH','H','Y']
stacked_images=['JH','HKs']

# run programme
dup_instances = check_for_duplicate_IDs(cats_to_check, cat_dir=baseDir+catDir) 
# dup_instances is an arr where each element represents the number of times each duplicate is found

if len(dup_instances) == 0:
    print("No duplicate IDs found. Updating catalogues...")
    update_catalogues(det_filters, original_cat, stacked_images, 'limiting_depths.txt', base_dir=baseDir, cat_dir=catDir, quick_test=True, lim_depths_dir=baseDir+'depths/COSMOS/results/', input_dir=baseDir+catDir, output_dir=baseDir+catDir)


elif np.any(dup_instances) == 0:
        print(f"{dup_instances} duplicates found...")

