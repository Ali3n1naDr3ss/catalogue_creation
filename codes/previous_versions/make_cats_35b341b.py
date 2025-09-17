#!/usr/bin/env python3

""" make_catalogue.py

Aim: This code is a wrapper for my other codes, with the goal of doing a high-z selection for all the images available.

Created: Thursday 2nd Nov 2017
Modified: from nake_catalogue.py 19/08/2025

"""

###################### Import useful modules #########################
import numpy as np
from datetime import datetime
from cat_codes_35b341b import *
from joblib import Parallel, delayed

se = False           # run source extractor in dual image mode
check = False       # check all the catalogues have been run and have the same length
mask = False        # mask these catalogues
fluxes = True      # convert into fluxes and errors
zptshift = False    # correct the zeropoints slightly (as derived from NA fitting)
combine = False     # combine 123 regions into one big catalogue, removing duplicates

######################### Set-up ####################################
detectionFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']

fields = ['COSMOS']

# which filters to include in dual image mode
reqFilters = ['Y', 'J', 'H', 'K', 'JH', 'HK', 'HSC_G', 'HSC_R', 'HSC_I', 'HSC_Z', 'HSC_Y']

n_jobs = len(detectionFilters)
# either set this to ['all'] to run all filters present in images.lis

# required aperture diameter
#reqMag = '2.0' # or 'MAG_AUTO' etc.

queue = 'none'  #'none'
overwrite = False # False should be default

####################################################################
# required aperture diameters to run through
# change at your peril!
apDiametersAS = [1.8, 2.0, 3.0, 4.0, 5.0]
IRACapDiametersAS = [2.8, 3.8, 5.8, 9.8, 11.6]

# these are the apertues that are included in the final full catalogue
requiredCats =  np.array(['1.8', '2.0', 'MAG_AUTO'])#, 'FLUX_RADIUS'])
requiredSpit = np.array(['2.8', '2.8', 'NONE'])#, 'NONE'])

# just run the first catalogue for now
requiredCats = requiredCats[[0,2]]
requiredSpit = requiredSpit[[0,2]]

# combine the three catalogues into one
field = fields[0][:-1]

if field == 'XMM':
    combinedDir = '../catalogues/XMMFULL/'
elif field == 'CDFS':
    combinedDir = '../catalogues/CDFSFULL/'

# the final name for the catalogue
#basename = '_DR2_MASKVISTA_'
#basename = '_DR2_UNMASKED_'
basename = '_DR6_'

# Different basename for HSC SSP DR3.
#basename = '_DR3_MASKVISTADET_'

if combine:
    combine_xmmcdfs(detectionFilters[0], basename, requiredCats,field, outDir = combinedDir, requiredSpit = requiredSpit) #, cuts = [23.8, -99.0, -99.0])
    exit()

############################### Loop ################################
## Loop through the different fields
# first do the source selection
for ff, fieldName in enumerate(fields):

    print('#############################################')
    print("Analysing field ", fieldName)
    # Run the function
    if se:
        run_all_source_extractors(fieldName=fieldName,
            detectionFilters=detectionFilters,
            apDiametersAS=apDiametersAS, queue="none", memory=7, reqFilters=reqFilters,
            IRACapDiametersAS=IRACapDiametersAS, overwrite=overwrite, se=True,
            n_jobs=len(detectionFilters))

if se:
    print("Finished se = True at: ",  datetime.now().strftime("%d/%m %H:%M:%S"), '\n')
    exit()

def run_masking(fields, detectionFilters, requiredCats, requiredSpit,
                baseDir='/raid/scratch/hullyott/cataloguing/DepthsTestDir/',
                check=True, mask=True, n_jobs=n_jobs):

    jobs = []  # for parallel apply_mask calls

    for ff, fieldName in enumerate(fields):
        print('#############################################')
        print("Analysing field ", fieldName)

        for df, detectionFilt in enumerate(detectionFilters):

            if check:
                print("check = True", '\n')
                # check everything is in place
                check_all_files(fieldName, detectionFilt)

            if mask:
                print("mask = True", '\n')
                combine_cats(fieldName, detectionFilt, requiredCats,
                             requiredSpit=requiredSpit,
                             imageDir=baseDir + 'data/')

                # Now collect all apply_mask jobs
                for ai, apD in enumerate(requiredCats):
                    if requiredSpit[ai] == 'NONE':
                        spitstring = ''
                        optstring = apD.replace('.','_')
                        
                    else:
                        requiredSpit[ai].replace('.','_')
                        spitstring = '_IRAC{0}as'.format(requiredSpit[ai])
                        optstring = '{0}as'.format(apD.replace('.','_'))

                    # 1st: unmasked ? masked
                    inputCat = (
                        baseDir + 'data/catalogues/COSMOS/det_{0}/'
                        'COSMOS_DR2_UNMASKED_{0}_{2}{3}.fits'
                    ).format(detectionFilt, fieldName, optstring, spitstring)
                    jobs.append(delayed(apply_mask)(inputCat, fieldName, save=True, hsc=True))

                    # 2nd: apply K cut
                    inputCat = (
                        baseDir + 'data/catalogues/COSMOS/det_{0}/'
                        'COSMOS_DR2_MASKED_{0}_{2}{3}.fits'
                    ).format(detectionFilt, fieldName, optstring, spitstring)
                    jobs.append(delayed(apply_mask)(inputCat, fieldName, save=True, hsc=True, cutband='K'))

                    # 3rd: if detection filter isn't K, also cut by detection filter
                    if detectionFilt != 'K':
                        inputCat = (
                            baseDir + 'data/catalogues/COSMOS/det_{0}/'
                            'COSMOS_DR2_MASKVISTA_{0}_{2}{3}.fits'
                        ).format(detectionFilt, fieldName, optstring, spitstring)
                        jobs.append(delayed(apply_mask)(inputCat, fieldName, save=True, cutband=detectionFilt, hsc=True))

                print("masking jobs queued for", fieldName, detectionFilt)

    # Run all masking jobs in parallel
    if jobs:
        print(f"Running {len(jobs)} masking jobs in parallel with n_jobs={n_jobs}...")
        results = Parallel(n_jobs=n_jobs, verbose=10)(jobs)
        return results

    return None

if mask:
    print("Beginning mask = True at: ",  datetime.now().strftime("%d/%m %H:%M:%S"), '\n')
    results = run_masking(
    fields, detectionFilters, requiredCats, requiredSpit,
    baseDir='/raid/scratch/hullyott/cataloguing/DepthsTestDir/',
    check=True, mask=True, n_jobs=n_jobs)  
                
if check | mask:
    exit()    
    
# Fluxes!
if fluxes:
    print("Beginning fluxes = True at: ",  datetime.now().strftime("%d/%m %H:%M:%S"), '\n')
    print("fluxes = True", '\n')
    #if detectionFilt != 'Ks':
    catString = 'DR2_MASKVISTADET_'
    #else:
        #catString = 'DR2_MASKVISTA_'

  #  XMM1_DR2_MASKED_HSC-G_1.8as_206testpsf.fits
        
    spawn_flux_errors(detectionFilters, fields, requiredCats, requiredSpit, queue = queue, catString = catString, short_test=False)
print("Finished fluxes = True at: ",  datetime.now().strftime("%d/%m %H:%M:%S"), '\n')
"""
HOLLY 28/03/25
Uses SE python module, circle_sum, to perform aperture 
photometry using subpixel sampling to improve accuracy
(how does this work?). circle_sum returns the error of
the sum (by doing what? see sep files). Where the
returned flux is negative, fluxes are re-cast as 1.0.
Next, fluxes are converted to magnitudes using the
standard relation and provided zeropoint from images.lis.
Bad pixels are identified and replaced with zero flux,
or large magnitudes (99). Fluxes are converted to erg/s/cm2/Hz
or microJankskys. Minimum errors are enforced so that
fluxes remain positive. Final flux and error are scaled by
the flux enclosed by the aperture, this ensures total flux
of the object is calculated while minimising the effects
of noise. 
"""



#if zptshift:

    # I want to apply the zeropoint shifts to
    # 1) the fluxes and errors as a multiplicative factor
    
    
    # 2) to the mag file for clarity of what is going on!
    

    
