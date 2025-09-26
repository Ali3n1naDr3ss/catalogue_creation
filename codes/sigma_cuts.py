# do sigma cuts 
# 24/09/2025 18:14

import os
import time
import numpy as np
import pandas as pd
import astropy.units as u
from datetime import datetime
from astropy.table import Table, vstack

# 5-sigma detections should be based on global depths
# 2-sigma non-detections should be based on local depths

# global depths for each filter are stored here: /raid/scratch/hullyott/cataloguing/current/depths/COSMOS/results

# locate master z8 file
# z=8 must have a 5 sigma in J or H and JH, HK? #TODO: check JH, HK criteria
# calculate SNR with signal==filt_magnitude, noise==global_noise_1.8as
# determine wheter SNR > 5
# write new cat with only the SNR > 5 == True detections

baseDir = '/raid/scratch/hullyott/cataloguing/current/'
#masterPath = os.path.join(baseDir, 'data/catalogues/COSMOS_experimental/z8_master.fits')
# TODO: use python pipeline master cat when the bug is fixed. for now, use the catalogue found by topcat
# topcat master cat
# (J not in JH) = unique J
# (H not in J)
# (H not in JH)
# ((H not in J) concat (H not in JH))
# (H not in ((H not in J) concat (H not in JH))) = unique H
# JH concat unique J concat unique H

depthBaseDir = os.path.join(baseDir, 'depths/COSMOS/results/')

# get depths in JH, J, H, HK
# find SNR in each
# determine 5 sigma threshold
# determine if it survieves 
# write new 


def find_SNR(path, depthBaseDir):

    filt = 'JH'

    table = Table.read(path)
    signal = table[filt]

    noisePath = depthBaseDir + filt + '_300.txt'


    with open(noisePath) as noise:
        for row_idx, row in enumerate(noise):
            if '1.8as' in row:
                cols = row.split('\t')
                for i in cols:
                    if i == '1.8as':
                        col_idx = cols.index(i)

            if 'full' in row and 'global' in row:
                this_row = row_idx
                print(row_idx, ',', col_idx)

    noise = pd.read_csv(noisePath, delim_whitespace=True)
    noise = noise.loc[this_row, "1.8as"]
    
########## RETURN HERE ######################
    snr = signal/noise # brighter is smaller
    # have to convert mags to flux density? or use flux counts (which are in the master z=8 table)? probs not counts bc how do you convert global noise into counts??? - > convert JH mags to f_nu and global depth to f_nu
########## RETURN HERE ######################
#    table = Table.read(fileDir+filter_name+file_suffix+'.txt', format='ascii.commented_header', delimiter='\t')
#    table = table['ap', '1.8as','type'] # only include aperture size of 1.8 arcsec

    print(snr)
    print(len(snr))
    # mask out snr<5
    mask = snr<5
    print(snr[mask])

find_SNR(masterPath, depthBaseDir)


