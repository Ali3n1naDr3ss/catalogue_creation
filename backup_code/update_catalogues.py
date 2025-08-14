# Aim: To perform various updates to the relevant catalogue(s) e.g. adding a SNR column
# Based on V1.2
# This version aims to handle changing multiple catalogues at once
# 26/05/2025 17:17
from astropy.io import fits
import numpy as np
def make_data_dirs(det_filters, cat_dir=None, base_dir=None):
    """ Creates a list of directories for the magnitude catalogue for each detection filter."""
    if cat_dir == None or base_dir == None:    
        raise ValueError(">>>>>>>>>>>>>>>>>>>>>>> ERROR: cat_dir or base_dir not supplied")    
    cats_with_mags = []
    for filt in det_filters:
        cats_with_mags.append(base_dir+cat_dir+'det_'+filt+'/COSMOS_DR2_MASKED_'+filt+'_1.8as_IRAC2.8as.fits')
    return cats_with_mags

def open_fits(cat_dir, quick_test=True):
    print(">>>>>>>>>>>>>>>>>>>>>>> ", cat_dir)
    with fits.open(cat_dir) as hdul:
        data = hdul[1].data
        hdu = hdul[0]
        if quick_test:
            data = data[:20]  # Test on a short catalog (only 20 rows)
            print("\n ######################### PERFORMING A QUICK TEST ######################### \n")
            print(f"Opening catalogue: {cat_dir} \n")
    # TODO: what about when multiple cats in list 
    return data, hdu

def rename_cat_cols(cat_data):
    """ Renames columns in original catalogues according to the mapping defined here."""

    print(f"Renaming catalogues...\n")
    print("WARNING: Check the number-filter name mapping is the same for new cats 28/05/25")
    # mapping from original suffix (automatically assigned by SExtractor) to new suffix (detection filter)
    rename_map = {
        '_1': '_JH',
        '_2': '_HK',
        '_1a': '_J',
        '_2a': '_H'}
    orig_cols = cat_data.columns
    new_cols = []

    # do re-naming
    for col in orig_cols:
        name = col.name
        for orig_suffix, new_suffix in rename_map.items():
            if name.endswith(orig_suffix):
                name = name[:-len(orig_suffix)] + new_suffix
                break
        # Create new column with the updated name but same format and data
        new_cols.append(fits.Column(name=name, format=col.format, array=col.array))
    
    return new_cols
    
def calc_detection_significance(cat_with_mags, limiting_depths_filename, image_names, output_cat, lim_depths_dir=None, cat_dir=None, quick_test=True):
    """ Calculates the detection significance (SNR) for each object in the catalogue based on the limiting depths.
    Returns a list of new columns to be added to the catalogue."""
    if lim_depths_dir == None or cat_dir == None:    
        raise ValueError(">>>>>>>>>>>>>>>>>>>>>> lim_depths_dir or cat_dir not supplied")
    print(f"Calculating SNR for: {cat_with_mags}")    
    # retrieve cat col names for comparison with rows in limiting_depths.txt
    data, hdu  = open_fits(cat_with_mags, quick_test=quick_test)
    col_names = data.columns.names
    SNRs_for_stack = []
    snr_cols = []
    
    # do comparison
    with open(lim_depths_dir+limiting_depths_filename, "r") as depth_file:
        for line in depth_file:
            for stack in image_names:
                if stack in line:
                    for i in enumerate(col_names):
                        if stack == col_names[i[0]]:
                            line_stack = line.strip().split(' ')[0]  # don't overwrite `stack`
                            if stack != line_stack:
                                continue  # continue for instances such as CFHT: we don't have these mags in current cats
                            lim_range_str = line.strip().split(' ')[1]
                            if '-' in lim_range_str:
                                lim_range_0, lim_range_1 = lim_range_str.split('-')[0],lim_range_str.split('-')[1]
                                lim_range_0, lim_range_1 = float(lim_range_0), float(lim_range_1)
                                low_lim = min(lim_range_0, lim_range_1)
                            else:
                                low_lim = float(lim_range_str)
                            mags = data[stack] # measured magnitudes in cat
                            for mag in mags:
                                mag = float(mag)
                                # find SNR, detection level
                                # see: https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/photo.html
                                # m1 - m2 = -2.5log10(f1/f2), let m2, f2 = limiting mag, flux
                                # SNR = f/f_lim
                                SNR = 10**(0.4*(low_lim - mag))
                                SNRs_for_stack.append(SNR)
                            SNR_arr = np.array(SNRs_for_stack)
                            snr_cols.append(fits.Column(name=f"SNR_{stack}", format='D', array=SNR_arr))
    return snr_cols

def update_cat(cols_to_add, input_cat, output_name=None, input_dir=None, output_dir=None, overwrite=True, combine=True):
    """ Updates the input catalogue by adding new columns defined in `cols_to_add`."""
    if input_dir == None or output_dir == None:    
        raise ValueError(">>>>>>>>>>>>>>>>>>>> ERROR: input_dir or output_dir not supplied")
    # Open the original catalogue and extract columns
    with fits.open(input_dir+input_cat) as hdul:
        orig_data = hdul[1].data
        orig_cols = fits.ColDefs(orig_data) # definitions class

    if combine:
        # combine cols-to-add with original cols
        new_col_names = [col.name for col in cols_to_add]
        # Filter out any original columns with conflicting names
        filtered_orig_cols = [col for col in orig_cols if col.name not in new_col_names]
        # Combine columns
        all_cols = filtered_orig_cols + cols_to_add
    else:
        all_cols = cols_to_add

    # Create new HDU with new columns
    new_hdu = fits.BinTableHDU.from_columns(cols_to_add)
    new_col_names = new_hdu.data.columns.names
    # make new HDU for updated cols
    if output_name is None:
       output_dir = output_dir + input_cat.replace('.fits', '_renamed.fits')
    else:
       output_dir = output_dir + input_cat.replace(input_cat, output_name)

    print("Creating new HDU...")
    new_hdu = fits.PrimaryHDU()  # blank primary HDU
    table_hdu = fits.BinTableHDU.from_columns(all_cols)
    hdulist = fits.HDUList([new_hdu, table_hdu])

    print("Writing updated catalogue...")
    hdulist.writeto(output_dir, overwrite=overwrite)
    print(f"Saved updated catalogue to: {output_dir}")

def update_catalogues(det_filters, original_cat, stacked_images=['JH','HKs'], limiting_depths_filename='limiting_depths.txt', base_dir=None, cat_dir=None, lim_depths_dir=None, input_dir=None, output_dir=None, quick_test=True):
    """Runs the functions in the update_catalogues programme to create new, updated catalogues."""

    mag_dirs = make_data_dirs(det_filters, cat_dir, base_dir)
    original_data, original_hdu = open_fits(original_cat)
    new_cols = rename_cat_cols(original_data)
    update_cat(new_cols, original_cat, 'cat_for_testing_on_renamed.fits', input_dir=base_dir+cat_dir, output_dir=base_dir+cat_dir, combine=False)
    snr_cols = calc_detection_significance(mag_dirs[0], limiting_depths_filename, det_filters, 'cat_for_testing_on_renamed.fits', lim_depths_dir=lim_depths_dir, cat_dir=base_dir+cat_dir, quick_test=quick_test)
    update_cat(snr_cols, 'cat_for_testing_on.fits', 'cat_for_testing_on_renamed.fits', input_dir=base_dir+cat_dir, output_dir=base_dir+cat_dir) # first need to use renamed cat

    # TODO: this will probs be a very slow way to handle so much data - consider batches or alternative method    
    for mag_dir in mag_dirs:
        snr_cols = calc_detection_significance(mag_dir, limiting_depths_filename, det_filters, 'cat_for_testing_on_renamed.fits', lim_depths_dir=lim_depths_dir, cat_dir=cat_dir)
        update_cat(snr_cols, 'cat_for_testing_on_renamed.fits', 'cat_for_testing_on_SNR.fits', input_dir=input_dir, output_dir=output_dir)


#update_catalogues(['JH','H','Y'], ['JH','H'])

# run the stuff

#det_filters = ['JH','H','Y']
#renamed_cat = 'cat_for_testing_on_renamed.fits'
#limiting_depths = 'limiting_depths.txt'
#stacked_images = ['JH','H'] # H for testing on old copies (JH and HK mags are broken)
#quick_test = True

#mag_dirs = make_data_dirs(det_filters)

#original_data, original_hdu = open_fits(original_cat)
#new_cols = rename_cat_cols(original_data, original_hdu, cat_name=original_cat)
#update_cat(new_cols, original_cat, renamed_cat, combine=False)

#snr_cols = calc_detection_significance(mag_dirs[0], limiting_depths, det_filters, renamed_cat)
#update_cat(snr_cols, 'cat_for_testing_on_renamed.fits', output_name='cat_for_testing_on_SNR.fits') # first need to use renamed cat

#for mag_dir in mag_dirs:
    #snr_cols = calc_detection_significance(mag_dir, limiting_depths, det_filters, renamed_cat)
    #update_cat(snr_cols, 'cat_for_testing_on_SNR.fits', 'cat_for_testing_on_SNR.fits')



# old stuff
#cat_with_mags = 'det_H/COSMOS_DR2_MASKED_H_1.8as_IRAC2.8as.fits'




