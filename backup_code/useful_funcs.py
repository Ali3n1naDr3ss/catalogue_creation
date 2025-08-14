def check_file_exists(catdir):
    """returns list of file directories that lead to a cgs catalogue"""
    import os    
    dirs = os.listdir(catdir)    
    subdirs = []
    file_list = []

    for name in dirs:
        subdir = os.path.join(catdir, name)
        if os.path.isdir(subdir):
            subdirs.append(subdir) # check if dir exists

    for subdir in subdirs:
        files = os.listdir(subdir)
        for fname in files:
            if 'MASKVISTA' in fname and 'Ks' in fname and fname.endswith('cgs.fits'):
                full_path = os.path.join(subdir, fname) # get cgs cats
                file_list.append(full_path)
            elif 'MASKVISTADET' in fname and fname.endswith('cgs.fits'):
                full_path = os.path.join(subdir, fname) # get cgs cats
                file_list.append(full_path)
    return file_list



def get_col_data(cat_dir, col_name):
    with fits.open(cat_dir) as hdul:
        data = hdul[1].data
        col_names = data.columns.names
        data = hdul[1].data[col_name]


def open_fits():
    for cat_dir in cat_dir_list:
        with fits.open(cat_dir) as hdul:
            data = hdul[1].data
            if quick_test:
                data = data[:20]  # Test on a short catalog (only x rows)
        orig_col_names = data.columns.names
        orig_cols = raw_data.columns
