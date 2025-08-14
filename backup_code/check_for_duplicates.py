

def exclude_filters(cat_dir_list):
    cat_list = []
    for cat_dir in cat_dir_list:
        if '_Y_' not in cat_dir and 'BACKUP_' not in cat_dir:
            cat_list.append(cat_dir)
            print(f"Catalogue found: {cat_dir} \n")
    return cat_list

def check_duplicates(cat_dir_list, quick_test=True, cat_dir=None):
    """Checks catalogues for duplicated ID numbers."""
    if cat_dir==None:    
        raise ValueError(">>>>>>>>>>>>>>>>>>>>>>> ERROR: cat_dir not supplied")
    from collections import defaultdict
    all_ids = []
    all_radec = defaultdict(list)  # Temporarily stores RA/DEC for all IDs
    id_to_file = defaultdict(set)  # Maps each ID to the files it appears in

    print("Opening catalogues...\n")
    from astropy.io import fits
    for cat_dir in cat_dir_list:
        with fits.open(cat_dir) as hdul:
            data = hdul[1].data
            if quick_test:
                data = data[:20]  # Test on a short catalog (only x rows)
            for id_, ra, dec in zip(data['ID'], data['RA'], data['DEC']):
                all_radec[id_].append((ra, dec))  # Save every RA/DEC
                id_to_file[id_].add(cat_dir)      # Track which files the ID appears in
            all_ids.append(data['ID'])

    # Find duplicated IDs
    import numpy as np
    combined_cat = np.concatenate(all_ids)
    unique_vals, counts = np.unique(combined_cat, return_counts=True)
    dup_vals = unique_vals[counts > 1]
    dup_instances = counts[counts > 1]

    # Now create `id_to_radec` ONLY for duplicated IDs
    id_to_radec = {id_: all_radec[id_] for id_ in dup_vals}

    print(f"\nNumber of duplicate IDs: {len(dup_instances)}")
    if len(dup_instances) > 0:
        import warnings
        warnings.warn(
            f"\nDuplicated IDs: {dup_vals} repeated\n"
            f"{dup_instances} times, respectively \n"
        )
        for dup_id in dup_vals[:10]:  # Show first 10 examples
            print(f"ID: {dup_id} found in:\n  {'... '.join(id_to_file[dup_id])}\n")

    return id_to_file, id_to_radec, dup_vals, dup_instances

def check_for_duplicate_IDs(cats_to_check, cat_dir=None):
    """Runs the functions in the check_for_duplicates programme and returns number of duplicates found."""
    cat_list = exclude_filters(cats_to_check)
    id_to_file, id_to_radec, dup_vals, dup_instances = check_duplicates(cat_list, cat_dir=cat_dir)
    return dup_instances



