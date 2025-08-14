



ra_corner = 150   # in degrees
dec_corner = 2    # in degrees
cutout_size = 500    # arcseconds

filter_ = 'JH'
image_path = '/raid/scratch/data/COSMOS/UVISTA_'+filter_
image_type = '_DR6_backup.fits'
make_cutout(image_path+image_type, ra_corner, dec_corner, cutout_size, corner='lower-right')

image_type = '_DR6_wht_backup.fits'
make_cutout(image_path+image_type, ra_corner, dec_corner, cutout_size, corner='lower-right')


filter_ = 'J'
image_path = '/raid/scratch/data/COSMOS/UVISTA_'+filter_
image_type = '_DR6_backup.fits'
make_cutout(image_path+image_type, ra_corner, dec_corner, cutout_size, corner='lower-right')

image_type = '_DR6_wht_backup.fits'
make_cutout(image_path+image_type, ra_corner, dec_corner, cutout_size, corner='lower-right')



