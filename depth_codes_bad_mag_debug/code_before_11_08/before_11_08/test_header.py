
def test_imgHeader(imageName):
    from datetime import datetime

    print("Script run at:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    from astropy.io import fits
    import os
    # get the pixel scale
    hdulist = fits.open(imageName)
    imageHeader = hdulist[0].header
    imageHeaderKeys = list(imageHeader.keys())

    # ensure header of data file is as expected 
    # i.e. if using a copy/cutout, header should == header of original file

    imgBase = os.path.basename(imageName)
    if 'DR6_' in imgBase:
        # the original image to which we will compare
        originalImageName = imgBase.split('DR6_')[0] + 'DR6_backup.fits'
        originalHDUlist = fits.open('/raid/scratch/data/COSMOS/'+originalImageName)
        originalHeader = originalHDUlist[0].header
        originalImageHeaderKeys = list(originalHeader.keys())

        if imageHeaderKeys == originalImageHeaderKeys:
            print("ImageHeader matches expected format. Checking values...")
            
            for key in originalImageHeaderKeys:
                if originalHeader[key] == imageHeader[key]:
                    print(f"original == image value: {originalHeader[key] == imageHeader[key]}", key, originalHeader[key], imageHeader[key])
                else:
                    print("original == image value: {originalHeader[key] == imageHeader[key]}", key, originalHeader[key], imageHeader[key])
    
        elif imageHeaderKeys != originalImageHeaderKeys:

            imageHeaderKeys = set(imageHeaderKeys)
            originalImageHeaderKeys = set(originalImageHeaderKeys)
            uniqueImageHeaderKeys = sorted(imageHeaderKeys - originalImageHeaderKeys)
            uniqueOrigImageHeaderKeys = sorted(originalImageHeaderKeys - imageHeaderKeys)

            warnings.warn(">>>>>>> WARNING: imageHeader does not match expected Header. Known issues: incorrect calculation of pixScale. \n")
            print("Keys only in imageHeader: ",uniqueImageHeaderKeys)
            print("Keys only in originalImageHeader: ",uniqueOrigImageHeaderKeys)
    else:
        cd1Test = ['CD1_1' in imageHeaderKeys]
        print(f">>>> WARNING: Passing file to SE. Header includes CD1_1?: {cd1Test} Is this the right file?", '\n', imageName, '\n') 
        cdone_o = -3600.0*imageHeader['CD1_1'] # TODO - this might break - not tested yet

    if 'CD1_1' in imageHeader:
        cdone_o = -3600.0*imageHeader['CD1_1']

    return imageHeader


a = test_imgHeader('/raid/scratch/data/data_test/COSMOS_test/UVISTA_J_DR6_backup.fits')

b = test_imgHeader('/raid/scratch/data/COSMOS/UVISTA_JH_DR6_backup.fits')
import difflib

a_header = test_imgHeader('/raid/scratch/data/data_test/COSMOS_test/UVISTA_J_DR6_backup.fits')
b_header = test_imgHeader('/raid/scratch/data/COSMOS/UVISTA_JH_DR6_backup.fits')

a_str = repr(a_header)  # or str(a_header)
b_str = repr(b_header)

if a_str == b_str:
    print("Headers are identical.")
else:
    print("Headers differ:")
    #for line in difflib.ndiff(a_str.splitlines(), b_str.splitlines()):
        #print(line)
        



