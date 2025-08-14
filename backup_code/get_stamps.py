#!/usr/bin/env python3

""" get_stamps.py

Aim: Cut out stamps of the objects of interest.

Created: Tuesday 8th May 2018
To python3: 25th March 2020
Modified: Wed 5th March 2025

ds9 -mecube...
"""

###################### Import useful modules #########################
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import astropy.units as u

import sys
# insert at 1, 0 is the script path (or '' in REPL)
#sys.path.insert(1, '/mnt/vardy/vardygroupshare/olddata/codes/')
from new_catalogue_codes import which_field, cut_out

############################ Set-up code  ############################
# Edit this
# This catalogue must have ID, RA and DEC in it!
field = 'COSMOS' # or 'COSMOS'
#field = 'CDFS' # XMM

#catDir 
catDir = '/raid/scratch/analysis/previous_cats/'
catDir = '/raid/scratch/bowlerr/mphys/'

inputCat = 'Bowler2014_COSMOS.fits'
inputCat = 'close_pair_unique_z4_2to3arcsec.fits'
inputCat = 'missing_sources.csv'

requiredFilters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y', 'CFHT-g', 'CFHT-r', 'CFHT-iy', 'CFHT-z', 'Y', 'J', 'H', 'Ks']

requiredFilters = ['u', 'g', 'r', 'i', 'HSC-G','HSC-R', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
requiredFilters = ['HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'K', 'Ye', 'Je', 'He']

# can set to ['all'] for all available filters
#requiredFilters = [''Y']
stampSize = 15.0 ## arcseconds across
overwrite = True # overwrite stamps of the same name

#######################################################################
## Leave this
interactive = False
outputDir = '/raid/scratch/temporaryfiles/stamps/' 

######################### Read in catalogue ##########################
## want ra and dec and id
if inputCat[-4:] == 'fits':
   objectdata = Table.read(catDir + inputCat)
elif inputCat[-3:] == 'csv':
#   objectdata = np.genfromtxt(catDir + inputCat, names=True, dtype = None, \
#                              delimiter = ',')
   objectdata = Table.read(catDir + inputCat)
   
else:
   ## assume ascii format
   #objectdata = np.genfromtxt(catDir + inputCat, names=True, dtype = None)
   objectdata = Table.read(catDir + inputCat, format = 'ascii.commented_header')

   # I need to fix the ra/dec
   c = SkyCoord(ra=objectdata['RAs'], dec =objectdata['DECs'], unit=(u.hourangle, u.deg))
   radeg = c.ra.degree
   decdeg = c.dec.degree
   newcol = Column(radeg, name = 'RA')
   objectdata.add_column(newcol)
   newcol = Column(decdeg, name = 'DEC')
   objectdata.add_column(newcol)
   
#
print(objectdata)
indices = np.arange(0, len(objectdata), 100)

objectdata = objectdata[indices]

print("There are {0} objects to cut out.".format(len(objectdata)))
print(objectdata.colnames)
#objectdata.remove_column('RA')
#objectdata.remove_column('DEC')
#objectdata.rename_column('ALPHA_J2000', 'RA')
#objectdata.rename_column('DELTA_J2000', 'DEC')

#ii = objectdata['RA'] > 100
#objectdata = objectdata[ii]
#print(objectdata)

ra = objectdata['RA']
dec = objectdata['DEC']
#ra = objectdata['UV_ra']
#dec = objectdata['UV_dec']
IDhere = objectdata['ID']

## and extract basename for directory name
ee = inputCat.rfind('.')
basename = inputCat[:ee]

if outputDir != 'here':
   if os.path.isdir(outputDir) == False:
      os.system('mkdir ' + outputDir)
   stampDir = outputDir + basename
else:
   stampDir = basename
    
if os.path.isdir(stampDir) == False:
   os.system('mkdir ' + stampDir)
print("Stamps will be saved to ", stampDir) 

imageDir = '/raid/scratch/data/'

########################### Check catalogue ##########################
## Organise into the different fields.
fields = which_field(ra, dec)

## get the unique fields
uniqueFields = np.unique(fields)
print("The unique fields are ", uniqueFields)

## Work out which fields we are in!
for ff, fieldName in enumerate(uniqueFields):

   ## reduce the ra and dec and id down
   hh = (fields == fieldName)
   raField = ra[hh]
   decField = dec[hh]
   idField = IDhere[hh]
      
   ## Extract the information here
   dirHere = imageDir + fieldName + '/'
   #print("The directory is ", dirHere)
   
   ## Find the required images and their locations
   ## Read in the images.lis file here
   inputFile = dirHere + 'images.lis'
   if os.path.isfile(inputFile):
      print("Reading in images...")
   else:
      print("No ", inputFile, " file exists!  Exiting...")
      exit()
        
   imagedata = Table.read(inputFile, format = 'ascii.commented_header')
   availableFilters = imagedata['Name']
   print("The available filters are ", availableFilters)
   print("The requiredFilters ", requiredFilters)
   
   if requiredFilters[0] == 'all':
      requiredFilters = availableFilters
   
   ## loop through the required filters!
   ## first check we can actually find all of them here
   for rf, reqFilt in enumerate(requiredFilters):
      ii = (reqFilt == availableFilters)
      
      if np.sum(ii) < 1:
         print("Missing filter ", reqFilt)
#         print("Available filters are ", availableFilters)
      else:
         print("Found filter ", reqFilt, imagedata['Image'][ii][0])
   
   if interactive:
      input("Check filters and press enter to continue... ")
      
   ## Now loop through the images again
   for rf, reqFilt in enumerate(requiredFilters):
      
      ii = (reqFilt == availableFilters)

      if np.any(ii) == False:
         print("ERROR cannot find this filter in the available filters.", reqFilt)
      else:
         
         #print imagedata['directory'][ii]
         #print imagedata['Image'][ii]
         #print imagedata['directory'][ii]
         #print imagedata['directory'][ii][0] == "here"
         
         ## Cut these out!
         if imagedata['directory'][ii][0] == "here":
            dataDir = dirHere
         else:
            dataDir = imagedata['directory'][ii] 
            #print("data dir = ", dataDir)
            dataDir = dataDir[0]
         
         imageName = imagedata['Image'][ii]
         print("Image name ", imageName[0])
         
         fitsMEF = cut_out(dataDir + imageName[0], raField, decField, idField, filterName = fieldName + '_' + reqFilt, stampSizeAS =  stampSize, header = False, overwrite = True, outDir = stampDir)
      

   print("Now extracting what I need!")
   ## Sort out all the different results
   for oi, obj in enumerate(idField):

      # first check if this exists yet
      fitsName = stampDir + '/ID{0}.fits'.format(obj)
#      print(os.path.isfile(fitsName), overwrite)
      if (os.path.isfile(fitsName) == False) or overwrite:

         print('Making mef')
         ## Make a MEF for this object
         new_hdul = fits.HDUList()
         new_hdul.append(fits.PrimaryHDU())
         ## And through the filters
         for rf, reqFilt in enumerate(requiredFilters):
            
            ## extract the data for this object, and append
            imageName = stampDir + '{0}_{1}_tmp.fits'.format(fieldName, reqFilt)
            print(os.path.isfile(imageName), imageName)
 #           exit()
            if os.path.isfile(imageName):
               hdu_image = fits.open(imageName, memmap=True)
               #hdu = hdu_image[0]
               #      imageData = hdu_image['{0}'.format(obj)].data
               #print("Slow", imageData[0:5, 0:5])
               imageData = hdu_image[oi+1].data
               #               print("Fast", imageData[0:5,0:5])

               hdrhere = hdu_image[oi+1].header
               
               nnhdu = fits.ImageHDU(imageData, header = hdrhere, name = reqFilt.strip())
               nnhdu.header.set('EXTNAME', reqFilt.strip())
               new_hdul.append(nnhdu)
               
               hdu_image.close()
            
         new_hdul.writeto(fitsName, overwrite = True)
         #         print("Fits saved to ", fitsName)
         new_hdul.close()
         
         if (oi % 100) == True:
            print("At object {0}/{1}".format(oi, len(idField)))
   ## tidy up
#   new_hdul.close()
os.system('rm '+stampDir +'*_tmp.fits')
    
print("Finished cutting out!")
print(stampDir)
