
def _run_one(filterName, detFilt, imagedata, apDiameterAS, IRACapDiametersAS, 
             inputSex, dirHere, queue, field, indiDir, overwrite, memory):


    print(f"Starting {filterName}, {detFilt}")

    if filterName[0:2] == 'ch':
        apDiameterASuse = IRACapDiametersAS
        print(f'Changing apertures for filter {filterName}, to {apDiameterASuse}')
    else:
        apDiameterASuse = apDiameterAS

    return run_se(
        detFilt, filterName, imagedata, apDiameterASuse, inputSex,
        dirHere, queue=queue, field=field, outputDir=indiDir,
        overwrite=overwrite, memory=memory)


def run_source_extractor(field, detFilt, apDiameterAS, queue = 'none', reqFilters = ['all'], excludeFilters = ['ch1o', 'ch2o', 'ugriy', 'GRI', 'GR'], interactive = True, overwrite = False, imageDir = baseDir + 'data/', memory = 10.0, IRACapDiametersAS = [2.8, 3.8, 5.8, 9.8, 11.6]):
    from concurrent.futures import ProcessPoolExecutor, as_completed

    print("Running run_source_extractor: ", '\n')
    # Note: removed 'YJ' from excludeFilters.
    # Note: removed 'JH', 'HKs' from excludeFilters.

    # standard inputs etc

    inputSex = baseDir + 'data/'+'bertin_config/video_mine.sex'
    #os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor/2.25.0/share/sextractor'
    os.environ['EXTRACTOR_DIR'] = '/usr/local/sextractor-2.25.0/config'
    # define the directories
    ## Extract the information here
    dirHere = imageDir + field + '/'
    print("The directory is ", dirHere)
    
    catDir = imageDir + 'catalogues/{0}/'.format(field)
    indiDir = catDir + 'det_{0}/'.format(detFilt)
    
    if os.path.isdir(catDir) and os.path.isdir(indiDir):
        print("Directories created.")
    else:
        os.system('mkdir ' + catDir)
        os.system('mkdir ' + indiDir)

    
    ###################################################
    ################# Sort filters out ################
    ## Find the required images and their locations
    ## Read in the images.lis file here
    inputFile = dirHere + 'images.lis'
    if os.path.isfile(inputFile):
        print("Reading in images...")
    else:
        print("No ", inputFile, " file exists!  Exiting...")
        exit()
        
    imagedata = Table.read(inputFile, format = 'ascii.commented_header')
    availableFilters = np.array(imagedata['Name'])
    print("The available filters are ", availableFilters)
    #for i, ii in enumerate(availableFilters):
    #    print(ii)
    # now remove the undesirable filters!
    
    excludeFilters = np.array(excludeFilters)
    
    if excludeFilters.size > 1:
        keep = np.ones(availableFilters.size, dtype=bool)
        for bf, badFilt in enumerate(excludeFilters):
            ii = np.where((availableFilters == badFilt) & (badFilt != detFilt))
            keep[ii] = False
            print("Excluding filter ", badFilt)
            
        imagedata = imagedata[keep]
        
    availableFilters = imagedata['Name']

    if reqFilters[0] == 'all':
        reqFilters = np.copy(availableFilters)

    else:
        print("Only running for a selection of filters ")
        availableFilters = np.array(reqFilters)

    # check!
    if interactive:
        print("I will now run SE for {0} filters in field {1}, on queue = {2}.".format(availableFilters.size, field, queue))
        print("The filters are ", availableFilters)
        #input("Press enter to continue.")
    
        
    # check the detection filter exists in catalogue
    detIndex = (availableFilters == detFilt)
    print('detIndex', detIndex)
    
    if np.sum(detIndex) < 1:
        print("The input detection filter does not exist.")
        print("Input = ", detFilt, " out of possible filters: ", availableFilters)
        
    else:
        detIndex = detIndex[0]
        
    # Change the deblending...

    #with ProcessPoolExecutor() as executor:
    from concurrent.futures import ThreadPoolExecutor, as_completed

    with ThreadPoolExecutor(max_workers=11) as executor: 
        futures = {
            executor.submit(
                _run_one, f, detFilt, imagedata, apDiameterAS,
                IRACapDiametersAS, inputSex, dirHere, queue,
                field, indiDir, overwrite, memory): f for f in availableFilters}

        for fut in as_completed(futures):
            filterName = futures[fut]
            try:
                fut.result()
                print(f"Finished {filterName}")
            except Exception as e:
                print(f"Error with filter {filterName}: {e}")

    return
