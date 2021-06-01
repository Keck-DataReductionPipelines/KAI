from astropy.io import fits, ascii
from astropy.table import Table
import numpy as np
from astropy.stats import sigma_clipped_stats
import photutils.detection as sf
import microlens.jlu.align_flystar as af
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import os
from datetime import date


def euclidean(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def removeEdges(imgArr, finalTable, cutoff = 100, x = "xcentroid", y = "ycentroid", scale = 0.009942):
    """Ignores stars close to the edge, i.e., within 100 * scale pixels by default."""
    cutoff = (cutoff * scale) / 0.009942
    finalTable.remove_rows([i for i in range(len(finalTable[x])) if finalTable[x][i] < cutoff])
    finalTable.remove_rows([i for i in range(len(finalTable[y])) if finalTable[y][i] < cutoff])
    finalTable.remove_rows([i for i in range(len(finalTable[x])) if finalTable[x][i] > len(imgArr[0]) - cutoff])
    finalTable.remove_rows([i for i in range(len(finalTable[y])) if finalTable[y][i] > len(imgArr) - cutoff])
    return finalTable




def psfRadius(starPeak, imagePeak, scale = 0.009942):
    """
    Arbitrary function that ties the predicted radius of the star to the peak brightness in the image.
    returns a radius: 10 < (Ratio of "magnitudes" cubed * 80) < 60 and scaled acc. to plate scale in arcsec
    """
    return ((max(min(80 * ((np.log10(starPeak) / np.log10(imagePeak)) ** 3), 60), 10) / scale) * 0.00942)




def isBrightSpeckleStar(starPeak, imagePeak):
    """ Can this star have halo speckles? If it's half the peak brightness, it is assumed yes. """
    return (starPeak > (0.5 * imagePeak))



def isNearBrightStar(nearStarPeak, mainPeak):
    """ Is the main star near another bright star? """
    return (0.2 * mainPeak < nearStarPeak)



def isInPSFRadius(starCoo, mainCoo, psfr):
    """ Is star within the psf radius from earlier? """
    return (euclidean(starCoo[0], starCoo[1], mainCoo[0], mainCoo[1]) < psfr)



def isInSecondaryRadius(starCoo, mainCoo, psfr):
    """ Is star within the secondary radius, defined as 1.35x the psf radius? (arbitrarily) """
    return (euclidean(starCoo[0], starCoo[1], mainCoo[0], mainCoo[1]) < psfr * 1.35)



def getTableIndex(mainIndex, starIndex):
    """ Helper function to get star being considered in removeCloseStars. """
    return mainIndex + starIndex + 1


def removeCloseStars(imgArr, newTable, pixelDist = 10, x = "xcentroid", y = "ycentroid", verbose = False, scale = 0.009942):
    """ 
    In order of brightest stars in table, checks if it's a valid star, and tracks secondaries near it.
    Removes all invalid stars from table, i.e, those that are likely halo speckles, or too close for valid psf generation.   
    Edge stars are already removed from table using the removeEdges() method.
    
    Terminology?
    main - Star that is being considered for psf-ness. 
    star - Star that is being checked for nearness to main
    psf radius - Radius (arbitrary function tied to brightness) that is checked for halo speckles/close bright stars
    secondary radius - Radius in which, if star found, mark as secondary
    
    Pseudocode:
    
    Sort the table by brightest
    imagePeak = brightest star peak value in image
    invalidStarIndices = indices of stars to be removed from table
    
    for each star in table, that is not marked invalid yet:
        set this star to the main star
        find the psf radius of this star
        
        for all star-lookalikes less bright than main star:
        
            if star-lookalike is within psf radius (arbitrary function tied to brightness):
            
                if main star can have halo speckles (arbitrary brightness limit tied to imagePeak):
                
                    if lookalike is too bright to be a halo speckle (arbitrary brightness limit tied to mainPeak):
                        add main star to invalidIndices (another bright star is too near it!)
                        add bright star to invalidIndices (too close to main star)
                        
                    else lookalike is probably halo speckle:
                        add lookalike to invalidIndices (too close to main star, not actually a star)
                        keep the main star!!
                        
                else main star cannot have halo speckles:
                
                    if there is a BRIGHT star-lookalike nearby:
                        its definitely a star and not halo speckle AND star is too close!
                        add main star to invalidIndices (another bright star is too near it!)
                        add bright star to invalidIndices (too close to main star)
                        
            else if star-lookalike is within secondary radius (arbitrary function tied to brightness):
                set the "psf?" value of this star in table to 0 (secondary)
                
    from the table, remove all the invalid stars!
    return table
    
    """
    
    newTable.sort("peak", reverse = True)
    imagePeak = newTable["peak"][0]
    
    mainIndex = 0
    invalidStarIndices = []
    
    while mainIndex < len(newTable):
        if mainIndex not in invalidStarIndices:
            
            mainCoo = (newTable[x][mainIndex], newTable[y][mainIndex], newTable["peak"][mainIndex])
            psfr = psfRadius(mainCoo[2], imagePeak, scale)
            
            if verbose: print("Found PSF Star: " + str(np.round(mainCoo, 2)))
                
            for starIndex, starCoo in enumerate(zip(newTable[x][mainIndex + 1:], newTable[y][mainIndex + 1:], newTable["peak"][mainIndex + 1:])):
                if isInPSFRadius(starCoo, mainCoo, psfr):
                    if isBrightSpeckleStar(mainCoo[2], imagePeak):
                        if isNearBrightStar(starCoo[2], mainCoo[2]):
                            if verbose: print("Ignoring PSF Star Due To Bright Nearby Star: " + str(np.round(starCoo, 2)))
                            invalidStarIndices.append(mainIndex)
                        if verbose: print("Ignoring Halo Speckle/Bright Nearby Star: " + str(np.round(starCoo, 2)))
                        invalidStarIndices.append(getTableIndex(mainIndex, starIndex))
                    else:
                        if isNearBrightStar(starCoo[2], mainCoo[2]):
                            if verbose: print("Ignoring PSF Star Due To Bright Nearby Star: " + str(np.round(starCoo, 2)))
                            invalidStarIndices.append(mainIndex)
                            if verbose: print("Ignoring Nearby Star: " + str(np.round(starCoo, 2)))
                            invalidStarIndices.append(getTableIndex(mainIndex, starIndex))
                elif isInSecondaryRadius(starCoo, mainCoo, psfr):
                    if verbose: print("Adding Possible Secondary Star: " + str(np.round(starCoo, 2)))
                    newTable["psf"][getTableIndex(mainIndex, starIndex)] = 0
                    newTable["mx"][getTableIndex(mainIndex, starIndex)] = mainCoo[0]
                    newTable["my"][getTableIndex(mainIndex, starIndex)] = mainCoo[1]             
            if verbose: print("*****")
        mainIndex += 1
        
    newTable.remove_rows(invalidStarIndices)
    return newTable
    



def getStats(imgArr):
    """ Returns mean median and mode of pixel values in image array. """
    mean, median, std = sigma_clipped_stats(imgArr, sigma=3.0)
    return mean, median, std



def getNStarsHelper(imgArr, expStars = 10, counter = 0, starsFound = 0, fwhm = 4.5, scale = 0.009942):
    """ 
    Attempts to get N stars. If stars aren't found, returns what it can and presents the error. 
    If so, starfinding may have to be done manually.
    """
    
    std = getStats(imgArr)[2]
    med = getStats(imgArr)[1]
    thr = 0
    finderClass = sf.DAOStarFinder(fwhm = fwhm, 
                                   threshold = max(thr, std * 5), 
                                   sharphi = 0.4,
                                   roundlo = -0.2, 
                                   exclude_border = True)
    newTable = finderClass.find_stars(imgArr - med)
    newTable["psf"] = 1
    newTable["mx"] = 0.0
    newTable["my"] = 0.0
    if len(newTable["peak"]) < expStars:
        print("Pre-Filtering PSF Stars Found: " + str(len(newTable["peak"])))
        print("More Stars Required. Widening Search Parameters.")
        nowStarsFound = len(newTable["peak"])
        newTable = getNStarsHelper(imgArr, expStars, counter + 1, nowStarsFound, fwhm)
    newTable = removeEdges(imgArr, newTable, scale = scale)
    newTable = removeCloseStars(imgArr, newTable, verbose = False, scale = scale)
    if ((len(newTable["peak"]) < expStars) and (counter < 20)):
        print("PSF Stars Found: " + str(len(newTable["peak"])))
        print("More Stars Required. Widening Search Parameters.")
        nowStarsFound = len(newTable["peak"])
        if (nowStarsFound < starsFound) or (std * 5 > thr):
            print("Failed to find sufficient stars. Do starfinding manually or reduce PSF stars to be found.".upper())
            return newTable
        newTable = getNStarsHelper(imgArr, expStars, counter + 1, nowStarsFound, fwhm)
    if counter >= 20:
        print("PSF Stars Found: " + str(len(newTable["peak"])))
        print("Failed to find sufficient stars. Do starfinding manually or reduce PSF stars to be found.".upper())
        return newTable
    return newTable


def getNStars(imgArr, numStars = 10, fwhm = 4.5, scale = 0.009942):
    """ 
    Reformats table of chosen psf stars and adds the relevant secondaries to it.
    """
    newTable = getNStarsHelper(imgArr, numStars, fwhm = fwhm, scale = scale)
    print("PSF Stars Found: " + str(len(newTable["peak"])) + "/" + str(numStars))
    print("Returning Best: " + str(min(len(newTable["peak"]), numStars)))
    newTable.sort("peak", reverse = True)
    secondaryTable = newTable[:]
    nonSecIdx = []
    for psfIdx in range(len(secondaryTable)):
        if (secondaryTable["psf"][psfIdx] != 0):
            nonSecIdx.append(psfIdx)
    secondaryTable.remove_rows(nonSecIdx)
    removed = []
    for i in range(len(newTable)):
        if newTable[i]["psf"] == 0:
            removed.append(i)
    newTable.remove_rows(removed)
    newTable.remove_rows(range(numStars, len(newTable)))
    goodSecIdx = []
    for i in range(len(secondaryTable)):
        for j in range(len(newTable)):
            if (secondaryTable["mx"][i] == newTable["xcentroid"][j]) and (secondaryTable["my"][i] == newTable["ycentroid"][j]):
                goodSecIdx.append(i)
                break
    for i in goodSecIdx:
        newTable.add_row(secondaryTable[i])
    return newTable

def generate_list(imgPath, numPSFStars = 7, fwhm = 4.5, scale = 0.009942, targName = ""):
    """
    Plots the table generated by PSFListGen, and prints the list of stars in a format that can be used.
    @params
    imgPath - path to .fits image, Eg. /g/lu/data/microlens/18jun22/combo/mag18jun22_ob040361_kp.fits
    numPSFStars - number of stars to be found
    scale - defaults to nirc2 scale. Enter in arcsec/pixel.
    targName - name of target (can be anything, not important)
    fwhm - full-width half maximum variable for finicky images
    """
    
    imgArr = fits.getdata(imgPath)
    newTable = getNStars(imgArr, numStars = numPSFStars, fwhm = fwhm, scale = scale)
    finalTable = newTable["xcentroid", "ycentroid", "peak", "psf", "mx", "my"]
    finalTable.rename_column("xcentroid", "x")
    finalTable.rename_column("ycentroid", "y")
    finalTable.rename_column("psf", "m")
    finalTable["x"] = np.round(finalTable["x"], 2)
    finalTable["y"] = np.round(finalTable["y"], 2)
    finalTable["peak"] = np.round(finalTable["peak"], 2)
    finalTable["mx"] = np.round(finalTable["mx"], 2)
    finalTable["my"] = np.round(finalTable["my"], 2)
    finalTable["name"] = ["loremipsum" for i in finalTable["x"]]
    finalTable["name"][0] = targName
    modded_plot_starlist_on_image_pixel(finalTable, imgPath, targName,
                                    flip = False, label = False, magCut = 400000, verbose = False)
    print("*************************************************")
    print("PSF List Details: ")
    finalTable.remove_column("name")
    finalTable.pprint()
    print("*************************************************")
    print("PSFListGen Identified the Following List of Stars:")
    print("[")
    for i in zip(finalTable["x"], finalTable["y"], finalTable["m"]):
        print("   " + str(list(i)) + ",")
    print("]")
    print("*************************************************")
    

def plot_starlist_on_image_arcsec(starList, imagePath, refCoo, scale = (9.942/1000), magCut = 23, label = True, verbose=True, flip = False):
    """
    Plot a NIRC2 image and overlay a starlist. Input in relative arcsec coordinates.
    @params
    starList - any astropy Table with a "name" (including target),
                                        "m", 
                                        "x" (relative arcsec coordinates), and 
                                        "y" (relative arcsec coordinates) column. 
    imagePath - path to image from root.
    refCoo - reference pixel coordinates of a star in image file.
    scale - plate scale. Defaults to NIRC2 plate scale.
    magCut - highest magnitude of star allowed to be plotted.
    label - if True, adds labels. Defaults to true.
    verbose - if True, prints the astropy table inputted. Defaults to true.
    flip - if True, flips image in x (helps in case image is backwards).
    
    Ensure that the image passed in is at a Position Angle of 0.
    
    The image is NOT Flipped on plotting by default. 
    In order to flip image, set flip = True.
    """
    
    
    # Initializing the Image and starList columns 
    img = fits.getdata(imagePath)
    
    try:
        xCoordList = starList['x']
        yCoordList = starList['y']
        mList = starList['m']
    except KeyError:
        raise KeyError("Starlist must have columns named 'x' and 'y', in pixel coordinates, and magnitude 'm'.")
        
    if label:
        try:
            nameList = starList['name']
        except KeyError:
            raise KeyError("Starlist must have a column named 'name' if label = True. Else, set label = False.")
    
   
    # Get image dimensions
    x_axis = np.arange(img.shape[0], dtype=float) #[0, 1, ...., 1166]
    y_axis = np.arange(img.shape[1], dtype=float)
    x_axis = ((x_axis + 1 - refCoo[0]) * scale)
    y_axis = ((y_axis + 1 - refCoo[1]) * scale)
    
    extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
    
    if flip:
        img = np.flip(img, axis = 1)
        x_axis *= -1
        extent = [x_axis[-1], x_axis[0], y_axis[0], y_axis[-1]]
    
    # Plotting prerequisites
    vmin = 10
    vmax = 1e5
    norm = LogNorm(vmin, vmax)
    
    if verbose:
        print(starList)
    
    # Plotting the image
    idx = np.where(mList < magCut)[0]
    plt.figure(figsize=(10,8))

    plt.imshow(img, cmap='gist_heat_r', norm=norm, extent=extent, origin = "lower")
    

    # Scatter the sources in starlist.
    sc = plt.scatter(xCoordList[idx], yCoordList[idx], c = mList[idx], s = 300, linewidth = 3)
    sc.set_facecolor('none')
    
    # Label the sources.
    if label:
        for i, txt in enumerate(nameList[idx]):
            plt.annotate(txt, (xCoordList[i], yCoordList[i] + 0.3))
    
    # Add titles, etc.
    plt.colorbar(label='Starlist Magnitude (mag)')
    plt.xlabel('Distance from Reference Star (arcsec)')
    plt.ylabel('Distance from Reference Star (arcsec)')
    plt.axis('equal')
    if label:
        targIndex = np.where(xCoordList == 0.0)[0]
        plt.title(nameList[targIndex][0].upper())
    else:
        plt.title("Target at " + str(refCoo))
    
    return 

def plot_starlist_on_image_pixel(starList, imagePath, refStar, scale = (9.942/1000), magCut = 23, label = True, verbose=True, flip = False, axes = "pixel"):
    """
    Plot a NIRC2 image and overlay a starlist. Input in pixel coordinates.
    @params
    starList - any astropy Table with a "name" (including target),
                                        "m", 
                                        "x" (pixel coordinates), and 
                                        "y" (pixel coordinates) column. 
    imagePath - path to image from root.
    refStar - name of reference star in astropy table, generally the target.
    scale - plate scale. Defaults to NIRC2 plate scale.
    magCut - highest magnitude of star allowed to be plotted.
    label - if True, adds labels. Defaults to true.
    verbose - if True, prints the astropy table inputted. Defaults to true.
    flip - if True, flips image in x (helps in case image is backwards).
    axes - determines units on the axes of plot (pixel or arcsec). Defaults to pixel.
    
    Ensure that the image passed in is at a Position Angle of 0.
    
    The image is NOT Flipped on plotting by default. 
    In order to flip image, set flip = True.
    
    """
    
    # Initializing the Image and starList columns 
    img = fits.getdata(imagePath)
    
    try:
        xCoordList = starList['x']
        yCoordList = starList['y']
        mList = starList['m']
        nameList = starList['name']
    except KeyError:
        raise KeyError("Starlist must have columns named 'name', 'x' and 'y', in pixel coordinates, and magnitude 'm'.")
    
    # Find the target and get its coordinates.
    targetIndex = np.where(nameList == refStar)[0]
    if len(targetIndex) == 0:
        raise IndexError("Target does not match any name in starlist.")
    
    # Get image dimensions and make relative to reference
    x_axis = np.arange(img.shape[0], dtype=float)
    y_axis = np.arange(img.shape[1], dtype=float)
    x_axis = (x_axis -  (xCoordList[targetIndex] - 1))
    y_axis = (y_axis - (yCoordList[targetIndex] - 1))

    
    # Make starList coordinates relative to reference
    xCoordList = (xCoordList - xCoordList[targetIndex])
    yCoordList = (yCoordList - yCoordList[targetIndex])
    
    # Handle if axes are asked in asrcsec
    if axes == "arcsec":
        x_axis *= scale
        y_axis *= scale
        xCoordList *= scale
        yCoordList *= scale
    
    # Extent of image to beplotted in imshow
    extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
    
    # Flips image in case it's backwards
    if flip:
        x_axis *= -1
        img = np.flip(img, axis = 1)
        extent = [x_axis[-1], x_axis[0], y_axis[0], y_axis[-1]]
    
    # Plotting prerequisites
    vmin = 10
    vmax = 1e5
    norm = LogNorm(vmin, vmax)
    
    # Verbose option
    if verbose:
        print(starList)
        
    
    # Plot the image
    idx = np.where(mList < magCut)[0]
    plt.figure(figsize=(10,8))

    plt.imshow(img, cmap='gist_heat_r', norm=norm, extent=extent, origin = "lower")
    

    # Scatter the sources in starlist.
    sc = plt.scatter(xCoordList[idx], yCoordList[idx], c = mList[idx], s = 300, linewidth = 3)
    sc.set_facecolor('none')
    
    # Label the sources.
    if label:
        for i, txt in enumerate(nameList[idx]):
            if axes == "pixel":
                plt.annotate(txt, (xCoordList[i], yCoordList[i] + (0.3 / scale)))
            elif axes == "arcsec":
                plt.annotate(txt, (xCoordList[i], yCoordList[i] + 0.3))
            else:
                raise ValueError("'axes' variable must be 'pixel' or 'arcsec'")
    
    # Plot titles, etc.
    plt.colorbar(label='Starlist Magnitude (mag)')
    plt.xlabel('Distance from Reference Star (' + axes + ")")
    plt.ylabel('Distance from Reference Star (' + axes + ")")
    plt.axis('equal')
    plt.title(refStar.upper())
    
    
    return


def modded_plot_starlist_on_image_pixel(starList, imagePath, refStar, scale = (9.942/1000), magCut = 23, label = True, verbose=True, flip = False, axes = "pixel"):
    """
    Plot a NIRC2 image and overlay a starlist. Input in pixel coordinates.
    @params
    starList - any astropy Table with a "name" (including target),
                                        "m", 
                                        "x" (pixel coordinates), and 
                                        "y" (pixel coordinates) column. 
    imagePath - path to image from root.
    refStar - name of reference star in astropy table, generally the target.
    scale - plate scale. Defaults to NIRC2 plate scale.
    magCut - highest magnitude of star allowed to be plotted.
    label - if True, adds labels. Defaults to true.
    verbose - if True, prints the astropy table inputted. Defaults to true.
    flip - if True, flips image in x (helps in case image is backwards).
    axes - determines units on the axes of plot (pixel or arcsec). Defaults to pixel.
    
    Ensure that the image passed in is at a Position Angle of 0.
    
    The image is NOT Flipped on plotting by default. 
    In order to flip image, set flip = True.
    
    """
    
    # Initializing the Image and starList columns 
    img = fits.getdata(imagePath)
    
    try:
        xCoordList = starList['x']
        yCoordList = starList['y']
        mList = starList['m']
        nameList = starList['name']
    except KeyError:
        raise KeyError("Starlist must have columns named 'name', 'x' and 'y', in pixel coordinates, and magnitude 'm'.")
    
    # Find the target and get its coordinates.
    targetIndex = np.where(nameList == refStar)[0]
    if len(targetIndex) == 0:
        raise IndexError("Target does not match any name in starlist.")
    
    # Get image dimensions and make relative to reference
    x_axis = np.arange(img.shape[0], dtype=float)
    y_axis = np.arange(img.shape[1], dtype=float)
    x_axis = (x_axis)
    y_axis = (y_axis)

    
    # Make starList coordinates relative to reference
    xCoordList = (xCoordList)
    yCoordList = (yCoordList)
    
    # Handle if axes are asked in arcsec
    if axes == "arcsec":
        x_axis *= scale
        y_axis *= scale
        xCoordList *= scale
        yCoordList *= scale
    
    # Extent of image to be plotted in imshow
    extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
    
    # Flips image in case it's backwards
    if flip:
        x_axis *= -1
        img = np.flip(img, axis = 1)
        extent = [x_axis[-1], x_axis[0], y_axis[0], y_axis[-1]]
    
    # Plotting prerequisites
    vmin = 10
    vmax = 1e5
    norm = LogNorm(vmin, vmax)
    
    # Verbose option
    if verbose:
        print(starList)
        
    
    # Plot the image
    idx = np.where(mList < magCut)[0]
    plt.figure(figsize=(10,8))

    plt.imshow(img, cmap='gist_heat_r', norm=norm, extent=extent, origin = "lower")
    plt.colorbar(label = "Raw Pixel Value")
    # Scatter the sources in starlist.
    sc = plt.scatter(xCoordList[idx], yCoordList[idx], color = list(zip(1.0 - mList[idx], [0.0 for i in mList[idx]], mList[idx])) , s = 300, linewidth = 3)
    sc.set_facecolor('none')
    
    # Label the sources.
    if label:
        for i, txt in enumerate(nameList[idx]):
            if axes == "pixel":
                plt.annotate(txt, (xCoordList[i], yCoordList[i] + (0.3 / scale)))
            elif axes == "arcsec":
                plt.annotate(txt, (xCoordList[i], yCoordList[i] + 0.3))
            else:
                raise ValueError("'axes' variable must be 'pixel' or 'arcsec'")
    
    # Plot titles, etc.
    plt.xlabel('Physical Image X (' + axes + ")")
    plt.ylabel('Physical Image Y (' + axes + ")")
    plt.axis('equal')
    plt.title(refStar.upper())
    
    
    return

def prepStarfinder(dirName, targetName, targetCoords, psfStars, filterName):
        """Creates a _psf.list file and saves it in the source_list/ directory."""
        
        # Covering possible variable inputs and initializing
        targetName = targetName.lower()
        filterName = filterName.capitalize()
        targetCoords = targetCoords[:]
        psfListLocation = dirName + targetName + '_psf.list'
        year = date.today().strftime("%Y" + ".0")
        
        # Modifying the starlist provided into an astropy table
        starno = 0
        for star in psfStars:
            star[0] = "{:.3f}".format(((star[0] - targetCoords[0]) * (9.942 / 1000) * -1))
            star[1] = "{:.3f}".format(((star[1] - targetCoords[1]) * (9.942 / 1000)))
            star.insert(0, "1.000")
            if starno > 0:
                star.insert(0, "S" + str(starno).zfill(3))
            else:
                star.insert(0, targetName)
            star.insert(4, "-")
            star.insert(4, year)
            star.insert(4, "0.000")
            star.insert(4, "0.000")
            starno += 1
            
        psfStars = np.array(psfStars)
        psfStarTable = Table(psfStars, names = ('#Name', 'Kp', 'Xarc',
                                                      'Yarc', 'Vx', 'Vy',
                                                      't0', 'Filt', 'PSF?'))
        
        ascii.write(psfStarTable, psfListLocation, 
                    format = 'fixed_width', 
                    delimiter = '  ', 
                    bookend = False, 
                    overwrite = True)