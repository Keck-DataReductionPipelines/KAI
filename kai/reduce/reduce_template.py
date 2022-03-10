##################################################
#
# General Notes:
# -- python uses spaces to figure out the beginnings
#    and ends of functions/loops/etc. So make sure
#    to preserve spacings properly (indent). This
#    is easy to do if you use emacs with python mode
#    and color coding.
# -- You will probably need to edit almost every
#    single line of the go() function.
# -- If you need help on the individual function calls,
#    then in the pyraf prompt, import the module and
#    then print the documentation for that function:
#    --> print kai.kailog.__doc__
#    --> print range.__doc__
#
##################################################

# Import python and iraf modules
from pyraf import iraf as ir
import numpy as np
import os, sys
import glob

# Import our own custom modules
from kai.reduce import calib
from kai.reduce import sky
from kai.reduce import data
from kai.reduce import util
from kai.reduce import dar
from kai.reduce import kai_util
from kai import instruments


##########
# Change the epoch, instrument, and distortion solution.
##########
epoch = '19apr21'
nirc2 = instruments.NIRC2()  # don't necessarily need this as NIRC2 is the default.

##########
# Make electronic logs
#    - run this first thing for a new observing run.
##########
def makelog():
    """Make an electronic log from all the files in the ../raw/ directory.
    The file will be called kai.log and stored in the same directory.

    @author Jessica Lu
    @author Sylvana Yelda
    """
    kai_util.makelog('../raw', instrument=nirc2)

    # If you are reducing OSIRIS, you need to flip the images first. 
    # raw_files = glob.glob('../raw/i*.fits')
    # osiris.flip_images(raw_files)

    # Download weather data we will need.
    dar.get_atm_conditions('2020')

    return
    

###############
# Analyze darks
###############
def analyze_darks():
    """Analyze the dark_calib results
    """
    util.mkdir('calib')
    os.chdir('calib')

    first_dark = 16
    calib.analyzeDarkCalib(first_dark)  # Doesn't support OSIRIS yet

    os.chdir('../')

##########
# Reduce
##########
def go():
    """Do the full data reduction.

    @author Jessica Lu
    @author Sylvana Yelda
    """

    ####################
    #
    # Calibration files:
    #     everything created under calib/
    #
    ####################
    # Darks - created in subdir darks/
    # Files n0001 - n0009
    #  - darks needed to make bad pixel mask
    #  - store the resulting dark in the file name that indicates the
    #    integration time (2.8s) and the coadds (10ca).
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 
    darkFiles = list(range(1, 9+1))
    calib.makedark(darkFiles, 'dark_2.8s_10ca.fits', instrument=nirc2)

    # Flats - created in subdir flats/
    # Files n0037,n0039,n0041,n0043,n0045: lamps off for flats at K
    # Files n0038,n0040,n0042,n0044,n0046: lamps on for flats at K
    offFiles = [37, 39, 41, 43, 45]
    onFiles = [38, 40, 42, 44, 46]
    calib.makeflat(onFiles, offFiles, 'flat_kp.fits', instrument=nirc2)
    
    # Masks (assumes files were created under calib/darks/ and calib/flats/)
    calib.makemask('dark_2.8s_10ca.fits', 'flat_kp.fits',
                   'supermask.fits', instrument=nirc2)


    ####################
    #
    # Galactic Center
    #
    ####################

    ##########
    # K-band reduction
    ##########
    #util.mkdir('kp')
    #os.chdir('kp')

    # Nite 1:
    #    SCI frames (don't forget to add 1 at end): 108-237
    #    reference star position in first image:  [407.77, 673.87]
    #    use a different Strehl star at position: [521.79, 398.92]
    #    Sky frames (don't forget to add 1 at end): 252-261
    #
    #    -- If you have more than one position angle, make sure to
    #       clean them seperatly.
    #    -- Strehl and Ref src should be the pixel coordinates of a bright
    #       (but non saturated) source in the first exposure of sci_files.
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 
    sci_files1 = list(range(108, 237+1))
    refSrc1 = [407.77, 673.87]
    strSrc1 = [521.79, 398.92]
    sky_files1 = list(range(252, 261+1))
    #
    sky.makesky(sky_files1, 'nite1', 'kp', instrument=nirc2)
    data.clean(sci_files1, 'nite1', 'kp', refSrc1, strSrc1, instrument=nirc2) 
    

    # Nite 2:
    #    SCI frames (don't forget to add 1 at end): 1108-1237
    #    reference star position in first image:  [407.77, 673.87]
    #    use a different Strehl star at position: [521.79, 398.92]
    #    Sky frames (don't forget to add 1 at end): 1252-1261 and 1400-1405
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 
    sci_files2 = list(range(1108, 1237+1))
    refSrc2 = [407.77, 673.87]
    strSrc2 = [521.79, 398.92]
    sky_files2 = list(range(1252, 1261+1)) + list(range(1400, 1405+1))
    #
    sky.makesky(sky_files2, 'nite2', 'kp', instrument=nirc2)
    data.clean(sci_files2, 'nite2', 'kp', refSrc2, strSrc2, instrument=nirc2) 

    # Combine:
    #    Combine all of them together (from both nights).
    #    Do frame selection (trim=1).
    #    Weight by Strehl (weight='strehl').
    #    Make 3 submaps.
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 
    sci_files = sci_files1 + sci_files2
    data.calcStrehl(sci_files, 'kp', instrument=nirc2)
    data.combine(sci_files, 'kp', '06junlgs', trim=1, weight='strehl',
                   submaps=3, instrument=nirc2)

    os.chdir('../')

    
    ############
    # Lp
    ############

    #util.mkdir('lp')
    #os.chdir('lp')

    # Nite 1
    sky_files = list(range(38, 53+1))
    sci_files = list(range(87, 94+1)) + list(range(115,122+1))
    refSrc = [695., 543.]
    strSrc = [680., 839.]

    sky.makesky_lp(sky_files, 'nite1', 'lp')
    data.clean(sci_files, 'nite1', 'lp', refSrc, strSrc) 
    data.calcStrehl(sci_files, 'lp')
    data.combine(sci_files, 'lp', '06junlgs', 
                   trim=0, weight=None, submaps=3)

    os.chdir('../')

    #
    # Could do other wavelengths or different fields (e.g. Arches).
    # Can pass in a field flag, for example:
    # data.clean(arch_files,'nite1','kp',refSrc,strSrc,
    #              field='f1')
    # data.combine(arch_files,'kp','06maylgs1',trim=1,weight='strehl'
    #                field='f1',submaps=3)
    #

############
# Jackknife
############
def jackknife():
    """
    Perform the Jackknife data reduction. The methodology is as follows:
    For N total images, there are N stacked images created by data.combine(),
    each of these combo frames consists of N-1 images. The suffix '_{i}' is
    given to each jackknife combo frame (i.e. mag27maylgs_12_ks.fits would be
    the 12th jackknife combo frame created). Submaps flag has been turned off
    for this function, as it is not needed.
    """

    ##########
    # Ks-band reduction
    ##########

    # Nite 1
    sci_files1 = list(range(173, 177+1))
    sky_files1 = list(range(206, 215+1))
    refSrc1 = [385., 440.] #This is the target nearest to center

    sky.makesky(sky_files1, 'nite1', 'ks', instrument=nirc2)
    data.clean(sci_files1, 'nite1', 'ks', refSrc1, refSrc1, instrument=nirc2)


    # Nite 2
    sci_files2 = list(range(195, 203+1))
    sky_files2 = list(range(206, 215+1))
    refSrc2 = [387., 443.] #This is the target nearest to center

    sky.makesky(sky_files2, 'nite2', 'ks', instrument=nirc2)
    data.clean(sci_files2, 'nite2', 'ks', refSrc2, refSrc2, instrument=nirc2)
    #-----------------

    sci_files = sci_files1 + sci_files2

    for i in enumerate(sci_files, start=1):
        jack_list = sci_files[:]
        jack_list.remove(i[1])
        data.calcStrehl(jack_list, 'ks', instrument=nirc2)
        data.combine(jack_list, 'ks', '27maylgs', trim=1, weight='strehl',
                    instrument=nirc2, outSuffix='_' + str(i[0]))
        os.chdir('reduce')
