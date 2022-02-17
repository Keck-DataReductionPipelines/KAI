# Copied from /u/jlu/data/microlens/20aug22os/reduce/reduce.py

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
epoch = '20aug22os'
osiris = instruments.OSIRIS()

##########
# Make electronic logs
#    - run this first thing for a new observing run.
##########
def makelog_and_prep_images():
    """Make an electronic log from all the files in the ../raw/ directory.
    The file will be called kai.log and stored in the same directory.

    @author Jessica Lu
    @author Sylvana Yelda
    """
    kai_util.makelog('../raw', instrument=osiris)

    # If you are reducing OSIRIS, you need to flip the images first. 
    raw_files = glob.glob('../raw/i*.fits')
    osiris.flip_images(raw_files)

    # Download weather data we will need.
    dar.get_atm_conditions('2020')

    return
    

###############
# Analyze darks
###############
# def analyze_darks():
#     """Analyze the dark_calib results
#     """
#     util.mkdir('calib')
#     os.chdir('calib')
# 
#     first_dark = 16
#     calib.analyzeDarkCalib(first_dark)  # Doesn't support OSIRIS yet
# 
#     os.chdir('../')

##########
# Reduce
##########
def go_calib():
    """Do the calibration reduction.

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
    #  - darks needed to make bad pixel mask
    #  - store the resulting dark in the file name that indicates the
    #    integration time (2.8s) and the coadds (10ca).
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 
    darkFiles = ['i200822_s003{0:03d}_flip'.format(ii) for ii in range(23, 27+1)]
    calib.makedark(darkFiles, 'dark_39.832s_1ca_6rd.fits', instrument=osiris)

    darkFiles = ['i200822_s003{0:03d}_flip'.format(ii) for ii in range(28, 32+1)]
    calib.makedark(darkFiles, 'dark_5.901s_1ca_4rd.fits', instrument=osiris)

    darkFiles = ['i200822_s020{0:03d}_flip'.format(ii) for ii in range(2, 10+1)]
    calib.makedark(darkFiles, 'dark_11.802s_4ca_4rd.fits', instrument=osiris)

    darkFiles = ['i200822_s021{0:03d}_flip'.format(ii) for ii in range(2, 10+1)]
    calib.makedark(darkFiles, 'dark_5.901s_8ca_1rd.fits', instrument=osiris)

    # Flats - created in subdir flats/
    offFiles = ['i200822_s003{0:03d}_flip'.format(ii) for ii in range(3, 21+1, 2)]
    onFiles  = ['i200822_s003{0:03d}_flip'.format(ii) for ii in range(4, 22+1, 2)]
    calib.makeflat(onFiles, offFiles, 'flat_kp_tdOpen.fits', instrument=osiris)

    # Masks (assumes files were created under calib/darks/ and calib/flats/)
    calib.makemask('dark_39.832s_1ca_6rd.fits', 'flat_kp_tdOpen.fits',
                   'supermask.fits', instrument=osiris)

def go():
    """
    Do the full data reduction.
    """
    ##########
    #
    # MB19284
    #
    ##########

    ##########
    # Kp-band reduction
    ##########

    target = 'mb19284'
    sci_files = ['i200822_a011{0:03d}_flip'.format(ii) for ii in range(2, 5+1)]
    sci_files += ['i200822_a012{0:03d}_flip'.format(ii) for ii in range(2, 25+1)]
    sky_files = ['i200822_a018{0:03d}_flip'.format(ii) for ii in range(2, 6+1)]
    refSrc = [917.75, 1033.5] # This is the target
    # Alternative star to try (bright star to bottom of target): [1015, 581.9]
    
    sky.makesky(sky_files, target, 'kp_tdOpen', instrument=osiris)
    data.clean(sci_files, target, 'kp_tdOpen', refSrc, refSrc, field=target, instrument=osiris)
    data.calcStrehl(sci_files, 'kp_tdOpen', field=target, instrument=osiris)
    data.combine(sci_files, 'kp_tdOpen', epoch, field=target,
                     trim=0, weight='strehl', submaps=3, instrument=osiris)

    ##########
    #
    # KB200101
    #
    ##########

    ##########
    # Kp-band reduction
    ##########

    #    -- If you have more than one position angle, make sure to
    #       clean them seperatly.
    #    -- Strehl and Ref src should be the pixel coordinates of a bright
    #       (but non saturated) source in the first exposure of sci_files.
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 
    target = 'kb200101'
    sci_files = ['i200822_a014{0:03d}_flip'.format(ii) for ii in range(2, 28+1)]
    sci_files += ['i200822_a015{0:03d}_flip'.format(ii) for ii in range(2, 5+1)]
    sci_files += ['i200822_a016{0:03d}_flip'.format(ii) for ii in range(2, 5+1)]
    sky_files = ['i200822_a017{0:03d}_flip'.format(ii) for ii in range(2, 6+1)]
    refSrc = [975, 1006] # This is the target
    # Alternative star to try (bright star to right of target): [1158, 994]
    
    sky.makesky(sky_files, target, 'kp_tdOpen', instrument=osiris)
    data.clean(sci_files, target, 'kp_tdOpen', refSrc, refSrc, field=target, instrument=osiris)
    data.calcStrehl(sci_files, 'kp_tdOpen', field=target, instrument=osiris)
    data.combine(sci_files, 'kp_tdOpen', epoch, field=target,
                     trim=1, weight='strehl', submaps=3, instrument=osiris)

############
# Jackknife
############
def jackknife():
    """
    Perform the Jackknife data reduction. The methodology is as follows:
    For N total images, there are N stacked images created by data.combine(),
    each of these combo frames consists of N-1 images. The suffix '_{i}' is
    given to each jackknife combo frame (i.e. mag_OB06284_17_kp_tdOpen.fits 
    would be 17th jackknife combo frame created). Submaps flag has been turned
    off for this function, as it is not needed.
    """
    ##########
    #
    # OSIRIS Format
    #
    ##########

    ##########
    # Kp-band reduction
    ##########

    # Nite 1
    target = 'OB06284'
    sci_files1 = ['i200810_a004{0:03d}_flip'.format(ii) for ii in range(2, 26+1)]
    sky_files1 = ['i200810_a007{0:03d}_flip'.format(ii) for ii in range(2, 6+1)]
    refSrc1 = [1071., 854.] # This is the target

    sky.makesky(sky_files1, target, 'kp_tdOpen', instrument=osiris)
    data.clean(sci_files1, target, 'kp_tdOpen', refSrc, refSrc1, field=target, instrument=osiris)

    # Nite 2
    sci_files2 = ['i200810_a004{0:03d}_flip'.format(ii) for ii in range(2, 26+1)]
    sky_files2 = ['i200810_a007{0:03d}_flip'.format(ii) for ii in range(2, 6+1)]
    refSrc2 = [1071., 854.] # This is the target

    sky.makesky(sky_files2, target, 'kp_tdOpen', instrument=osiris)
    data.clean(sci_files2, target, 'kp_tdOpen', refSrc, refSrc2, field=target, instrument=osiris)
    #-----------------

    sci_files = sci_files1 + sci_files2

    for i in enumerate(sci_files, start=1):
        jack_list = sci_files[:]
        jack_list.remove(i[1])
        data.calcStrehl(jack_list, 'kp_tdOpen', field=target, instrument=osiris)
        data.combine(jack_list, 'kp_tdOpen', epoch, field=target,
                     trim=0, weight='strehl', instrument=osiris, outSuffix=str(i[0]))
        os.chdir('reduce')
