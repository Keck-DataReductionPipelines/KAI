import os, sys
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import stats
import math
from pyraf import iraf as ir
from . import kai_util
from kai.reduce import util
from kai import instruments
from kai import strehl
import time
import pdb
import numpy as np
from . import dar
from . import bfixpix
import subprocess
import copy
import shutil
from datetime import datetime

module_dir = os.path.dirname(__file__)

# Remember to change these if you are going to use the wide camera.
# You can change them in your reduce.py file after importing data.py
# Narrow Camera
distCoef = ''
# distXgeoim = module_dir + '/distortion/nirc2_narrow_xgeoim.fits'
# distYgeoim = module_dir + '/distortion/nirc2_narrow_ygeoim.fits'
#changing to test new distortion solutions right now
#distXgeoim = '/g/lu/data/m53/2015may/dist_sol/Leg_6r9_may_X.fits'
#distYgeoim = '/g/lu/data/m53/2015may/dist_sol/Leg_6r9_may_Y.fits'
# Wide Camera
#distCoef = module_dir + '/distortion/coeffs/nirc2_cubic_wide'
#distXgeoim = ''
#distYgeoim = ''
# Wide Camera Hai Fu
#distCoef = ' '
#distXgeoim = module_dir + '/distortion/nirc2_wide_X_distortion.fits'
#distYgeoim = module_dir + '/distortion/nirc2_wide_Y_distortion.fits'

supermaskName = 'supermask.fits'
outputVerify = 'ignore'


def clean(files, nite, wave, refSrc, strSrc, badColumns=None, field=None,
          skyscale=False, skyfile=None, angOff=0.0, cent_box=12,
          fixDAR=True,
          raw_dir=None, clean_dir=None,
          instrument=instruments.default_inst, check_ref_loc=True):
    """
    Clean near infrared NIRC2 or OSIRIS images.

    This program should be run from the reduce/ directory.
    Example directory structure is:
    calib/
        flats/
        flat_kp.fits
        flat.fits (optional)
        masks/
        supermask.fits
    kp/
        sci_nite1/
        sky_nite1/
        sky.fits

    All output files will be put into clean_dir (if specified, otherwise
    ../clean/) in the following structure:
    kp/
        c*.fits
        distort/
        cd*.fits
        weight/
        wgt*.fits

    The clean directory may be optionally modified to be named
    <field_><wave> instead of just <wave>. So for instance, for Arches
    field #1 data reduction, you might call clean with: field='arch_f1'.

    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    nite : str
        Name for night of observation (e.g.: "nite1"), used as suffix
        inside the reduce sub-directories.
    wave : str
        Name for the observation passband (e.g.: "kp"), used as
        a wavelength suffix
    field : str, default=None
        Optional prefix for clean directory and final
        combining. All clean files will be put into <field_><wave>. You
        should also pass the same into combine(). If set to None (default)
        then only wavelength is used.
    skyscale : bool, default=False
        Whether or not to scale the sky files to the common median.
        Turn on for scaling skies before subtraction.
    skyfile : str, default=''
        An optional file containing image/sky matches.
    angOff : float, default = 0
        An optional absolute offset in the rotator
        mirror angle for cases (wave='lp') when sky subtraction is done with
        skies taken at matching rotator mirror angles.
    cent_box: int (def = 12)
        the box to use for better centroiding the reference star
    badColumns : int array, default = None
        An array specifying the bad columns (zero-based).
        Assumes a repeating pattern every 8 columns.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    clean_dir : str, optional
        Directory where clean files will be stored. By default,
        assumes that clean files will be stored in '../clean'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    
    # Make sure directory for current passband exists and switch into it
    util.mkdir(wave)
    os.chdir(wave)
    
    # Determine directory locatons
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')
    
    sciDir = waveDir + '/sci_' + nite + '/'
    util.mkdir(sciDir)
    ir.cd(sciDir)

    # Set location of raw data
    rawDir = rootDir + 'raw/'
    
    # Check if user has specified a specific raw directory
    if raw_dir is not None:
        rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')
    
    # Setup the clean directory
    cleanRoot = rootDir + 'clean/'
    
    # Check if user has specified a specific clean directory
    if clean_dir is not None:
        cleanRoot = util.trimdir(os.path.abspath(clean_dir) + '/')
    
    if field is not None:
        clean = cleanRoot + field + '_' + wave + '/'
    else:
        clean = cleanRoot + wave + '/'
    
    distort = clean + 'distort/'
    weight = clean + 'weight/'
    masks = clean + 'masks/'
    
    util.mkdir(cleanRoot)
    util.mkdir(clean)
    util.mkdir(distort)
    util.mkdir(weight)
    util.mkdir(masks)
    
    # Open a text file to document sources of data files
    data_sources_file = open(clean + 'data_sources.txt', 'a')
    
    try:
        # Setup flat. Try wavelength specific, but if it doesn't
        # exist, then use a global one.
        flatDir = redDir + 'calib/flats/'
        flat = flatDir + 'flat_' + wave + '.fits'
        if not os.access(flat, os.F_OK):
            flat = flatDir + 'flat.fits'

        # Bad pixel mask
        _supermask = redDir + 'calib/masks/' + supermaskName

        # Determine the reference coordinates for the first image.
        # This is the image for which refSrc is relevant.
        firstFile = instrument.make_filenames([files[0]], rootDir=rawDir)[0]
        hdr1 = fits.getheader(firstFile, ignore_missing_end=True)
        radecRef = [float(hdr1['RA']), float(hdr1['DEC'])]
        aotsxyRef = kai_util.getAotsxy(hdr1)

        # Setup a Sky object that will figure out the sky subtraction
        skyDir = waveDir + 'sky_' + nite + '/'
        skyObj = Sky(sciDir, skyDir, wave, scale=skyscale,
                     skyfile=skyfile, angleOffset=angOff,
                     instrument=instrument)

        # Prep drizzle stuff
        # Get image size from header - this is just in case the image
        # isn't 1024x1024 (e.g., NIRC2 sub-arrays). Also, if it's
        # rectangular, choose the larger dimension and make it square
        imgsizeX = float(hdr1['NAXIS1'])
        imgsizeY = float(hdr1['NAXIS2'])

        distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr1)
        if (imgsizeX >= imgsizeY):
            imgsize = imgsizeX
        else:
            imgsize = imgsizeY
        setup_drizzle(imgsize)

        ##########
        # Loop through the list of images
        ##########
        for f in files:
            # Define filenames
            _raw = instrument.make_filenames([f], rootDir=rawDir)[0]
            _cp = instrument.make_filenames([f])[0]
            _ss = instrument.make_filenames([f], prefix='ss')[0]
            _ff = instrument.make_filenames([f], prefix='ff')[0]
            _ff_f = _ff.replace('.fits', '_f.fits')
            _ff_s = _ff.replace('.fits', '_s.fits')
            _bp = instrument.make_filenames([f], prefix='bp')[0]
            _cd = instrument.make_filenames([f], prefix='cd')[0]
            _ce = instrument.make_filenames([f], prefix='ce')[0]
            _cc = instrument.make_filenames([f], prefix='c')[0]
            _wgt = instrument.make_filenames([f], prefix='wgt')[0]
            _statmask = instrument.make_filenames([f], prefix='stat_mask')[0]
            _crmask = instrument.make_filenames([f], prefix='crmask')[0]
            _mask = instrument.make_filenames([f], prefix='mask')[0]
            _pers = instrument.make_filenames([f], prefix='pers')[0]
            _max = _cc.replace('.fits', '.max')
            _coo = _cc.replace('.fits', '.coo')
            _rcoo = _cc.replace('.fits', '.rcoo')
            _dlog_tmp = instrument.make_filenames([f], prefix='driz')[0]
            _dlog = _dlog_tmp.replace('.fits', '.log')
            
            out_line = '{0} from {1} ({2})\n'.format(_cc, _raw,
                                                     datetime.now())
            data_sources_file.write(out_line)

            # Clean up if these files previously existed
            util.rmall([_cp, _ss, _ff, _ff_f, _ff_s, _bp, _cd, _ce, _cc,
                        _wgt, _statmask, _crmask, _mask, _pers, _max, _coo, _rcoo, _dlog])

            ### Copy the raw file to local directory ###
            ir.imcopy(_raw, _cp, verbose='no')
            
            ### Make persistance mask ###
            # - Checked images, this doesn't appear to be a large effect.
            #clean_persistance(_cp, _pers, instrument=instrument)

            ### Sky subtract ###
            # Get the proper sky for this science frame.
            # It might be scaled or there might be a specific one for L'.
            sky = skyObj.getSky(_cp)

            ir.imarith(_cp, '-', sky, _ss)

            ### Flat field ###
            ir.imarith(_ss, '/', flat, _ff)

            ### Make a static bad pixel mask ###
            # _statmask = supermask + bad columns
            clean_get_supermask(_statmask, _supermask, badColumns)

            ### Fix bad pixels ###
            # Produces _ff_f file
            bfixpix.bfixpix(_ff, _statmask)
            util.rmall([_ff_s])

            ### Fix cosmic rays and make cosmic ray mask. ###
            clean_cosmicrays(_ff_f, _crmask, wave)

            ### Combine static and cosmic ray mask ###
            # This will be used in combine later on.
            # Results are stored in _mask, _mask_static is deleted.
            clean_makemask(_mask, _crmask, _statmask, wave, instrument=instrument)

            ### Background Subtraction ###
            bkg = clean_bkgsubtract(_ff_f, _bp)

            ### Drizzle individual file ###
            clean_drizzle(distXgeoim, distYgeoim, _bp, _ce, _wgt, _dlog, fixDAR=fixDAR, instrument=instrument)

            ### Make .max file ###
            # Determine the non-linearity level. Raw data level of
            # non-linearity is 12,000 but we subtracted
            # off a sky which changed this level. The sky is
            # scaled, so the level will be slightly different
            # for every frame.
            nonlinSky = skyObj.getNonlinearCorrection(sky)

            coadds = fits.getval(_ss, instrument.hdr_keys['coadds'])
            satLevel = (coadds*instrument.get_saturation_level()) - nonlinSky - bkg
            file(_max, 'w').write(str(satLevel))

            ### Rename and clean up files ###
            ir.imrename(_bp, _cd)
            # util.rmall([_cp, _ss, _ff, _ff_f])

            ### Make the *.coo file and update headers ###
            # First check if PA is not zero
            hdr = fits.getheader(_raw, ignore_missing_end=True)
            phi = instrument.get_position_angle(hdr)

            clean_makecoo(_ce, _cc, refSrc, strSrc, aotsxyRef, radecRef,
                          instrument=instrument, check_loc=check_ref_loc,
                          cent_box=cent_box)

            ### Move to the clean directory ###
            util.rmall([clean + _cc, clean + _coo, clean + _rcoo,
                        distort + _cd, weight + _wgt,
                        clean + _ce, clean + _max,
                        masks + _mask, _ce])

            os.rename(_cc, clean + _cc)
            os.rename(_cd, distort + _cd)
            os.rename(_wgt, weight + _wgt)
            os.rename(_mask, masks + _mask)
            os.rename(_max, clean + _max)
            os.rename(_coo, clean + _coo)
            os.rename(_rcoo, clean + _rcoo)

            # This just closes out any sky logging files.
            #skyObj.close()
        data_sources_file.close()
    finally:
        # Move back up to the original directory
        #skyObj.close()
        ir.cd('../')

    # Change back to original directory
    os.chdir('../')

def clean_get_supermask(_statmask, _supermask, badColumns):
    """
    Create temporary mask for each individual image that will contain the
    supermask plus the designated bad columns.

    _statmask -- output file containing supermask + bad columns
    """

    maskFits = fits.open(_supermask)

    # Check that we have some valid bad columns.
    if badColumns != None and len(badColumns) != 0:
        for cc in badColumns:
            if (cc < 0):
                continue

            # Make column index from 0-512 n steps of 8
            colIndex = np.arange(cc, 512, 8)
            maskFits[0].data[0:512,colIndex] = 1

    # Save to a temporary file.
    maskFits[0].writeto(_statmask, output_verify=outputVerify)

def clean_makemask(_mask, _mask_cosmic, _mask_static, wave,
                   instrument=instruments.default_inst):
    """
    _mask -- output name for final mask
    _mask_cosmic -- should contain cosmic ray mask
    _mask_static -- should contain supermask + bad columns

    Output:
    _mask is created to be supermask + bad columns + cosmic rays
    _mask will have 0=bad and 1=good pixels (as drizzle expects)
    _mask can be directly passed into drizzle
    """
    # Get the masks to combine
    staticMask = fits.getdata(_mask_static)
    cosmicMask = fits.getdata(_mask_cosmic)

    mask = staticMask + cosmicMask

    # check subarray
    if (instrument.name == 'NIRC2') and ('lp' in wave or 'ms' in wave) and (mask.shape[0] > 512):
        _lpmask = module_dir + '/masks/nirc2_lp_edgemask.fits'
        lpmask = fits.getdata(_lpmask)
        mask += lpmask

    # Set to 0 or 1 -- note they are inverted
    weightone = (mask == 0)
    weightzero = (mask != 0)

    # Drizzle expects 0 = bad, 1 = good pixels.
    outMask = np.zeros(mask.shape)
    outMask[weightone] = 1
    outMask[weightzero] = 0

    # Trim 12 rows from top and bottom b/c the distortion solution
    # introduces a torque to the image.
    if (instrument.name == 'NIRC2'):
        outMask[1012:1024,0:1024] = 0
        outMask[0:12,0:1024] = 0

    # Write out to file
    fits.writeto(_mask, outMask, output_verify=outputVerify)
    #outMask[0].writeto(_mask, output_verify=outputVerify)


def clean_lp(files, nite, wave, refSrc, strSrc, angOff, skyfile):
    """
    Only here for backwards compatability.
    You should use clean() instead.
    """
    clean(files, nite, wave, refSrc, strSrc,
          angOff=angOff, skyfile=skyfile)

def combine(files, wave, outroot, field=None, outSuffix=None,
            trim=False, weight=None, fwhm_max=0, submaps=0,
            fixDAR=True, mask=True,
            clean_dirs=None, combo_dir=None,
            instrument=instruments.default_inst):
    """
    Accepts a list of cleaned images and does a weighted combining after
    performing frame selection based on the Strehl and FWHM.
    
    Each image must have an associated *.coo file which gives the rough
    position of the reference source.
    
    Parameters
    ----------
    files : list of int
        Integer list of the files to include in combine. Does not require
        padded zeros.
    wave : str
        Name for the observation passband (e.g.: "kp", "lp", or "h"), used as
        a wavelength suffix
    outroot : str
        The output root name (e.g. '06jullgs'). The final combined file names
        will be <outroot>_<field>_<outSuffix>_<wave>.
        The <field> and <outSuffix> keywords are optional.
        
        Examples:
        06jullgs_kp for outroot='06jullgs' and wave='kp'
        06jullgs_arch_f1_kp for adding field='arch_f1'
    field : str, default=None
        Optional field name. Used to get to clean directory and also affects
        the final output file name.
    outSuffix : str
        Optional suffix used to modify final output file name.
        Can use suffix to indicate a night of observation (e.g.: "nite1").
    trim : bool, default=False
        Optional file trimming based on image quality. Default
        is False. Set to True to turn trimming on.
    weight : str, default=None
        Optional weighting. Set to 'strehl' to weight by Strehl, as found in
        strehl_source.txt file.
        OR set to a file name with the first column being the file name
        (e.g., c0021.fits) and the second column being the weight. Weights will
        be renormalized to sum to 1.0.
        Default = None, no weighting.
    fwhm_max : float, default=0
        The maximum allowed FWHM for keeping frames when trimming is turned on.
    submaps : int, default=0
        Set to the number of submaps to be made (def=0).
    fixDAR : bool, default=True
    mask : bool, default=True
    clean_dirs : list of str, optional
        List of directories where clean files are stored. Needs to be same
        length as files list. If not specified, by default assumes that
        clean files are stored in '../clean'.
    combo_dir : str, optional
        Directory where combo files will be stored. By default,
        assumes that combo files will be stored in '../combo'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    
    # Setup some files and directories
    redDir = util.getcwd()
    rootDir = util.trimdir( os.path.abspath(redDir + '../') + '/')
    
    # Determine clean directory and add field and suffixes to outroot
    cleanRoot = rootDir + 'clean/'
    
    if field is not None:
        cleanDir = cleanRoot + field + '_' + wave + '/'
        outroot += '_' + field
    else:
        cleanDir = cleanRoot + wave + '/'
    
    # If clean directories are specified for each file,
    # first tack on the field and wave to each path
    if clean_dirs is not None:
        # If incorrect number of clean directories specified, raise ValueError
        if len(clean_dirs) != len(files):
            err_str = 'Length of clean_dirs needs to match number of files, '
            err_str += str(len(files))
            
            raise ValueError(err_str)
        
        # Tack on field and wave to each path
        for clean_dir_index in range(len(clean_dirs)):
            cleanRoot = util.trimdir(
                            os.path.abspath(clean_dirs[clean_dir_index] + '/'))
            
            if field is not None:
                clean_dirs[clean_dir_index] = cleanRoot + '/' + field +\
                                              '_' + wave + '/'
            else:
                clean_dirs[clean_dir_index] = cleanRoot + '/' + wave + '/'
    
    if (outSuffix != None):
        outroot += '_' + outSuffix
    
    # Set up combo directory. This is the final output directory.
    comboDir = rootDir + 'combo/'
    
    if combo_dir is not None:
        comboDir = util.trimdir(os.path.abspath(combo_dir) + '/')
    
    util.mkdir(comboDir)
    
    
    # Make strings out of all the filename roots.
    allroots = instrument.make_filenames(files, prefix='')
    allroots = [aa.replace('.fits', '') for aa in allroots]
    
    # If clean directories were specified for each file, copy over
    # clean files to a common clean directory into new combo directory
    if clean_dirs is not None:
        # Save and make new clean directory, inside the new combo directory
        cleanDir = comboDir + 'clean/'
        
        field_wave_suffix = ''
        
        if field is not None:
            field_wave_suffix = field + '_' + wave + '/'
        else:
            field_wave_suffix = wave + '/'
        
        cleanDir += field_wave_suffix
        
        util.mkdir(cleanDir)
        util.mkdir(cleanDir + 'distort/')
        util.mkdir(cleanDir + 'masks/')
        util.mkdir(cleanDir + 'weight/')
        
        # Determine all unique clean directories, which we'll be sourcing from
        (unique_clean_dirs,
         unique_clean_dirs_index) = np.unique(clean_dirs, return_inverse=True)
        
        c_lis_file = open(cleanDir + 'c.lis', 'w')
        data_sources_file = open(cleanDir + 'data_sources.txt', 'w')
        
        # Go through each clean file and copy over the data files
        for cur_file_index in range(len(files)):
            cur_file_root = allroots[cur_file_index]
            
            source_clean_dir = unique_clean_dirs[
                unique_clean_dirs_index[cur_file_index]]
            source_file_root = cur_file_root
            
            # Change first digit of file names to be index of clean dir
            # i.e.: unique 1000s place digit for each night going into combo
            allroots[cur_file_index] =\
                str(unique_clean_dirs_index[cur_file_index]) + cur_file_root[1:]
            
            dest_clean_dir = cleanDir
            dest_file_root = allroots[cur_file_index]
            
            # Copy data files
            shutil.copy(source_clean_dir + 'c' + source_file_root + '.fits',
                        dest_clean_dir + 'c' + dest_file_root + '.fits')
            
            shutil.copy(source_clean_dir + 'c' + source_file_root + '.max',
                        dest_clean_dir + 'c' + dest_file_root + '.max')
            
            shutil.copy(source_clean_dir + 'c' + source_file_root + '.coo',
                        dest_clean_dir + 'c' + dest_file_root + '.coo')
            
            shutil.copy(source_clean_dir + 'distort/cd' + source_file_root + '.fits',
                        dest_clean_dir + 'distort/cd' + dest_file_root + '.fits')
            
            shutil.copy(source_clean_dir + 'masks/mask' + source_file_root + '.fits',
                        dest_clean_dir + 'masks/mask' + dest_file_root + '.fits')
            
            shutil.copy(source_clean_dir + 'weight/wgt' + source_file_root + '.fits',
                        dest_clean_dir + 'weight/wgt' + dest_file_root + '.fits')
            
            # Append file to c.lis and text list of data sources
            c_lis_file.write(dest_clean_dir + 'c' + dest_file_root + '.fits\n')
            
            out_line = '{0} from {1}{2} ({3})\n'.format(
                'c' + dest_file_root + '.fits',
                source_clean_dir, 'c' + source_file_root + '.fits',
                datetime.now())
            data_sources_file.write(out_line)
        
        c_lis_file.close()
        data_sources_file.close()
        
        # Copy over strehl source list(s) from clean directories
        
        # Need to rename file names in list to new names
        out_strehl_file = open(cleanDir + 'strehl_source.txt', 'w')
        
        # Go through each clean directory's strehl_source file
        for cur_clean_dir_index in range(0, len(unique_clean_dirs)):
            
            # Open existing Strehl file in clean directory
            with open(unique_clean_dirs[cur_clean_dir_index] +
                      'strehl_source.txt', 'r') as in_strehl_file:
                
                for line in in_strehl_file:
                    # Check for header line
                    if line[0] == '#':
                        # Don't skip header if it is first clean directory
                        if cur_clean_dir_index == 0:
                            out_strehl_file.write(line)
                            
                        # Otherwise skip to next line
                        continue
                    
                    # Correct file names and write to overall strehl file
                    corrected_line = 'c' + str(cur_clean_dir_index) + line[2:]
                    out_strehl_file.write(corrected_line)
        
        out_strehl_file.close()
    
    # Make a deep copy of all the root filenames    
    roots = copy.deepcopy(allroots) # This one will be modified by trimming
    
    # This is the output root filename
    _out = comboDir + 'mag' + outroot + '_' + wave
    _sub = comboDir + 'm' + outroot + '_' + wave
    
    ##########
    # Determine if we are going to trim and/or weight the files
    # when combining. If so, then we need to determine the Strehl
    # and FWHM for each image. We check strehl source which shouldn't
    # be saturated. *** Hard coded to strehl source ***
    ##########

    # Load the strehl_source.txt file
    if ((weight is not None) or
        os.path.exists(os.path.join(cleanDir,'strehl_source.txt'))):
        strehls, fwhm = loadStrehl(cleanDir, roots)
    else:
        # if the file doesn't exist don't use
        print('combine: the strehl_source file does not exist: '+os.path.join(cleanDir,'strehl_source.txt'))

        # fill out some variables for later use
        strehls = np.zeros(len(roots))-1.0
        fwhm = np.zeros(len(roots)) -1.0
        trim = False
    

    # Default weights
    # Create an array with length equal to number of frames used,
    # and with all elements equal to 1/(# of files)
    weights = np.array( [1.0/len(roots)] * len(roots) )

    ##########
    # Trimming
    ##########
    if trim:
        roots, strehls, fwhm, weights = trim_on_fwhm(roots, strehls, fwhm,
                                                     fwhm_max=fwhm_max)

    ##########
    # Weighting
    ##########
    if weight == 'strehl':
        weights = weight_by_strehl(roots, strehls)

    if ((weight is not None) and (weight is not 'strehl')):
        # Assume weight is set to a filename
        if not os.path.exists(weight):
            raise ValueError('Weights file does not exist, %s' % weight)

        print('Weights file: ', weight)

        weights = readWeightsFile(roots, weight)

    # Determine the reference image
    # refImage_index = 0    # Use the first image from night
    refImage_index = np.argmin(fwhm)    # Use the lowest FWHM frame
    refImage = cleanDir + 'c' + roots[refImage_index] + '.fits'
    print('combine: reference image - %s' % refImage)

    ##########
    # Write out a log file. With a list of images in the
    # final combination.
    ##########
    combine_log(_out, roots, strehls, fwhm, weights)

    # See if all images are at same PA, if not, rotate all to PA = 0
    # temporarily. This needs to be done to get correct shifts.
    diffPA = combine_rotation(cleanDir, roots, instrument=instrument)

    # Make a table of coordinates for the reference source.
    # These serve as initial estimates for the shifts.
    #combine_ref(_out + '.coo', cleanDir, roots, diffPA, refImage_index)
    combine_coo(_out + '.coo', cleanDir, roots, diffPA, refImage_index)

    # Keep record of files that went into this combine
    combine_lis(_out + '.lis', cleanDir, roots, diffPA)

    # Register images to get shifts.
    shiftsTab = combine_register(_out, refImage, diffPA)

    # Determine the size of the output image from max shifts
    xysize = combine_size(shiftsTab, refImage, _out, _sub, submaps)

    ##########
    # Sort frames -- recall that submaps assume sorted by FWHM.
    ##########
    roots, strehls, fwhm, weights, shiftsTab = sort_frames(roots, strehls, fwhm, weights, shiftsTab)

    # Combine all the images together.
    combine_drizzle(xysize, cleanDir, roots, _out, weights, shiftsTab,
                    wave, diffPA, fixDAR=fixDAR, mask=mask, instrument=instrument)

    # Now make submaps
    if (submaps > 0):
        combine_submaps(xysize, cleanDir, roots, _sub, weights,
                        shiftsTab, submaps, wave, diffPA, fixDAR=fixDAR,
                        mask=mask, instrument=instrument)

    # Remove *.lis_r file & rotated rcoo files, if any - these
    # were just needed to get the proper shifts for xregister
    _lisr = _out + '.lis_r'
    util.rmall([_lisr])
    for i in range(len(allroots)):
        _rcoo = cleanDir + 'c' + str(allroots[i]) + '.rcoo'
        util.rmall([_rcoo])
    
    # Change back to original directory
    os.chdir('../')

def rot_img(root, phi, cleanDir):
    """Rotate images to PA=0 if they have a different PA from one
    another. If the entire data set is taken at a single PA, leave
    it as is. Do this only if set includes various PAs.
    """
    pa = str(phi)
    ir.unlearn('rotate')
    ir.rotate.verbose = 'no'
    ir.rotate.boundary = 'constant'
    ir.rotate.constant = 0
    ir.rotate.interpolant = 'spline3'

    inCln = cleanDir + 'c' + root + '.fits'
    outCln = cleanDir + 'r' + root + '.fits'

    util.rmall([outCln])

    if (phi != 0):
        print('Rotating frame: ',root)
        ir.rotate(inCln, outCln, pa)
    else:
        ir.imcopy(inCln, outCln, verbose='no')

    return

def gcSourceXY(name, label_file='/Users/jlu/data/gc/source_list/label.dat'):
    """
    Queries label.dat for the xy offset from Sgr A* (in arcsec)
    for the star given as an input
    
    Parameters
    ----------
    name : str
        Name of a star (e.g. 'irs16NE')
    label_file : str, default='/Users/jlu/data/gc/source_list/label.dat'
        Full path of label.dat file to search
    
    Returns
    -------
    pos : float list (2 elements)
        x and y offset from Sgr A* in arcsec
    """

    # Read in label.dat
    table = Table.read(label_file, format='ascii')
    cols = list(table.columns.keys())

    nameCol = table[cols[0]]
    names = [n.strip() for n in nameCol]

    try:
        id = names.index(name)

        x = table[cols[2]][id]
        y = table[cols[3]][id]
    except ValueError as e:
        print('Could not find source ' + name + ' in label.dat.')
        x = 0
        y = 0

    return [x,y]


def calcStrehl(files, wave,
               clean_dir=None, field=None,
               instrument=instruments.default_inst):
    """
    Make Strehl and FWHM table on the strehl source for all
    cleaned files.

    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    wave : str
        Name for the observation passband (e.g.: "kp"), used as
        a wavelength suffix
    field : str, default=None
        Optional prefix for clean directory and final
        combining. All clean files will be put into <field_><wave>. You
        should also pass the same into combine(). If set to None (default)
        then only wavelength is used.
    clean_dir : str, optional
        Directory where clean files will be stored. By default,
        assumes that clean files will be stored in '../clean'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    
    # Make sure directory for current passband exists and switch into it
    util.mkdir(wave)
    os.chdir(wave)
    
    # Determine directory locatons
    waveDir = util.getcwd()
    redDir = util.trimdir( os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir( os.path.abspath(redDir + '../') + '/')
    
    
    # Setup the clean directory
    cleanRoot = rootDir + 'clean/'
    
    # Check if user has specified a specific clean directory
    if clean_dir is not None:
        cleanRoot = util.trimdir(os.path.abspath(clean_dir) + '/')
    
    if field is not None:
        cleanDir = cleanRoot + field + '_' + wave + '/'
    else:
        cleanDir = cleanRoot + wave + '/'
    
    
    # Make a list of all the images
    clis_file = cleanDir + 'c.lis'
    strehl_file = cleanDir + 'strehl_source.txt'
    util.rmall([clis_file, strehl_file])

    clean_files = instrument.make_filenames(files, rootDir=cleanDir, prefix='c')

    # Keep a record of the cleaned files. 
    _clis = open(clis_file, 'w')
    for c in clean_files:
        _clis.write('%s\n' % c)
    _clis.close()

    # Calculate Strehl, FWHM
    strehl.calc_strehl(clean_files, strehl_file, instrument=instrument)
    
    # Check that the number of lines in the resulting strehl file
    # matches the number of images we have. If not, some of the images
    # are bad and were dropped.
    strehlTable = Table.read(strehl_file, format='ascii', header_start=None)
    cols = list(strehlTable.columns.keys())

    if len(clean_files) != len(strehlTable):
        print(len(clean_files), len(strehlTable))
        # Figure out the dropped files.
        droppedFiles = []
        for cc in clean_files:
            root = os.path.split(cc)[-1]
            
            foundIt = False
            for ss in strehlTable[cols[0]]:
                if root in ss:
                    foundIt = True
                    continue
            if foundIt == False:
                droppedFiles.append(root)

        raise RuntimeError('calcStrehl: Strehl widget lost files: ',
                           droppedFiles)
    
    # Switch back to parent directory
    os.chdir('../')

def weight_by_strehl(roots, strehls):
    """
    Calculate weights based on the strehl of each image.
    This does some intelligent handling for REALLY bad data quality.
    """
    # Set negative Strehls to the lowest detected strehl.
    bidx = (np.where(strehls <= 0))[0]
    gidx = (np.where(strehls > 0))[0]
    if len(bidx) > 0:
        badroots = [roots[i] for i in bidx]
        print('Found files with incorrect Strehl. May be incorrectly')
        print('weighted. Setting weights to minimum weight. ')
        print('\t' + ','.join(badroots))

    strehl_min = strehls[gidx].min()
    strehls[bidx] = strehl_min

    # Now determine a fractional weight
    wgt_tot = sum(strehls)
    weights = strehls / wgt_tot

    return weights


def trim_on_fwhm(roots, strehls, fwhm, fwhm_max=0):
    """
    Take a list of files and trim based on the FWHM. All files that have a
    FWHM < 1.25 * FWHM.min()
    are kept.

    The returned arrays contain only those files that pass the above criteria.
    """
    # Trim level (fwhm) can be passed in or determined
    # dynamically.
    
    if (fwhm_max == 0):
        # Determine the minimum FWHM
        idx = np.where(fwhm > 0)
        fwhm_min = fwhm[idx].min()

        # Maximum allowed FWHM to keep frame
        fwhm_max = 1.25 * fwhm_min

    # Pull out those we want to include in the combining
    keep = np.where((fwhm <= fwhm_max) & (fwhm > 0))[0]
    strehls = strehls[keep]
    fwhm = fwhm[keep]
    roots = [roots[i] for i in keep]
    weights = np.array( [1.0/len(roots)] * len(roots) )

    print('combine: Keeping %d frames with FWHM < %4.1f' \
        % (len(roots), fwhm_max))

    return (roots, strehls, fwhm, weights)


def readWeightsFile(roots, weightFile):
    """
    Expects a file of the format:
    column1 = file name (e.g. c0001.fits).
    column2 = weights.
    """

    weightsTable = trim_table_by_name(roots, weightFile)

    weights = weightsTable['col2']
    
    # Renormalize so that weights add up to 1.0
    weights /= weights.sum()
    
    # Double check that we have the same number of
    # lines in the weightsTable as files.
    if (len(weights) != len(roots)):
        print('Wrong number of lines in  ' + weightFile)

    return weights


def loadStrehl(cleanDir, roots):
    """
    Load Strehl and FWHM info. The file format will be
    column1 = name of cleaned fits file (e.g. c0001.fits).
              Expects single character before a 4 digit number.
    column2 = strehl
    column3 = RMS error (nm)
    column4 = FWHM (mas)
    column5 = MJD (UT)
    """
    _strehl = cleanDir + 'strehl_source.txt'
    
    # Read in file and get strehls and FWHMs
    strehlTable = trim_table_by_name(roots, _strehl)
    strehls = strehlTable['col2']
    fwhm = strehlTable['col4']
    
    # Double check that we have the same number of
    # lines in the strehlTable as files.
    if (len(strehls) != len(roots)):
        print('Wrong number of lines in  ' + _strehl)

    return (strehls, fwhm)

def trim_table_by_name(outroots, tableFileName):
    """
    Takes a list of values (listed in tableFileName) and trim them down based on
    the desired output list of root files names (outroots).
    """
    table = Table.read(tableFileName, format='ascii', header_start=None)
    
    good = np.zeros(len(table), dtype=bool)
    
    for rr in range(len(outroots)):
        for ii in range(len(table)):
            if outroots[rr] in table[ii][0]:
                good[ii] = True

    newtable = table[good]

    return newtable


def combine_drizzle(imgsize, cleanDir, roots, outroot, weights, shifts,
                    wave, diffPA, fixDAR=True, mask=True,
                    instrument=instruments.default_inst):
    _fits = outroot + '.fits'
    _tmpfits = outroot + '_tmp.fits'
    _wgt = outroot + '_sig.fits'
    _dlog = outroot + '_driz.log'
    _maxFile = outroot + '.max'

    util.rmall([_fits, _tmpfits, _wgt, _dlog])

    # Make directory for individually drizzled and pre-shifted images.
    #util.mkdir(cleanDir + 'shifted')

    # Prep drizzle stuff
    setup_drizzle(imgsize)

    # BUG: with context... when many files are drizzled
    # together, a new, bigger context file is created, but this
    # fails with a Bus error.
    #ir.drizzle.outcont = _ctx
    ir.drizzle.outcont = ''
    ir.drizzle.fillval = 0.

    satLvl_combo = 0.0

    # Set a cleanDir variable in IRAF. This avoids the long-filename problem.
    ir.set(cleanDir=cleanDir)
    
    # Variable to store weighted sum of MJDs
    mjd_weightedSum = 0.0

    # Get the distortion maps for this instrument.
    hdr0 = fits.getheader(cleanDir + 'c' + roots[0] + '.fits')
    distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr0)
    
    print('combine: drizzling images together')
    f_dlog = open(_dlog, 'a')
    for i in range(len(roots)):
        # Cleaned image
        _c = cleanDir + 'c' + roots[i] + '.fits'
        _c_ir = _c.replace(cleanDir, 'cleanDir$')

        # Cleaned but distorted image
        _cd = cleanDir + 'distort/cd' + roots[i] + '.fits'
        _cdwt = cleanDir + 'weight/cdwt.fits'
        _cd_ir = _cd.replace(cleanDir, 'cleanDir$')
        _cdwt_ir = _cdwt.replace(cleanDir, 'cleanDir$')

        util.rmall([_cdwt])

        # Multiply each distorted image by it's weight
        ir.imarith(_cd_ir, '*', weights[i], _cdwt_ir)

        # Fix the ITIME header keyword so that it matches (weighted).
        # Drizzle will add all the ITIMEs together, just as it adds the flux.
        itime = fits.getval(_cdwt, instrument.hdr_keys['itime'])
        itime *= weights[i]
        fits.setval(_cdwt, instrument.hdr_keys['itime'], value=itime)

        # Get pixel shifts
        xsh = shifts[i][1]
        ysh = shifts[i][2]

        # Read in PA of each file to feed into drizzle for rotation
        hdr = fits.getheader(_c, ignore_missing_end=True)
        phi = instrument.get_position_angle(hdr)
        if (diffPA == 1):
            ir.drizzle.rot = phi

        if (fixDAR == True):
            darRoot = _cdwt.replace('.fits', 'geo')
            (xgeoim, ygeoim) = dar.darPlusDistortion(_cdwt, darRoot,
                                                     xgeoim=distXgeoim,
                                                     ygeoim=distYgeoim,
                                                     instrument=instrument)

            xgeoim = xgeoim.replace(cleanDir, 'cleanDir$')
            ygeoim = ygeoim.replace(cleanDir, 'cleanDir$')
            ir.drizzle.xgeoim = xgeoim
            ir.drizzle.ygeoim = ygeoim
        else:
            ir.drizzle.xgeoim = distXgeoim
            ir.drizzle.ygeoim = distYgeoim

        
        # Read in MJD of current file from FITS header
        mjd = float(hdr['MJD-OBS'])
        mjd_weightedSum += weights[i] * mjd
        
        # Drizzle this file ontop of all previous ones.
        f_dlog.write(time.ctime())

        if (mask == True):
            _mask = 'cleanDir$masks/mask' + roots[i] + '.fits'
        else:
            _mask = ''
        ir.drizzle.in_mask = _mask
        ir.drizzle.outweig = _wgt
        ir.drizzle.xsh = xsh
        ir.drizzle.ysh = ysh
        ir.drizzle.outnx = imgsize
        ir.drizzle.outny = imgsize
        
        print('Drizzling: ', roots[i])
        print('     xsh = {0:8.2f}'.format( xsh ))
        print('     ysh = {0:8.2f}'.format( ysh ))
        print('  weight = {0:8.2f}'.format( weights[i] ))
        print('   outnx = {0:8d}'.format( imgsize ))
        ir.drizzle(_cdwt_ir, _tmpfits, Stdout=f_dlog)

        # Read .max file with saturation level for final combined image
        # by weighting each individual satLevel and summing.
        # Read in each satLevel from individual .max files
        _max = cleanDir + 'c' + roots[i] + '.max'
        #getsatLvl = Table.read(_max, format='ascii.no_header') #changed from , header_start=None
        #satLvl = getsatLvl[0][0]
        getsatLvl = open(_max)
        satLvl = float(getsatLvl.read())
        getsatLvl.close()
        satLvl_wt = satLvl * weights[i]
        satLvl_combo += satLvl_wt

    f_dlog.close()

    print('satLevel for combo image = ', satLvl_combo)
    # Write the combo saturation level to a file
    _max = open(_maxFile, 'w')
    _max.write('%15.4f' % satLvl_combo)
    _max.close()

    # Clean up the drizzled image of any largely negative values.
    # Don't do this! See how starfinder handles really negative pixels,
    # and if all goes well...don't ever correct negative pixels to zero.
    fits_f = fits.open(_tmpfits)
    
    tmp_stats = stats.sigma_clipped_stats(fits_f[0].data,
                                          sigma_upper=1, sigma_lower=10,
                                          iters=5)
    sci_mean = tmp_stats[0]
    sci_stddev = tmp_stats[2]

    # Find and fix really bad pixels
    idx = np.where(fits_f[0].data < (sci_mean - 10*sci_stddev))
    fits_f[0].data[idx] = sci_mean - 10*sci_stddev

    # Set the ROTPOSN value for the combined image.
    if (diffPA == 1):
        phi = 0.7
        fits_f[0].header.set('ROTPOSN', "%.5f" % phi,
                              'rotator user position')

    # Add keyword with distortion image information
    fits_f[0].header.set('DISTCOEF', "%s" % distCoef,
                          'Distortion Coefficients File')
    fits_f[0].header.set('DISTORTX', "%s" % distXgeoim,
                          'X Distortion Image')
    fits_f[0].header.set('DISTORTY', "%s" % distYgeoim,
                          'Y Distortion Image')

    # Fix the DATASEC header keyword, if it exists.
    if 'DATASEC' in fits_f[0].header:
        fits_f[0].header['DATASEC'] = '[1:{0:d},1:{0:d}]'.format(imgsize)
    
    # Calculate weighted MJD and store in header
    mjd_weightedMean = mjd_weightedSum / np.sum(weights)
    time_obs = Time(mjd_weightedMean, format='mjd')
    
    fits_f[0].header.set('MJD-OBS', mjd_weightedMean, 'Weighted modified julian date of combined observations')
    
    ## Also update date field in header
    fits_f[0].header.set('DATE', '{0}'.format(time_obs.fits), 'Weighted observation date')
    
    
    # Save to final fits file.
    fits_f[0].writeto(_fits, output_verify=outputVerify)
    util.rmall([_tmpfits, _cdwt])

def combine_submaps(imgsize, cleanDir, roots, outroot, weights,
            shifts, submaps, wave, diffPA, fixDAR=True, mask=True,
            instrument=instruments.default_inst):
    """
    Assumes the list of roots are pre-sorted based on quality. Images are then
          divided up with every Nth image going into the Nth submap.

    mask: (def=True) Set to false for maser mosaics since they have only
          one image at each positions. Masking produces artifacts that
          Starfinder can't deal with.
    """

    extend = ['_1', '_2', '_3']
    _out = [outroot + end for end in extend]
    _fits = [o + '.fits' for o in _out]
    _tmp = [o + '_tmp.fits' for o in _out]
    _wgt = [o + '_sig.fits' for o in _out]
    _log = [o + '_driz.log' for o in _out]
    _max = [o + '.max' for o in _out]

    util.rmall(_fits + _tmp + _wgt + _log + _max)

    # Prep drizzle stuff
    setup_drizzle(imgsize)
    print('Drizzle imgsize = ', imgsize)
    ir.drizzle.outcont = ''

    satLvl_tot = np.zeros(submaps, dtype=float)
    satLvl_sub = np.zeros(submaps, dtype=float)

    print('combine: drizzling sub-images together')
    f_log = [open(log, 'a') for log in _log]

    # Final normalization factor
    weightsTot = np.zeros(submaps, dtype=float)
    
    # Array to store weighted sum of MJDs in each submap
    mjd_weightedSums = np.zeros(submaps, dtype=float)

    # Get the distortion maps for this instrument.
    hdr0 = fits.getheader(cleanDir + 'c' + roots[0] + '.fits')
    distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr0)

    # Set a cleanDir variable in IRAF. This avoids the long-filename problem.
    ir.set(cleanDir=cleanDir)

    for i in range(len(roots)):
        # Cleaned image
        _c = cleanDir + 'c' + roots[i] + '.fits'
        _c_ir = _c.replace(cleanDir, 'cleanDir$')

        # Cleaned but distorted image
        _cd = cleanDir + 'distort/cd' + roots[i] + '.fits'
        cdwt = cleanDir + 'weight/cdwt.fits'
        _cd_ir = _cd.replace(cleanDir, 'cleanDir$')
        _cdwt_ir = cdwt.replace(cleanDir, 'cleanDir$')

        # Multiply each distorted image by it's weight
        util.rmall([cdwt])

        ir.imarith(_cd, '*', weights[i], _cdwt_ir)
        
        # Fix the ITIME header keyword so that it matches (weighted).
        # Drizzle will add all the ITIMEs together, just as it adds the flux.
        itime = fits.getval(cdwt, instrument.hdr_keys['itime'])
        itime *= weights[i]
        fits.setval(cdwt, instrument.hdr_keys['itime'], value=itime)
        
        # Get pixel shifts
        xsh = shifts[i][1]
        ysh = shifts[i][2]
        
        # Determine which submap we should be drizzling to.
        sub = int(i % submaps)
        fits_im = _tmp[sub]
        wgt = _wgt[sub]
        log = f_log[sub]
        
        # Read in PA of each file to feed into drizzle for rotation
        hdr = fits.getheader(_c,ignore_missing_end=True)
        phi = instrument.get_position_angle(hdr)
        if (diffPA == 1):
            ir.drizzle.rot = phi

        # Calculate saturation level for submaps
        # by weighting each individual satLevel and summing.
        # Read in each satLevel from individual .max files
        max_indiv = cleanDir + 'c' + roots[i] + '.max'
        satfile = open(max_indiv)
        satLvl = float(satfile.read()) #changed to simple i/o because the astropy table was breaking for a textfile with a single entry
        #getsatLvl = Table.read(max_indiv, format='ascii', header_start=None)
        #satLvl = getsatLvl[0][0]
        satLvl_wt = satLvl * weights[i]
        satLvl_tot[sub] += satLvl_wt
        
        # Add up the weights that go into each submap
        weightsTot[sub] += weights[i]
        
        satLvl_sub[sub] = satLvl_tot[sub] / weightsTot[sub]
        
        if (fixDAR == True):
            darRoot = cdwt.replace('.fits', 'geo')
            print('submap: ',cdwt)
            (xgeoim, ygeoim) = dar.darPlusDistortion(cdwt, darRoot,
                                                     xgeoim=distXgeoim,
                                                     ygeoim=distYgeoim,
                                                     instrument=instrument)
            xgeoim = xgeoim.replace(cleanDir, 'cleanDir$')
            ygeoim = ygeoim.replace(cleanDir, 'cleanDir$')
            ir.drizzle.xgeoim = xgeoim
            ir.drizzle.ygeoim = ygeoim
        else:
            ir.drizzle.xgeoim = distXgeoim
            ir.drizzle.ygeoim = distYgeoim

        # Read in MJD of current file from FITS header
        mjd = float(hdr['MJD-OBS'])
        mjd_weightedSums[sub] += weights[i] * mjd
        
        # Drizzle this file ontop of all previous ones.
        log.write(time.ctime())

        if (mask == True):
            _mask = 'cleanDir$masks/mask' + roots[i] + '.fits'
            #_mask = cleanDir + 'masks/mask' + roots[i] + '.fits'
        else:
            _mask = ''
        ir.drizzle.in_mask = _mask
        ir.drizzle.outweig = wgt
        ir.drizzle.xsh = xsh
        ir.drizzle.ysh = ysh

        ir.drizzle(_cdwt_ir, fits_im, Stdout=log)
    
    # Calculate weighted MJDs for each submap
    mjd_weightedMeans = mjd_weightedSums / weightsTot
    submaps_time_obs = Time(mjd_weightedMeans, format='mjd')
    
    for f in f_log:
        f.close()
        
    print('satLevel for submaps = ', satLvl_sub)
    # Write the saturation level for each submap to a file
    for l in range(submaps):
        _maxsub = open(_max[l], 'w')
        _maxsub.write('%15.4f' % satLvl_sub[l])
        _maxsub.close()

    for s in range(submaps):
        fits_f = fits.open(_tmp[s])
        
        # Clean up the drizzled image of any largely negative values.
        # Don't do this! See how starfinder handles really negative pixels,
        # and if all goes well...don't ever correct negative pixels to zero.
        tmp_stats = stats.sigma_clipped_stats(fits_f[0].data,
                                          sigma_upper=1, sigma_lower=10,
                                          iters=5)
        sci_mean = tmp_stats[0]
        sci_stddev = tmp_stats[2]

        # Find and fix really bad pixels
        idx = np.where(fits_f[0].data < (sci_mean - 10*sci_stddev))
        fits_f[0].data[idx] = 0.0

        # Normalize properly
        fits_f[0].data = fits_f[0].data / weightsTot[s]

        # Fix the ITIME header keyword so that it matches (weighted).
        itime = fits_f[0].header.get('ITIME')
        itime /= weightsTot[s]
        #fits_f[0].header.update('ITIME', '%.5f' % itime)
        fits_f[0].header['ITIME'] = ('%.5f' % itime)

        # Set the ROTPOSN value for the combined submaps.

        fits_f[0].header.set('ITIME', '%.5f' % itime)
        
        # Set the ROTPOSN value for the combined submaps. 

        if (diffPA == 1):
            phi = 0.7
            fits_f[0].header.set('ROTPOSN', "%.5f" % phi,
                                  'rotator user position')

        # Add keyword with distortion image information
        fits_f[0].header.set('DISTORTX', "%s" % distXgeoim,
                              'X Distortion Image')
        fits_f[0].header.set('DISTORTY', "%s" % distYgeoim,
                              'Y Distortion Image')
        
        # Store weighted MJDs in header
        fits_f[0].header.set('MJD-OBS', mjd_weightedMeans[s], 'Weighted modified julian date of combined observations')
    
        ## Also update date field in header
        fits_f[0].header.set('DATE', '{0}'.format(submaps_time_obs[s].fits), 'Weighted observation date')
        
        # Write out final submap fits file
        fits_f[0].writeto(_fits[s], output_verify=outputVerify)
    
    
    util.rmall(_tmp)
    util.rmall([cdwt])


def combine_rotation(cleanDir, roots, instrument=instruments.default_inst):
    """
    Determine if images are different PAs. If so, then
    temporarily rotate the images for xregister to use
    in order to get image shifts that are fed into drizzle.

    WARNING: If multiple PAs are found, then everything
    is rotated to PA = 0.
    """
    diffPA = 0

    clean_files = instrument.make_filenames(roots, rootDir=cleanDir, prefix='c')

    for cc in range(len(clean_files)):
        hdr = fits.getheader(clean_files[cc], ignore_missing_end=True)
        phi = instrument.get_position_angle(hdr)

        if cc == 0:
            phiRef = phi
            
        diff = phi - phiRef

        if (diff != 0.0):
            print('Different PAs found')
            diffPA = 1
            break

    if (diffPA == 1):
        for cc in range(len(clean_files)):
            hdr = fits.getheader(clean_files[cc], ignore_missing_end=True)
            phi = instrument.get_position_angle(hdr)
            rot_img(roots[cc], phi, cleanDir)

    return (diffPA)

def sort_frames(roots, strehls, fwhm, weights, shiftsTab):
    sidx = np.argsort(fwhm)

    # Make sorted lists.
    strehls = strehls[sidx]
    fwhm = fwhm[sidx]
    weights = weights[sidx]
    roots = [roots[i] for i in sidx]
    shiftsX = shiftsTab['col1']
    shiftsX = shiftsX[sidx]
    shiftsY = shiftsTab['col2']
    shiftsY = shiftsY[sidx]

    # Move all the ones with fwhm = -1 to the end
    gidx = (np.where(fwhm > 0))[0]
    bidx = (np.where(fwhm <= 0))[0]
    goodroots = [roots[i] for i in gidx]
    badroots = [roots[i] for i in bidx]
    if len(bidx) > 0:
        print('Found files with incorrect FWHM. They may be rejected.')
        print('\t' + ','.join(badroots))

    strehls = np.concatenate([strehls[gidx], strehls[bidx]])
    fwhm = np.concatenate([fwhm[gidx], fwhm[bidx]])
    weights = np.concatenate([weights[gidx], weights[bidx]])
    shiftsX = np.concatenate([shiftsX[gidx], shiftsX[bidx]])
    shiftsY = np.concatenate([shiftsY[gidx], shiftsY[bidx]])
    roots = goodroots + badroots

    newShiftsTab = shiftsTab.copy()
    for rr in range(len(newShiftsTab)):
        newShiftsTab[rr][0] = roots[rr]
        newShiftsTab[rr][1] = shiftsX[rr]
        newShiftsTab[rr][2] = shiftsY[rr]

    return (roots, strehls, fwhm, weights, newShiftsTab)


def combine_ref(coofile, cleanDir, roots, diffPA, refImage_index=0):
    """
    Pulls reference star coordinates from image header keywords.
    """
    # Delete any previously existing file
    util.rmall([coofile])

    cFits = [cleanDir + 'c' + r + '.fits' for r in roots]

    _allCoo = open(coofile, 'w')

    # write reference source coordinates
    hdr = fits.getheader(cFits[refImage_index],ignore_missing_end=True)
    _allCoo.write(' ' + hdr['XREF'] + '   ' + hdr['YREF'] + '\n')

    # write all coordinates, including reference frame
    for i in range(len(roots)):
        hdr = fits.getheader(cFits[i],ignore_missing_end=True)
        _allCoo.write(' ' + hdr['XREF'] + '   ' + hdr['YREF'] + '\n')

    _allCoo.close()


def combine_coo(coofile, cleanDir, roots, diffPA, refImage_index=0):
    """
    Pulls reference star coordinates from *.coo files.
    """
    # Delete any previously existing file
    util.rmall([coofile])

    # If images were rotated because of differing PAs, make a
    # different input list
    if (diffPA == 1):
        cCoos = [cleanDir + 'c' + r + '.rcoo' for r in roots]
    else:
        cCoos = [cleanDir + 'c' + r + '.coo' for r in roots]

    # Need to make table of coordinates of a reference source. These
    # will be used as initial estimates of the shifts (they don't necessarily
    # need to be real sources).
    _allCoo = open(coofile, 'w')

    # First line must be the coordinates in the reference image
    _allCoo.write(open(cCoos[refImage_index], 'r').read())

    # Now loop through all files (including the reference) and print
    # coordinates of same reference source.
    for i in range(len(roots)):
        _allCoo.write(open(cCoos[i], 'r').read())

    _allCoo.close()


def combine_lis(outfile, cleanDir, roots, diffPA):
    # Delete previously existing file
    util.rmall([outfile])

    cFits = [cleanDir + 'c' + r + '.fits' for r in roots]

    # Write all the files to a list
    f_lis = open(outfile, 'w')
    f_lis.write('\n'.join(cFits) + '\n')
    f_lis.close()

    # If images were rotated because of differing PAs, make a
    # different input list for xregister (to get shifts)
    if (diffPA == 1):
        rFits = [cleanDir + 'r' + r + '.fits' for r in roots]
        out = outfile + '_r'
        f_lis = open(out, 'w')
        f_lis.write('\n'.join(rFits) + '\n')
        f_lis.close()

def combine_register(outroot, refImage, diffPA):
    shiftFile = outroot + '.shifts'
    util.rmall([shiftFile])

    # xregister parameters
    ir.immatch
    ir.unlearn('xregister')
    ir.xregister.coords = outroot + '.coo'
    ir.xregister.output = ''
    ir.xregister.append = 'no'
    ir.xregister.databasefmt = 'no'
    ir.xregister.verbose = 'no'
    ir.xregister.xwindow='30'
    ir.xregister.ywindow='30'
    ir.xregister.correlation='fourier'
    ir.xregister.function='centroid'

    print('combine: registering images')
    if (diffPA == 1):
        input = '@' + outroot + '.lis_r'
    else:
        input = '@' + outroot + '.lis'

    hdu = fits.open(refImage)
    nx = hdu[0].header['NAXIS1']
    ny = hdu[0].header['NAXIS2']

    regions = '['+str(nx/2-nx/4)+':'+str(nx/2+nx/4)+','+str(ny/2-ny/4)+':'+str(ny/2+ny/4)+']'
    #regions = '[*,*]'
    # print 'input = ', input
    print('xwindow,ywindow',ir.xregister.xwindow,ir.xregister.ywindow)
    print('refImage = ', refImage)
    print('regions = ', regions)
    print('shiftFile = ', shiftFile)

    fileNames = Table.read(input[1:], format='ascii.no_header') # removed , header_start=None
    fileNames = np.array(fileNames)
    fileNames = np.array(fileNames, dtype='S')
    coords = Table.read(outroot + '.coo', format='ascii', header_start=None)
    shiftsTable_empty = np.zeros((len(fileNames), 3), dtype=float)
    shiftsTable = Table(shiftsTable_empty, dtype=('S50', float, float)) #dtype=(float, float, 'S50')

    for ii in range(len(fileNames)):
        inFile = fileNames[ii]

        tmpCooFile = outroot + '_tmp.coo'
        _coo = open(tmpCooFile, 'w')
        _coo.write('%.2f  %.2f\n' % (coords[0][0], coords[0][1]))
        _coo.write('%.2f  %.2f\n' % (coords[ii+1][0], coords[ii+1][1])) #Changed from [0][ii+1] to [ii+1][0]
        _coo.close()

        util.rmall([shiftFile])
        print('inFile = ', inFile)
        ir.xregister.coords = tmpCooFile
        ir.xregister(inFile, refImage, regions, shiftFile)

        # # Read in the shifts file. Column format is:
        # # Filename.fits  xshift  yshift
        _shifts = Table.read(shiftFile, format='ascii.no_header')
        shiftsTable[ii][0] = _shifts[0][0]
        shiftsTable[ii][1] = _shifts[0][1]
        shiftsTable[ii][2] = _shifts[0][2]


    util.rmall([shiftFile])
    shiftsTable.write(shiftFile, format = 'ascii')

    return (shiftsTable)


def combine_log(outroot, roots, strehls, fwhm, weights):
    _log = outroot + '.log'
    util.rmall([_log])

    f_log = open(_log, 'w')
    for i in range(len(roots)):
        f_log.write('c%s %6.2f %5.2f %6.3f\n' %
                        (roots[i], fwhm[i], strehls[i], weights[i]))

    f_log.close()

def combine_size(shiftsTable, refImage, outroot, subroot, submaps):
    """Determine the final size of the fully combined image. Use the
    shifts stored in the shiftsTable.

    @param shiftsTable: Table with x and y shifts for each image
    @type shiftsTable: ascii table
    @param refImage: The reference image from which the shifts are
        calculated from.
    @type refImage: string
    @param outroot: The name of the file for which shift information
        will be stored. The filename will be <outroot>.coo.
    @type outroot: string
    @param subroot: Same as outroot but for submaps
    @type subroot: string
    @param submaps: number of submaps
    @type sbumaps: int
    """
    x_allShifts = shiftsTable['col1']
    y_allShifts = shiftsTable['col2']


    xhi = abs(x_allShifts.max())
    xlo = abs(x_allShifts.min())
    yhi = abs(y_allShifts.max())
    ylo = abs(y_allShifts.min())

    # Make sure to include the edges of all images.
    # Might require some extra padding on one side.
    maxoffset = max([xlo, xhi, ylo, yhi])

    orig_img = fits.getdata(refImage)
    orig_size = (orig_img.shape)[0]
    padd = int(np.floor(orig_size * 0.01))

    # Read in the reference star's position in the ref image and translate
    # it into the coordinates of the final main and sub maps.
    hdr = fits.getheader(refImage, ignore_missing_end=True)
    xrefSrc = float(hdr['XREF'])
    yrefSrc = float(hdr['YREF'])

    xrefSrc = xrefSrc + (maxoffset + padd)
    yrefSrc = yrefSrc + (maxoffset + padd)

    cooMain = [outroot + '.coo']
    cooSubs = ['%s_%d.coo' % (subroot, i) for i in range(submaps+1)]
    cooAll = cooMain + cooSubs

    util.rmall(cooAll)
    for coo in cooAll:
        _allCoo = open(coo, 'w')
        _allCoo.write('%9.3f %9.3f\n' % (xrefSrc, yrefSrc))
        _allCoo.close()

    xysize = int(float(orig_size) + ((maxoffset + padd) * 2.0))
    print('combine: Size of output image is %d' % xysize)

    return xysize

def setup_drizzle(imgsize):
    """Setup drizzle parameters for NIRC2 data.
    @param imgsize: The size (in pixels) of the final drizzle image.
    This assumes that the image will be square.
    @type imgsize: int
    @param mask: The name of the mask to use during
    drizzle.
    @param type: str
    """
    # Setup the drizzle parameters we will use
    ir.module.load('stsdas', doprint=0, hush=1)
    ir.module.load('analysis', doprint=0, hush=1)
    ir.module.load('dither', doprint=0, hush=1)
    ir.unlearn('drizzle')
    ir.drizzle.outweig = ''
    ir.drizzle.in_mask = ''
    ir.drizzle.wt_scl = 1
    ir.drizzle.outnx = imgsize
    ir.drizzle.outny = imgsize
    ir.drizzle.pixfrac = 1
    ir.drizzle.kernel = 'lanczos3'
    ir.drizzle.scale = 1
    ir.drizzle.shft_un = 'input'
    ir.drizzle.shft_fr = 'output'
    ir.drizzle.align = 'center'
    ir.drizzle.expkey = 'ITIME'
    ir.drizzle.in_un = 'counts'
    ir.drizzle.out_un = 'counts'

def clean_drizzle(xgeoim, ygeoim, _bp, _cd, _wgt, _dlog, fixDAR=True, instrument=instruments.default_inst):
    # Get the distortion maps for this instrument.
    hdr = fits.getheader(_bp)
    distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr)
    
    if (fixDAR == True):
        darRoot = _cd.replace('.fits', 'geo')

        (xgeoim, ygeoim) = dar.darPlusDistortion(_bp, darRoot, xgeoim, ygeoim, instrument=instrument)

        ir.drizzle.xgeoim = xgeoim
        ir.drizzle.ygeoim = ygeoim
    else:
        ir.drizzle.xgeoim = distXgeoim
        ir.drizzle.ygeoim = distYgeoim

    ir.drizzle(_bp, _cd, outweig=_wgt, Stdout=_dlog)

def clean_cosmicrays(_ff, _mask, wave):
    """Clean the image of cosmicrays and make a mask containing the location
    of all the cosmicrays. The CR masks can later be used in combine() to
    keep cosmicrays from being included.

    @param _ff: Flat fielded file on which to fix cosmic rays. A new
        image will be created with the _f appended to it.
    @type _ff: string
    @param _mask: The filename used for the resulting mask.
    @type _mask: string
    @parram wave: The filter of the observations (e.g. 'kp', 'lp'). This
        is used to determine different thresholds for CR rejection.
    @type wave: string
    """
    # Determine the threshold at which we should start looking
    # for cosmicrays. Need to figure out the mean level of the
    # background.
    ff_img = fits.getdata(_ff)
    tmp_stats = stats.sigma_clipped_stats(ff_img,
                                          sigma_upper=2, sigma_lower=5,
                                          iters=5)
    mean = tmp_stats[0]
    stddev = tmp_stats[2]

    # CR candidates are those that exceed surrounding pixels by
    # this threshold amount.
    crthreshold = 5.0*stddev

    fluxray = 13.
    if 'h' in wave:
        fluxray = 10.
    if 'kp' in wave:
        fluxray = 13.
    if 'lp' in wave:
        fluxray = 10.0
    if 'ms' in wave:
        fluxray = 10.0

    ir.module.load('noao', doprint=0, hush=1)
    ir.module.load('imred', doprint=0, hush=1)
    ir.module.load('crutil', doprint=0, hush=1)
    ir.unlearn('cosmicrays')

    ir.cosmicrays(_ff, ' ', crmasks=_mask, thresho=crthreshold,
                  fluxrat=fluxray, npasses=10., window=7,
                  interac='no', train='no', answer='NO')

    ir.imcopy(_mask+'.pl', _mask, verbose='no')
    if os.path.exists(_mask + '.pl'): os.remove(_mask + '.pl')

def clean_cosmicrays2(_ff, _ff_cr, _mask, wave,
                      instrument=instruments.default_inst):
    """Clean the image of cosmicrays and make a mask containing the location
    of all the cosmicrays. The CR masks can later be used in combine() to
    keep cosmicrays from being included.

    @param _ff: Flat fielded file on which to fix cosmic rays. A new
        image will be created with the _f appended to it.
    @type _ff: string
    @param _ff_cr: Output image with cosmicrays fixed.
    @type _ff_cr: string
    @param _mask: The filename used for the resulting mask.
    @type _mask: string
    @parram wave: The filter of the observations (e.g. 'kp', 'lp'). This
        is used to determine different thresholds for CR rejection.
    @type wave: string
    """
    # Determine the threshold at which we should start looking
    # for cosmicrays. Need to figure out the mean level of the
    # background.
    ff_img, ff_hdr = fits.getdata(_ff, header=True)
    mean, median, stddev = stats.sigma_clipped_stats(ff_img,
                                                     sigma_upper=2, sigma_lower=5,
                                                     iters=5)

    # Get the instrument gain
    gain = instrument.get_gain()

    sampmode = ff_hdr[instrument.hdr_keys['sampmode']]
    numreads = ff_hdr[instrument.hdr_keys['nfowler']]
    
    if sampmode == 2:
        readnoise = 60
    else:
        readnoise = 15.0 * (16.0 / numreads)**0.5

    from jlu.util import cosmics
    c = cosmics.cosmicsimage(ff_img, gain=gain, readnoise=readnoise,
                             sigclip=10, sigfrac=0.5, objlim=5.0)
    c.run(maxiter=3)
    fits.writeto(_ff_cr, c.cleanarray, ff_hdr,
                   clobber=True, output_verify=outputVerify)
    fits.writeto(_mask, np.where(c.mask==True, 1, 0), ff_hdr,
                   clobber=True, output_verify=outputVerify)

    return

def clean_persistance(_n, _pers, instrument=instruments.default_inst):
    """
    Make masks of the persistance to be used in combining the images
    later on.
    """
    # Read in image
    fits_f = fits.open(_n)
    img = fits_f[0].data

    # Define the high pixels
    persPixels = where(img > instrument.get_saturation_level())

    # Set saturated pixels to 0, good pixels to 1
    fits_f[0].data[persPixels] = 0
    fits_f[0].data = fits_f[0].data / fits_f[0].data

    # Save to an image
    fits_f[0].writeto(_pers, output_verify=outputVerify)


def clean_bkgsubtract(_ff_f, _bp):
    """Do additional background subtraction of any excess background
    flux. This isn't strictly necessary since it just removes a constant."""
    # Open the image for processing.
    fits_f = fits.open(_ff_f)

    # Calculate mean and STD for science image
    tmp_stats = stats.sigma_clipped_stats(fits_f[0].data,
                                          sigma_upper=1, sigma_lower=10,
                                          iters=5)
    sci_mean = tmp_stats[0]
    sci_stddev = tmp_stats[2]

    # Excess background flux at (mean - 2*std)
    bkg = sci_mean - (2.0 * sci_stddev)
    #print 'Bkg mean = %5d +/- %5d   bkg = %5d  Name = %s' % \
    #      (sci_mean, sci_stddev, bkg, _ff_f)

    # Open old, subtract BKG

    # Find really bad pixels
    idx = np.where(fits_f[0].data < (sci_mean - 10*sci_stddev))

    # Subtract background
    fits_f[0].data -= bkg

    # Fix really bad negative pixels.
    fits_f[0].data[idx] = 0.0

    # Write to new file
    fits_f[0].writeto(_bp, output_verify=outputVerify)

    # Return the background we subtracted off
    return bkg

def clean_makecoo(_ce, _cc, refSrc, strSrc, aotsxyRef, radecRef,
                  instrument=instruments.default_inst, check_loc=True,
                  update_fits=True,cent_box=12):
    """Make the *.coo file for this science image. Use the difference
    between the AOTSX/Y keywords from a reference image and each science
    image to tell how the positions of the two frames are related.

    @param _ce: Name of the input cleaned file.
    @type _ce: string
    @param _cc: Name of the output header modified image.
    @type _cc: string
    @param refSrc: Array with the X/Y positions of the reference source.
        This will be put into the image header and the *.coo file.
    @type refSrc: array of floats with length=2 [x, y]
    @param strSrc: Array with the X/Y positions of the strehl source.
        This will be put into the image header.
    @type strSrc: array of floats with length=2 [x, y]
    @param aotsxyRef: The AOTSX/Y header values from the reference image.
    @type aotsxyRef: array of floats with length=2 [x, y]
    @param radecRef: The RA/DEC header values from the reference image.
    @type radecRef: array of floats with length=2 [x, y]

    check_loc (bool):  If True the reference source is recentered for this frame.
                     Use False if the offsets are large enough to move the reference source off of the image
    update_fits : update the fits files with the reference pixel values
    cent_box : box size to center the source (default: 12)
    """

    hdr = fits.getheader(_ce, ignore_missing_end=True)

    radec = [float(hdr['RA']), float(hdr['DEC'])]
    aotsxy = kai_util.getAotsxy(hdr)

    # Determine the image's PA and plate scale
    phi = instrument.get_position_angle(hdr)
    scale = instrument.get_plate_scale(hdr)

    # Determine the instrument angle w.r.t. the AO bench.
    inst_angle = instrument.get_instrument_angle(hdr)

    # Calculate the pixel offsets from the reference image
    # We've been using aotsxy2pix, but the keywords are wrong
    # for 07maylgs and 07junlgs
    #d_xy = kai_util.radec2pix(radec, phi, scale, radecRef)
    d_xy = kai_util.aotsxy2pix(aotsxy, scale, aotsxyRef, inst_angle=inst_angle)

    # In the new image, find the REF and STRL coords
    xref = refSrc[0] + d_xy[0]
    yref = refSrc[1] + d_xy[1]
    xstr = strSrc[0] + d_xy[0]
    ystr = strSrc[1] + d_xy[1]
    print('clean_makecoo: xref, yref start = {0:.2f} {1:.2f}'.format(xref, yref))

    # re-center stars to get exact coordinates
    if check_loc:

        text = ir.imcntr(_ce, xref, yref, cbox=cent_box, Stdout=1)
        values = text[0].split()
        xref = float(values[2])
        yref = float(values[4])

        text = ir.imcntr(_ce, xstr, ystr, cbox=cent_box, Stdout=1)
        values = text[0].split()
        xstr = float(values[2])
        ystr = float(values[4])
        print('clean_makecoo: xref, yref final = {0:.2f} {1:.2f}'.format(xref, yref))

    # write reference star x,y to fits header
    if update_fits:
        fits_f = fits.open(_ce)
        fits_f[0].header.set('XREF', "%.3f" % xref,
                             'Cross Corr Reference Src x')
        fits_f[0].header.set('YREF', "%.3f" % yref,
                             'Cross Corr Reference Src y')
        fits_f[0].header.set('XSTREHL', "%.3f" % xstr,
                             'Strehl Reference Src x')
        fits_f[0].header.set('YSTREHL', "%.3f" % ystr,
                             'Strehl Reference Src y')
        fits_f[0].writeto(_cc, output_verify=outputVerify)

    file(_cc.replace('.fits', '.coo'), 'w').write('%7.2f  %7.2f\n' % (xref, yref))

    # Make a temporary rotated coo file, in case there are any data sets
    # with various PAs; needed for xregister; remove later
    xyRef_rot = kai_util.rotate_coo(xref, yref, phi)
    xref_r = xyRef_rot[0]
    yref_r = xyRef_rot[1]

    xyStr_rot = kai_util.rotate_coo(xstr, ystr, phi)
    xstr_r = xyStr_rot[0]
    ystr_r = xyStr_rot[1]

    file(_cc.replace('.fits', '.rcoo'), 'w').write('%7.2f  %7.2f\n' % (xref_r, yref_r))

    return

def mosaic_ref(outFile, cleanDir, roots, diffPA, instrument=instruments.default_inst):
    """Calculate an initial guess at the offsets between mosaic frames.
    using the AOTSX/Y keywords from a reference image and each science
    image to tell how the positions of the two frames are related.

    @param cleanDir: Name of the input cleaned file.
    @type cleanDir: string
    @param roots: List of root filenames
    @type roots: list of strings
    @param diffPA: 1 = found different PAs so use rot images.
    @type difPA: int
    """
    
    # Setup clean file lists.
    if (diffPA == 1):
        fileNames = instrument.make_filenames(roots, rootDir=cleanDir, prefix='r')
    else:
        fileNames = instrument.make_filenames(roots, rootDir=cleanDir, prefix='c')
        
    hdrRef = fits.getheader(fileNames[0], ignore_missing_end=True)
    aotsxyRef = kai_util.getAotsxy(hdrRef)

    # Determine the image's PA and plate scale
    phi = instrument.get_position_angle(hdrRef)
    scale = instrument.get_plate_scale(hdrRef)
    inst_angle = instrument.get_instrument_angle(hdrRef)
    print('inst_angle = ', inst_angle)

    _out = open(outFile, 'w')

    # First line of shifts file must be for a reference
    # image (assumed to be the first image).
    _out.write('%7.2f  %7.2f\n' % (0.0, 0.0))

    for rr in range(len(roots)):
        hdr = fits.getheader(fileNames[rr], ignore_missing_end=True)
        aotsxy = kai_util.getAotsxy(hdr)

        # Calculate the pixel offsets from the reference image
        # We've been using aotsxy2pix, but the keywords are wrong
        # for 07maylgs and 07junlgs
        d_xy = kai_util.aotsxy2pix(aotsxy, scale, aotsxyRef, inst_angle=inst_angle)

        _out.write('%7.2f  %7.2f\n' % (d_xy[0], d_xy[1]))

    _out.close()

    return

class Sky(object):
    def __init__(self, sciDir, skyDir, wave, scale=1,
                 skyfile='', angleOffset=0.0,
                 instrument=instruments.default_inst):
        # Setup some variables we will need later on
        self.sciDir = sciDir
        self.skyDir = skyDir
        self.wave = wave
        self.skyFile = skyfile
        self.scale = scale
        self.angleOffset = angleOffset

        self.instrument = instrument

        self.defaultSky = skyDir + 'sky_' + wave + '.fits'

        if (wave == 'lp' or wave == 'ms'):
            self.__initLp__()

        # This will be the final returned skyname
        self.skyName = skyDir + 'sky_scaled.fits'

    def __initLp__(self):
        print('Initializing Lp Sky skyfile=%s' % (self.skyFile))

        # Read skies from manual sky file (format: raw_science   sky)
        if (self.skyFile):
            skyTab = Table.read(self.skyDir + self.skyFile,
                                format='ascii', header_start=None)
            self.images = skyTab[skyTab.colnames[0]]
            skies = skyTab[skyTab.colnames[1]]

            skyAng = np.zeros([len(skies)], Float64)
            for i in range(0,len(skies)):
                sky = skies[i].strip()
                hdr = fits.getheader(self.skyDir + sky, ignore_missing_end=True)
                skyAng[i] = float(hdr['ROTPPOSN'])

        else:
            # Read in the sky table. Determine the effective K-mirror
            # angle for each sky.
            skyTab = Table.read(self.skyDir + 'rotpposn.txt',
                                format='ascii', header_start=None)
            skies = skyTab[skyTab.colnames[0]]
            skyAng = skyTab[skyTab.colnames[1]]

        # The optimal sky angle to use is skyAng = A + B*sciAng
        self.angFitA = self.angleOffset
        self.angFitB = 1.0

        # Open a log file that we will keep
        _skylog = self.sciDir + 'sci_sky_subtract.log'
        util.rmall([_skylog])
        f_skylog = open(_skylog, 'w')

        # Stuff we are keeping
        self.skyTab = skyTab
        self.skies = skies
        self.skyAng = skyAng
        self.f_skylog = f_skylog

    def getSky(self, _n):
        if (self.wave == 'lp' or self.wave == 'ms'):
            sky = self.getSkyLp(_n)
        else:
            sky = self.defaultSky

        # Edit the science image to contain the
        # original sky name that will be subtracted.
        skyOrigName = sky[sky.rfind('/')+1:]
        ir.hedit(_n, 'SKYSUB', skyOrigName, add='yes', show='no', verify='no')

        # Now scale the sky to the science image
        skyScale = self.scaleSky(_n, sky)

        return skyScale

    def scaleSky(self, _n, _sky):
        """Scale the mean level of the sky so that it matches the
        science image.

        @param _n: name of science frame
        @type _n: string
        @param _sky: name of sky frame
        @type _sky: string
        """
        util.rmall([self.skyName])

        # scale sky to science frame
        if self.scale:
            n_img = fits.getdata(_n)
            sci_stats = stats.sigma_clipped_stats(n_img,
                                                  sigma_upper=1, sigma_lower=10,
                                                  iters=20)
            sci_mean = sci_stats[0]

            sky_img = fits.getdata(_sky)
            sky_stats = stats.sigma_clipped_stats(sky_img,
                                                  sigma_upper=5, sigma_lower=15,
                                                  iters=5)
            sky_mean = sky_stats[0]
            
            
            fact = sci_mean/sky_mean
            #print 'scaleSky: factor = %5f  sci_mean = %5f  sky_mean = %5f' % \
            #      (fact, sci_mean, sky_mean)
            ir.imarith(_sky, '*', fact, self.skyName)
        else:
            ir.imcopy(_sky, self.skyName)

        return self.skyName


    def getSkyLp(self, _n):
        """Determine which sky we should use for L'. Does all the
        rotator mirror angle matching.

        @param _n: Name of science frame.
        @type _n: string
        @returns sky: name of sky file to use.
        @rtype sky: string
        """
        # Sky subtract
        # determine the best angle for sky or use manual file

        # -- Determine the rotpposn for this image
        sciAng_tmp = ir.hselect(_n, "ROTPPOSN", "yes", Stdout=1)
        sciAng = float(sciAng_tmp[0])

        # -- Determine the best sky rotpposn.
        skyBest = self.angFitA + (self.angFitB * sciAng)

        # -- Repair all angles to be between -180 and 180.
        if (skyBest > 180): skyBest -= 360.0
        if (skyBest < -180): skyBest += 360.0
        if (sciAng > 180): sciAng -= 360.0
        if (sciAng < -180): sciAng += 360.0

        if (self.skyFile):
            for i in range(0,len(self.images)):
                if (self.images[i] == _n):
                    skyidx = i
        else:
            # -- Determine which sky file to use
            diff = [abs(skyAngle - skyBest) for skyAngle in self.skyAng]
            skyidx = np.argmin(diff)

        sky = self.skyDir + self.skies[skyidx]

        print(('Science = ', _n))
        print(('Sky image = ', sky))

        foo = '%s - %s  %6.1f  %6.1f' % \
              (_n, self.skies[skyidx], sciAng, self.skyAng[skyidx])

        self.f_skylog.write( foo )

        return sky

    def getNonlinearCorrection(self, sky):
        """Determine the non-linearity level. Raw data level of
        non-linearity is 12,000 but we subtracted
        off a sky which changed this level. The sky is
        scaled, so the level will be slightly different
        for every frame.

        @param sky: File name of the sky used.
        @type sky: string
        @returns (sky_mean + sky_stddev) which is the value that should
            be subtracted off of the saturation count level.
        @rtype float
        """
        # Read in the FITS file
        sky_img, sky_hdr = fits.getdata(sky, header=True)

        # Get the sigma-clipped mean and stddev
        sky_stats = stats.sigma_clipped_stats(sky_img,
                                              sigma_upper=4, sigma_lower=4,
                                              iters=4)
        sky_mean = sky_stats[0]
        sky_stddev = sky_stats[2]

        # -- Log what we did
        if (self.wave == 'lp' or self.wave == 'ms'):
            foo = ' %7d %7d\n' % (sky_mean, sky_stddev)

            self.f_skylog.write( foo )

        return sky_mean + sky_stddev

    def close(self):
        """Close log files opened at init."""
        if (self.wave == 'lp' or self.wave == 'ms'):
            self.f_skylog.close()


def mosaic(files, wave, outroot, field=None, outSuffix=None,
            trim=0, weight=0, fwhm_max=0, submaps=0, fixDAR=True, maskSubmap=False,
            instrument=instruments.default_inst):
    """Accepts a list of cleaned images and does a weighted combining after
    performing frame selection based on the Strehl and FWHM.

    Each image must have an associated *.coo file which gives the rough
    position of the reference source.

    @param files: List of integer file numbers to include in combine.
    @type files: list of int
    @param wave: Filter of observations (e.g. 'kp', 'lp', 'h')
    @type wave: string
    @param outroot: The output root name (e.g. '06jullgs'). The final combined
        file names will be <outroot>_<field>_<wave>. The <field> keyword
        is optional.

        Examples:
        06jullgs_kp for outroot='06jullgs' and wave='kp'
        06jullgs_arch_f1_kp for adding field='arch_f1'
    @type outroot: string
    @kwparam field: Optional field name used to get to clean directory and
        also effects the final output file name.
    @type field: string
    @kwparam trim: Optional file trimming based on image quality. Default
        is 0. Set to 1 to turn trimming on.
    @kwparam outSuffix: Optional suffix used to modify final output file name.
    @type outSuffix: string
    @type trim: 0 or 1
    @kwparam weight: Optional weighting based on Strehl. Set to 1 to
        to turn file weighting on (default is 0).
    @type weight: 0 or 1
    @kwparam fwhm_max: The maximum allowed FWHM for keeping frames when
        trimming is turned on.
    @type fwhm_max: int
    @kwparam submaps: Set to the number of submaps to be made (def=0).
    @type submaps: int
    @kwparam mask: Set to false for maser mosaics; 06maylgs1 is an exception
    @type mask: Boolean
    """
    # Start out in something like '06maylgs1/reduce/kp/'
    # Setup some files and directories
    waveDir = util.getcwd()
    redDir = util.trimdir( os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir( os.path.abspath(redDir + '../') + '/')
    if (field != None):
        cleanDir = util.trimdir( os.path.abspath(rootDir +
                                                   'clean/' +field+
                                                   '_' +wave) + '/')
        outroot += '_' + field
    else:
        cleanDir = util.trimdir( os.path.abspath(rootDir +
                                                   'clean/' + wave) + '/')

    if (outSuffix != None):
        outroot += outSuffix

    # This is the final output directory
    comboDir = rootDir + 'combo/'
    util.mkdir(comboDir)

    # Make strings out of all the filename roots.
    roots = instrument.make_filenames(files, prefix='')

    # This is the output root filename
    _out = comboDir + 'mag' + outroot + '_' + wave
    _sub = comboDir + 'm' + outroot + '_' + wave

    ##########
    # Determine if we are going to trim and/or weight the files
    # when combining. If so, then we need to determine the Strehl
    # and FWHM for each image. We check strehl source which shouldn't
    # be saturated. *** Hard coded to strehl source ***
    ##########

    # Load Strehls and FWHM for sorting and trimming
    strehls, fwhm = loadStrehl(cleanDir, roots)

    # Default weights
    # Create an array with length equal to number of frames used,
    # and with all elements equal to 1/(# of files)
    weights = np.array( [1.0/len(roots)] * len(roots) )

    ##########
    # Trimming
    ##########
    if trim:
        roots, strehls, fwhm, weights = trim_on_fwhm(roots, strehls, fwhm,
                                                     fwhm_max=fwhm_max)

    ##########
    # Weighting
    ##########
    if weight == 'strehl':
        weights = weight_by_strehl(roots, strehls)

    if ((weight != None) and (weight != 'strehl')):
        # Assume weight is set to a filename
        if not os.path.exists(weight):
            raise ValueError('Weights file does not exist, %s' % weight)

        weights = readWeightsFile(roots, weight)

    # Determine the reference image
    refImage = cleanDir + 'c' + roots[0] + '.fits'
    print('combine: reference image - %s' % refImage)

    ##########
    # Write out a log file. With a list of images in the
    # final combination.
    ##########
    combine_log(_out, roots, strehls, fwhm, weights)

    # See if all images are at same PA, if not, rotate all to PA = 0
    # temporarily. This needs to be done to get correct shifts.
    print('Calling combine_rotation')
    diffPA = combine_rotation(cleanDir, roots, instrument=instrument)

    # Make a table of initial guesses for the shifts.
    # Use the header keywords AOTSX and AOTSY to get shifts.
    print('Calling mosaic_ref')
    mosaic_ref(_out + '.init.shifts', cleanDir, roots, diffPA, instrument=instrument)

    # Keep record of files that went into this combine
    print('Calling combine_lis')
    combine_lis(_out + '.lis', cleanDir, roots, diffPA)

    # Register images to get shifts.
    print('Calling mosaic_register')
    shiftsTab = mosaic_register(_out, refImage, diffPA)

    # Determine the size of the output image from max shifts
    print('Calling mosaic_size')
    xysize = mosaic_size(shiftsTab, refImage, _out, _sub, submaps)

    # Combine all the images together.
    print('Calling mosaic_drizzle')
    combine_drizzle(xysize, cleanDir, roots, _out, weights, shiftsTab,
                    wave, diffPA, fixDAR=fixDAR, instrument=instrument)

    # Now make submaps
    if (submaps > 0):
        combine_submaps(xysize, cleanDir, roots, _sub, weights,
                        shiftsTab, submaps, wave, diffPA,
                        fixDAR=fixDAR, mask=maskSubmap)

    # Remove *.lis_r file & rotated rcoo files, if any - these
    # were just needed to get the proper shifts for xregister
    _lisr = _out + '.lis_r'
    util.rmall([_lisr])
    for i in range(len(roots)):
        _rcoo = cleanDir + 'c' + str(roots[i]) + '.rcoo'
        util.rmall([_rcoo])


def mosaic_register(outroot, refImage, diffPA):
    """
    Register images for a mosaic. This only calculates the exact
    shifts between each image... it doesn't do the combining.

    @param outroot: The root for the output image. The resulting
    shifts will be written into a file called <outroot>.shifts
    @type outroot: string
    @param refImage: The name of the reference image.
    @type refImage: string
    """
    shiftFile = outroot + '.shifts'
    util.rmall([shiftFile])

    # xregister parameters
    ir.immatch
    ir.unlearn('xregister')
    ir.xregister.coords = outroot + '.init.shifts'
    ir.xregister.output = ''
    ir.xregister.append = 'no'
    ir.xregister.databasefmt = 'no'
    ir.xregister.verbose = 'yes'
    ir.xregister.correlation = 'fourier'
    ir.xregister.xwindow = '10'
    ir.xregister.ywindow = '10'

    print('combine: registering images')
    if (diffPA == 1):
        input = '@' + outroot + '.lis_r'
    else:
        input = '@' + outroot + '.lis'

    regions = '[*,*]'
    ir.xregister(input, refImage, regions, shiftFile)

    # Read in the shifts file. Column format is:
    # Filename.fits  xshift  yshift
    shiftsTable = Table.read(shiftFile, format='ascii', header_start=None)

    return (shiftsTable)


def mosaic_size(shiftsTable, refImage, outroot, subroot, submaps):
    """
    Determine the final size for the completed mosaic.

    @params shiftsTable: Table from mosaic_register containing the
    shifts for all the images.
    @type shiftsTable: string
    @param refImage: The first image used as  reference.
    @type refImage: string
    @param outroot: The root name for the resulting output file.
    @type outroot: string
    @param subroot:
    @type subroot: string
    @param submaps:
    @type submaps:
    """
    x_allShifts = shiftsTable['col1']
    y_allShifts = shiftsTable['col2']

    xhi = abs(x_allShifts.max())
    xlo = abs(x_allShifts.min())
    yhi = abs(y_allShifts.max())
    ylo = abs(y_allShifts.min())

    # Make sure to include the edges of all images.
    # Might require some extra padding on one side.
    maxoffset = max([xlo, xhi, ylo, yhi])

    orig_img = fits.getdata(refImage)
    orig_size = (orig_img.shape)[0]
    padd = int(np.floor(orig_size * 0.02))

    xref = x_allShifts[0]
    yref = y_allShifts[0]

    xref = xref + (maxoffset + padd)
    yref = yref + (maxoffset + padd)

    # Read in 16C's position in the ref image and translate
    # it into the coordinates of the final main and sub maps.
    hdr = fits.getheader(refImage,ignore_missing_end=True)
    xrefSrc = float(hdr['XREF'])
    yrefSrc = float(hdr['YREF'])

    xrefSrc = xrefSrc + (maxoffset + padd)
    yrefSrc = yrefSrc + (maxoffset + padd)

    cooMain = [outroot + '.coo']
    cooSubs = ['%s_%d.coo' % (subroot, i) for i in range(submaps+1)]
    cooAll = cooMain + cooSubs

    util.rmall(cooAll)
    for coo in cooAll:
        _allCoo = open(coo, 'w')
        _allCoo.write('%9.3f %9.3f\n' % (xrefSrc, yrefSrc))
        _allCoo.close()

    xysize = float(orig_size) + ((maxoffset + padd) * 2.0)
    print('combine: Size of output image is %d' % xysize)

    return xysize
