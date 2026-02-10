import os, sys
import pylab as plt
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import stats
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy.nddata import CCDData
from astropy.nddata import block_replicate
from astropy import units as u
from astropy.coordinates import SkyCoord
import ccdproc as ccdp
import math
import drizzle
import kai
from . import kai_util
from kai.reduce import util, lin_correction
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
import warnings
from datetime import datetime
from scipy.ndimage import shift
from scipy.optimize import least_squares
from scipy.ndimage import rotate
from scipy.fft import fftn, ifftn, fftshift, ifftshift
from scipy.signal import windows

from scipy import ndimage
from scipy.interpolate import griddata


module_dir = os.path.dirname(__file__)

supermaskName = 'supermask.fits'
outputVerify = 'ignore'

def clean(files, nite, wave, refSrc, strSrc,
        dark_frame=None,
        badColumns=None, field=None,
        skyscale=False, skyfile=None, angOff=0.0, cent_box=12,
        fixDAR=True, use_koa_weather=False,
        raw_dir=None, clean_dir=None,
        instrument=instruments.default_inst, check_ref_loc=True,
        ref_offset_method='aotsxy'):
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
    refSrc : [float, float]
        x and y coordinates for the reference source, provided as a list of two
        float coordinates.
    strSrc : [float, float]
        x and y coordinates for the Strehl source, provided as a list of two
        float coordinates.
    dark_frame : str, default=None
        File name for the dark frame in order to carry out dark correction.
        If not provided, dark frame is not subtracted and a warning is thrown.
        Assumes dark file is located under ./calib/darks/
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
    cent_box : int (def = 12)
        the box to use for better centroiding the reference star
    badColumns : int array, default = None
        An array specifying the bad columns (zero-based).
        Assumes a repeating pattern every 8 columns.
    fixDAR : boolean, default = True
        Whether or not to calculate DAR correction coefficients.
    use_koa_weather : boolean, default = False
        If calculating DAR correction, this keyword specifies if the atmosphere
        conditions should be downloaded from the KOA weather data. If False,
        atmosphere conditions are downloaded from the MKWC CFHT data.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    clean_dir : str, optional
        Directory where clean files will be stored. By default,
        assumes that clean files will be stored in '../clean'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    ref_offset_method : str, default='aotsxy'
        Method to calculate offsets from reference image.
        Options are 'aotsxy' or 'radec'.
        In images where 'aotsxy' keywords aren't reliable, 'radec' calculated
        offsets may work better.
    """

    # Determine directory locations
    redDir = os.getcwd() + '/'
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')

    # Set location of raw data
    rawDir = rootDir + 'raw/'
    # Check if user has specified a specific raw directory
    if raw_dir is not None:
        if raw_dir.startswith('/'):
            rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')
        else:
            rawDir = util.trimdir(os.path.abspath(redDir + raw_dir) + '/')

    waveDir = util.trimdir(os.path.abspath(redDir + wave) + '/')
    sciDir = util.trimdir(os.path.abspath(waveDir + '/sci_' + nite) + '/')

    # Make sure directory for current passband exists and switch into it
    util.mkdir(wave)
    os.chdir(wave)
    
    util.mkdir(sciDir)
    os.chdir(sciDir)

    # Setup the clean directory
    cleanRoot = rootDir + 'clean/'
    # Check if user has specified a specific clean directory
    if clean_dir is not None:
        if clean_dir.startswith('/'):
            cleanRoot = util.trimdir(os.path.abspath(clean_dir) + '/')
        else:
            cleanRoot = util.trimdir(os.path.abspath(redDir + clean_dir) + '/')

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
        radecRef = instrument.get_radec(hdr1)
        aotsxyRef = kai_util.getAotsxy(hdr1)

        if ref_offset_method == 'pcu':
            pcuxyRef = instrument.get_pcuxyRef(hdr1)
        else:
            pcuxyRef = None

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
        
        # Read in dark frame data
        # Throw warning if dark frame not provided for dark correction
        if dark_frame is not None:
            dark_file = redDir + '/calib/darks/' + dark_frame
            
            # Read in dark frame data
            dark_data = fits.getdata(dark_file, ignore_missing_end=True)
        else:
            warning_message = 'Dark frame not provided for clean().'
            warning_message += '\nCleaning without dark subtraction.'
        
            warnings.warn(warning_message)
        
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
            util.rmall([
                _cp, _ss, _ff, _ff_f, _ff_s, _bp, _cd, _ce, _cc,
                _wgt, _statmask, _crmask, _mask, _pers, _max, _coo,
                _rcoo, _dlog,
            ])

            ### Copy the raw file to local directory ###
            if os.path.exists(_cp): os.remove(_cp)
            # Save as the primary HDU since Osiris saves as SCI by default
            raw_file = fits.open(_raw, ignore_missing_end=True)
            cp_primary_hdu = fits.PrimaryHDU(data=raw_file[0].data, header=raw_file[0].header)
            cp_primary_hdu.name = 'PRIMARY'
            cp_primary_hdu.header.pop('EXTNAME', None)
            cp_primary_hdu.header.pop('EXTVER', None)
            cp_primary_hdu.writeto(_cp, output_verify='ignore')
            #shutil.copy(_raw, _cp)

            # FIXME Add KAI version to header
            #with fits.open(_cp, mode="update") as filehandle:
            #    filehandle[0].header['ORIGIN'] = 'KAI v' + kai.__version__
            
            

            ### Make persistance mask ###
            # - Checked images, this doesn't appear to be a large effect.
            #clean_persistance(_cp, _pers, instrument=instrument)
            
            # Dark correction
            if dark_frame is not None:
                with fits.open(_cp, mode='denywrite', output_verify = 'ignore', 
                ignore_missing_end=True) as cur_frame:
                    frame_data = cur_frame[0].data
                    frame_header = cur_frame[0].header
                frame_data = frame_data - dark_data
                frame_hdu = fits.PrimaryHDU(data=frame_data, 
                header=frame_header)
                frame_hdu.writeto(_cp, 
                output_verify='ignore', 
                overwrite=True)
            
            # Linearity correction
            if instrument == 'NIRC2':
                lin_correction.lin_correction(_cp, instrument=instrument)
            
            ### Sky subtract ###
            # Get the proper sky for this science frame.
            # It might be scaled or there might be a specific one for L'.
            sky = skyObj.getSky(_cp)

            util.imarith(_cp, '-', sky, _ss)
            
            # Check if sky subtraction is correct 
            # or if scale sky must be applied
            ss = fits.getdata(_ss)
            if np.median(ss) < -10:
                raise Exception('Sky subtraction caused negative image. Rerun clean() with skyscale = True')

            ### Flat field ###
            util.imarith(_ss, '/', flat, _ff)
            
            ### Make a static bad pixel mask ###
            # _statmask = supermask + bad columns
            clean_get_supermask(_statmask, _supermask, badColumns)

            ### Fix bad pixels ###
            # Produces _ff_f file
            bfixpix.bfixpix(_ff, _statmask)
            util.rmall([_ff_s])

            ### Fix cosmic rays and make cosmic ray mask. ###
            clean_cosmicrays(_ff_f, _crmask, wave, _supermask)

            ### Combine static and cosmic ray mask ###
            # This will be used in combine later on.
            # Results are stored in _mask, _mask_static is deleted.
            clean_makemask(_mask, _crmask, _statmask, wave, instrument=instrument)

            ### Background Subtraction ###
            bkg = clean_bkgsubtract(_ff_f, _bp)

            ### Drizzle individual file ###
            clean_drizzle(distXgeoim, distYgeoim, _bp, _ce, _wgt, _dlog,
                          fixDAR=fixDAR, instrument=instrument,
                          use_koa_weather=use_koa_weather)
            
            hdr = fits.getheader(_raw, ignore_missing_end=True)
            
            ### Make .max file ###
            # Determine the non-linearity level. Raw data level of
            # non-linearity is 12,000 but we subtracted
            # off a sky which changed this level. The sky is
            # scaled, so the level will be slightly different
            # for every frame.
            nonlinSky = skyObj.getNonlinearCorrection(sky)

            coadds = fits.getval(_ss, instrument.hdr_keys['coadds'])
            sat_level_units = instrument.get_saturation_level_units()
            if sat_level_units == 'DN':
                satLevel = (coadds*instrument.get_saturation_level(hdr)) - nonlinSky - bkg
                open(_max, 'w').write(str(satLevel))
            elif sat_level_units == 'DN/coadd':
                satLevel = (instrument.get_saturation_level(hdr)) - nonlinSky - bkg
                open(_max, 'w').write(str(satLevel))

            ### Rename and clean up files ###
            shutil.move(_bp, _cd)
            # util.rmall([_cp, _ss, _ff, _ff_f])

            ### Make the *.coo file and update headers ###
            # First check if PA is not zero
            phi = instrument.get_position_angle(hdr)

            clean_makecoo(_ce, _cc, refSrc, strSrc, aotsxyRef, radecRef,
                          instrument=instrument, check_loc=check_ref_loc,
                          cent_box=cent_box, offset_method=ref_offset_method,pcuxyRef=pcuxyRef)

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
        os.chdir('../')
        os.chdir(redDir)

    # Change back to original directory
    os.chdir(redDir)

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
            trim=False, weight=None, fwhm_max=0, strehl_trim_min = None,
            submaps=0,
            fixDAR=True, use_koa_weather=False,
            mask=True,
            clean_dirs=None, combo_dir=None,
            instrument=instruments.default_inst,
            debug_plot_correlation=False
           ):
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
        If set to default=0 and trim=True, then we use FWHM < 1.25 * FWHM.min().
    strehl_trim_min : float or None, default = None
        The minimum strehl value allowed when trimming is turned on. When a value
        is specified, both FWHM and strehl trimming are performed.
        Recommended *ONLY* when AO performance is bad leading to a tight core, but
        very poor PSF otherwise.
        If None, then will not trim on strehl. 
    submaps : int, default=0
        Set to the number of submaps to be made (def=0).
    fixDAR : boolean, default = True
        Whether or not to calculate and apply DAR correction coefficients.
    use_koa_weather : boolean, default = False
        If calculating DAR correction, this keyword specifies if the atmosphere
        conditions should be downloaded from the KOA weather data. If False,
        atmosphere conditions are downloaded from the MKWC CFHT data.
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
    waveDir = util.getcwd()
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
            
            try: 
                shutil.copy(source_clean_dir + 'c' + source_file_root + '.rcoo',
                        dest_clean_dir + 'c' + dest_file_root + '.rcoo')
            except:
                print('No rcoo files to copy')
            
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
        if strehl_trim_min is not None:
            roots, strehls, fwhm, weights = trim_on_strehl(roots, strehls, fwhm, strehl_trim_min)

    ##########
    # Weighting
    ##########
    if weight == 'strehl':
        weights = weight_by_strehl(roots, strehls)

    if ((weight is not None) and (weight != 'strehl')):
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
    diffPA, phis = combine_rotation(cleanDir, roots, instrument=instrument)

    # Make a table of coordinates for the reference source.
    # These serve as initial estimates for the shifts.
    #combine_ref(_out + '.coo', cleanDir, roots, diffPA, refImage_index)
    combine_coo(_out + '.coo', cleanDir, roots, diffPA, refImage_index)

    # Keep record of files that went into this combine
    combine_lis(_out + '.lis', cleanDir, roots, diffPA)


    # Register images to get shifts.
    shiftsTab = combine_register(_out, refImage, diffPA, instrument=instrument,
                                 plot_correlation=debug_plot_correlation)

    # Determine the size of the output image from max shifts
    xysize = combine_size(shiftsTab, refImage, _out, _sub, submaps, phis)

    ##########
    # Sort frames -- recall that submaps assume sorted by FWHM.
    ##########
    roots, strehls, fwhm, weights, shiftsTab = sort_frames(roots, strehls, fwhm, weights, shiftsTab)
               
    # Combine all the images together.
    combine_drizzle(xysize, cleanDir, roots, _out, weights, shiftsTab,
                    wave, diffPA, fixDAR=fixDAR, mask=mask, instrument=instrument,
                    use_koa_weather=use_koa_weather)

    # Now make submaps
    if (submaps > 0):
        combine_submaps(xysize, cleanDir, roots, _sub, weights,
                        shiftsTab, submaps, wave, diffPA, fixDAR=fixDAR,
                        mask=mask, instrument=instrument,
                        use_koa_weather=use_koa_weather)

    # Remove *.lis_r file & rotated rcoo files, if any - these
    # were just needed to get the proper shifts for xregister
    #_lisr = _out + '.lis_r'
    #util.rmall([_lisr])
    #for i in range(len(allroots)):
    #    _rcoo = cleanDir + 'c' + str(allroots[i]) + '.rcoo'
    #    util.rmall([_rcoo])
    
    # Change back to original directory
    os.chdir(redDir)


def rot_img(root, phi, cleanDir, edit_header_PA = False):
    """
    Rotate image with scipy.ndimage.rotate. Only the image is rotated.
    The header WCS is not modified. The output imageis pre-pended with 'r'
    instead of 'c'.

    Parameters
    ----------
    root : str
        Root name of the file to rotate. Do not include the prefix (e.g. 'c').
    phi : float
        Angle to rotate the input image to get to PA=0.
    cleanDir : str
        The clean directory to find the input image and save the output rotated image.
    edit_header_PA : bool
        Default is False.
    """

    inCln = cleanDir + 'c' + root + '.fits'
    outCln = cleanDir + 'r' + root + '.fits'

    in_img, in_hdr = fits.getdata(inCln, header=True)
    in_wcs, rot_mat, new_coord_mat = rotate_wcs(in_hdr, phi)
    
    if in_wcs.wcs.has_cd():
        new_cd = new_coord_mat
    elif in_wcs.wcs.has_pc():
        new_pc = new_coord_mat

    # Rotate the image
    print('Rotating frame: ',root)
    #in_img[np.where(np.isnan(in_img) == True)] = 0
    out_img = rotate(in_img, -phi, order=3, mode='constant', cval=0, reshape=False)
    
    
    out_hdr = copy.deepcopy(in_hdr)
    
    
    # Don't use the default to_header for these
    # params since it switches PC to CD
    # https://github.com/astropy/astropy/issues/1084
    if in_wcs.wcs.has_cd():  # CD matrix
        out_hdr['CD1_1'] = new_cd[0][0]
        out_hdr['CD1_2'] = new_cd[0][1]
        out_hdr['CD2_1'] = new_cd[1][0]
        out_hdr['CD2_2'] = new_cd[1][1]
    
    elif in_wcs.wcs.has_pc():  # PC matrix + CDELT
        out_hdr['PC1_1'] = new_pc[0][0]
        out_hdr['PC1_2'] = new_pc[0][1]
        out_hdr['PC2_1'] = new_pc[1][0]
        out_hdr['PC2_2'] = new_pc[1][1]
    
    #CCD to image transform
    ltm = np.dot(rot_mat, [[1, 0],[0, 1]])
    out_hdr['LTM1_1'] = ltm[0][0]
    out_hdr['LTM1_2'] = ltm[1][0]
    out_hdr['LTM2_1'] = ltm[0][1]
    out_hdr['LTM2_2'] = ltm[1][1]

    # CCD to image tranform part 2
    # Assuming rotation around center
    img_size_x = in_hdr['NAXIS1']
    img_size_y = in_hdr['NAXIS2']
    orig_arr = np.concatenate((img_size_x/2*np.ones((2,1)), img_size_y/2*np.ones((2,1))), axis =1)
    add_arr = np.array([img_size_x/2, img_size_y/2])
    rotated_arr = np.dot(-rot_mat, orig_arr) + add_arr
    out_hdr['LTV1'] = rotated_arr[0][0]
    out_hdr['LTV2'] = rotated_arr[1][1]
    
    # Deletes the incorrectly generated PC
    # if relevant (https://github.com/astropy/astropy/issues/1084)
    if not in_wcs.wcs.has_pc():
        try:
            del out_hdr['PC1_1']
            del out_hdr['PC1_2']
            del out_hdr['PC2_1']
            del out_hdr['PC2_2']
        except:
            pass

    if edit_header_PA:
        out_hdr['PA'] = out_hdr['ROTPOSN'] + phi
        
    fits.writeto(outCln, out_img, out_hdr, output_verify=outputVerify,
                    overwrite=True)

    return

def rotate_wcs(hdr, phi):
    in_wcs = wcs.WCS(hdr)

    # Rotate the WCS
    theta = np.deg2rad(phi)
    sina = np.sin(theta)
    cosa = np.cos(theta)
    rot_mat = np.array([[cosa, -sina],
                        [sina, cosa]])
    
    if in_wcs.wcs.has_cd():  # CD matrix
        new_cd = np.dot(rot_mat, in_wcs.wcs.cd)
        in_wcs.wcs.cd = new_cd
        in_wcs.wcs.set()
        return in_wcs, rot_mat, new_cd

    elif in_wcs.wcs.has_pc():  # PC matrix + CDELT
        new_pc = np.dot(rot_mat, in_wcs.wcs.get_pc())
        in_wcs.wcs.pc = new_pc
        in_wcs.wcs.set()
        return in_wcs, rot_mat, new_pc
        
    else:
        raise TypeError("Unsupported wcs type (only CD or PC matrix allowed)")

    

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
    os.chdir(redDir)

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

def trim_on_strehl(roots, strehls, fwhm, strehl_trim_min):
    """
    Take a list of files and trim based on the strehl. 

    The returned arrays contain only those files that pass the above criteria.
    """
    # Pull out those we want to include in the combining
    keep = np.where(strehls >= strehl_trim_min)[0]
    strehls = strehls[keep]
    fwhm = fwhm[keep]
    roots = [roots[i] for i in keep]
    weights = np.array( [1.0/len(roots)] * len(roots) )

    print('combine: Keeping %d frames with strehl > %4.1f' \
        % (len(roots), strehl_trim_min))

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
                    wave, diffPA, fixDAR=True, use_koa_weather=False,
                    mask=True, instrument=instruments.default_inst,
                   ):
    
    _fits = outroot + '.fits'
    _tmpfits = outroot + '_tmp.fits'
    _wgt = outroot + '_sig.fits'
    _dlog = outroot + '_driz.log'
    _maxFile = outroot + '.max'

    util.rmall([_fits, _tmpfits, _wgt, _dlog])

    satLvl_combo = 0.0

    # Variable to store weighted sum of MJDs
    mjd_weightedSum = 0.0

    # Get the distortion maps for this instrument.
    hdr0 = fits.getheader(cleanDir + 'c' + roots[0] + '.fits')
    distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr0)
    
    print('combine: drizzling images together')
    f_dlog = open(_dlog, 'a')
    kernel = 'lanczos3'
    driz = drizzle.resample.Drizzle(kernel = kernel,
                    out_shape = (imgsize, imgsize), #np.shape(cdwt_img),
                    fillval = 0
                    )

    imgsize_half = imgsize / 2.0
    
    f_dlog.write('- Base dir: ' + cleanDir + '\n')                 
    for i in range(len(roots)):
        f_dlog.write(time.ctime() + '\n')
        f_dlog.write('- {} is image {} to be drizzled'.format(roots[i], i) + '\n')
        
        # Cleaned image
        _c = cleanDir + 'c' + roots[i] + '.fits'

        # Cleaned but distorted image
        _cd = cleanDir + 'distort/cd' + roots[i] + '.fits'
        _cdwt = cleanDir + 'weight/cdwt' + roots[i] + '.fits'

        util.rmall([_cdwt])

        # Multiply each distorted image by it's weight
        fits_cd = fits.open(_cd)
        cdwt_img = fits_cd[0].data * weights[i]   # Keep image variable for drizzling later
        fits_cd[0].data = cdwt_img                # Also save the weighted image to a file.
        fits_cd[0].header[instrument.hdr_keys['itime']] *= weights[i]
        fits_cd.writeto(_cdwt, output_verify=outputVerify)

        # Padd size (must be same as in combine_size)
        imgsize_orig = np.max(cdwt_img.shape)
        padd_1percent = int(np.floor(imgsize_orig * 0.01))
        padd_bigger_img = (imgsize - imgsize_orig) / 2.0

        # Get pixel shifts
        xsh = shifts[i][1]
        ysh = shifts[i][2]

        # For the first image, read in the header, otherwise use
        # the loaded in header from the previously drizzled image
        if i == 0:
            hdr = fits.getheader(_cd, ignore_missing_end=True)

        # Read in PA of each file to feed into drizzle for rotation
        hdr_current_img = fits.getheader(_cd, ignore_missing_end=True)
        phi = instrument.get_position_angle(hdr_current_img)
        if (diffPA == 1):
             _, _, cd_mat = rotate_wcs(hdr_current_img, phi)
             #cd_mat = rot_img(roots[i], phi, cleanDir, return_cd_only = True)

        # Setup distortion maps to include both geometric and differential atmospheric
        # refraction corrections (if specified).
        if (fixDAR == True):
            darRoot = _cdwt.replace('.fits', 'geo')
            (_xgeoim, _ygeoim) = dar.darPlusDistortion(
                                   _cdwt, darRoot,
                                   xgeoim=distXgeoim,
                                   ygeoim=distYgeoim,
                                   instrument=instrument,
                                   use_koa_weather=use_koa_weather)
        else:
            _xgeoim = distXgeoim
            _ygeoim = distYgeoim

        f_dlog.write('- Input data image: clean' + _cdwt.split('/clean')[1] + '\n')
        f_dlog.write('- X-shift distortion image: clean' + _xgeoim.split('/clean')[1] + '\n')
        f_dlog.write('- Y-shift distortion image: clean' + _ygeoim.split('/clean')[1] + '\n')

        # Get exposure time
        itime_keyword = 'ITIME'
        exp_time = hdr_current_img[itime_keyword]

        # Read in MJD of current file from FITS header
        mjd = instrument.get_mjd(hdr)
        mjd_weightedSum += weights[i] * mjd
        
        # weight the image by multiplying by mask
        # this is what is said to be done by in_mask in the iraf version 
        # (https://ftp.eso.org/scisoft/scisoft4/sources/iraf/extern/eis/doc/drizzle.hlp.html)
        if (mask == True):
            _mask = cleanDir + 'masks/mask' + roots[i] + '.fits'
            mask_img = fits.getdata(_mask)
            wgt_in = np.ones(np.shape(mask_img))*mask_img
            f_dlog.write('- Mask image: clean' + _mask.split('/clean')[1] + '\n')
        else:
            wgt_in = np.ones(np.shape(cdwt_img))
        
        print('Drizzling: ', roots[i])
        print('     xsh = {0:8.2f}'.format( xsh ))
        print('     ysh = {0:8.2f}'.format( ysh ))
        print('  weight = {0:8.2f}'.format( weights[i] ))
        print('   outnx = {0:8d}'.format( imgsize ))

        # We tell it the input is distorted/shifted and we want to undistort it
        wcs_in = wcs.WCS(hdr_current_img)
        wcs_out = wcs.WCS(hdr_current_img)

        if (diffPA == 1):
            f_dlog.write('- Rotating image. phi = {} \n'.format(phi))
            wcs_in.wcs.cd = cd_mat

        f_dlog.write('- Shifting image. xshift = {0:8.2f}, yshift = {1:8.2f} \n'.format(xsh, ysh))
        wcs_in.wcs.crpix = [cdwt_img.shape[0]/2.0, cdwt_img.shape[1]/2.0]
        wcs_out.wcs.crpix = [(imgsize / 2.0) + xsh, (imgsize / 2.0) + ysh]
    
        xgeoim = fits.getdata(_xgeoim).astype('float32')
        ygeoim = fits.getdata(_ygeoim).astype('float32')
    
        xdist = wcs.DistortionLookupTable( xgeoim, [0, 0], [0, 0], [1, 1])
        ydist = wcs.DistortionLookupTable( ygeoim, [0, 0], [0, 0], [1, 1])
    
        wcs_in.cpdis1 = xdist
        wcs_in.cpdis2 = ydist

        pixmap = drizzle.utils.calc_pixmap(wcs_in, wcs_out)
        
        # Drizzle this file ontop of all previous ones
        
        # Catches case when exposure time is a fraction of a second
        if exp_time > 0 and exp_time < 1:
            wht_scale = exp_time
            pixfrac = 1.0
            driz.add_image(cdwt_img, pixmap = pixmap, 
                                weight_map = wgt_in,
                                exptime = 1,
                                xmax = int(imgsize),
                                ymax = int(imgsize),
                                wht_scale = wht_scale,
                                pixfrac = pixfrac,
                                in_units = 'counts')
        else:
            wht_scale = 1.0
            pixfrac = 1.0
            driz.add_image(cdwt_img, pixmap = pixmap,
                                weight_map = wgt_in,
                                exptime = exp_time,
                                xmax = int(imgsize),
                                ymax = int(imgsize),
                                wht_scale = wht_scale,
                                pixfrac = pixfrac,
                                in_units = 'counts')
        f_dlog.write('- Drizzling onto full output image. Kernel: ' + kernel + '\n')
    
        #swtich from output cps to counts by multiplying by total counts
        out_img = driz.out_img * driz._texptime
        img_hdu = fits.PrimaryHDU(data=out_img, header=hdr)

        # import pylab as plt
        # plt.figure(1)
        # plt.clf()
        # plt.imshow(out_img, vmin=0, vmax=100)
        # plt.show()
        

        # make header
        if i == 0:
            # set CRPIX by the first image
            shift_crpix1 = wcs_out.wcs.crpix[0] + (imgsize - np.shape(cdwt_img)[0])/2 + xsh
            shift_crpix2 = wcs_out.wcs.crpix[1] + (imgsize - np.shape(cdwt_img)[0])/2 + ysh
            img_hdu.header.set('CRPIX1', shift_crpix1 + wcs_in.cpdis1.get_offset(shift_crpix1, shift_crpix2))
            img_hdu.header.set('CRPIX2', shift_crpix2 + wcs_in.cpdis2.get_offset(shift_crpix1, shift_crpix2))
        img_hdu.header.set('D{0:03d}VER'.format(i + 1), 'DRIZZLE VERSION {}'.format(drizzle.__version__))
        img_hdu.header.set('D{0:03d}DATA'.format(i + 1), 'clean' + _cdwt.split('/clean')[1], 'Drizzle, input data image')
        img_hdu.header.set('D{0:03d}DEXP'.format(i + 1), exp_time, 'Drizzle, input image exposure time (s)')
        img_hdu.header.set('D{0:03d}OUDA'.format(i + 1), 'combo' + _tmpfits.split('/combo')[1], 'Drizzle, output data image')
        img_hdu.header.set('D{0:03d}OUWE'.format(i + 1), 'combo' + _wgt.split('/combo')[1], 'Drizzle, output weighting image')
        img_hdu.header.set('D{0:03d}OUCO'.format(i + 1), '', 'Drizzle, output context image')
        if (mask == True):
            img_hdu.header.set('D{0:03d}MASK'.format(i + 1), 'clean' + _mask.split('/clean')[1], 'Drizzle, input mask')
        else:
            img_hdu.header.set('D{0:03d}MASK'.format(i + 1), '', 'Drizzle, input mask')
        img_hdu.header.set('D{0:03d}WTSC'.format(i + 1), wht_scale, 'Drizzle, weighting factor for input image')
        img_hdu.header.set('D{0:03d}KERN'.format(i + 1), kernel, 'Drizzle, form of weight distribution kernel')
        img_hdu.header.set('D{0:03d}PIXF'.format(i + 1), pixfrac, 'Drizzle, linear size of drop')
        img_hdu.header.set('D{0:03d}COEF'.format(i + 1), '', 'Drizzle, coefficients file name')
        img_hdu.header.set('D{0:03d}XGIM'.format(i + 1), 'clean' + _xgeoim.split('/clean')[1], 'Drizzle, X distortion image name')
        img_hdu.header.set('D{0:03d}YGIM'.format(i + 1), 'clean' + _ygeoim.split('/clean')[1], 'Drizzle, Y distortion image name')
        img_hdu.header.set('D{0:03d}SCAL'.format(i + 1), 1, 'Drizzle, scale (pixel size) of output image')
        if (diffPA == 1):
            img_hdu.header.set('D{0:03d}ROT'.format(i + 1), phi, 'Drizzle, rotation angle, degrees anticlockwise')
        else:
            img_hdu.header.set('D{0:03d}ROT'.format(i + 1), 0, 'Drizzle, rotation angle, degrees anticlockwise')
        img_hdu.header.set('D{0:03d}XSH'.format(i + 1), xsh, 'Drizzle, X shift applied')
        img_hdu.header.set('D{0:03d}YSH'.format(i + 1), ysh, 'Drizzle, Y shift applied')
        img_hdu.header.set('D{0:03d}SFTU'.format(i + 1), 'pixels', 'Drizzle, units used for shifts (output or input)')
        img_hdu.header.set('D{0:03d}SFTF'.format(i + 1), 'pixels', 'Drizzle, frame in which shifts were applied') #this might be wrong
        img_hdu.header.set('D{0:03d}EXKY'.format(i + 1), itime_keyword, 'Drizzle, exposure keyword name in input image')
        img_hdu.header.set('D{0:03d}INUN'.format(i + 1), 'counts', 'Drizzle, units of input image - counts or cps')
        img_hdu.header.set('D{0:03d}OUUN'.format(i + 1), 'counts', 'Drizzle, units of output image - counts or cps')
        img_hdu.header.set('D{0:03d}FVAL'.format(i + 1), '0', 'Drizzle, fill value for zero weight output pixel')
                       
        img_hdu.writeto(_tmpfits, output_verify='ignore', 
                                    overwrite=True)

        hdr = img_hdu.header

        wgt_hdu = fits.PrimaryHDU(data=driz.out_wht, header=hdr)
        wgt_hdu.writeto(_wgt, output_verify='ignore', 
                                    overwrite=True)

        f_dlog.write('- Output data image: combo' + _tmpfits.split('/combo')[1] + '\n')
        f_dlog.write('- Output weight image: combo' + _wgt.split('/combo')[1] + '\n')

        # Read .max file with saturation level for final combined image
        # by weighting each individual satLevel and summing.
        # Read in each satLevel from individual .max files
        _max = cleanDir + 'c' + roots[i] + '.max'
        getsatLvl = open(_max)
        satLvl = float(getsatLvl.read())
        getsatLvl.close()
        satLvl_wt = satLvl * weights[i]
        satLvl_combo += satLvl_wt

    f_dlog.write('Writing final images')
    print(_cdwt)
    print(_tmpfits)

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
                                          maxiters=5)
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
        if 'PA' in fits_f[0].header:
            fits_f[0].header.set('PA', phi,
                              'PA set by KAI')

    # Add keyword with distortion image information
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
    
    fits_f[0].header.set(
        'MJD-OBS', mjd_weightedMean,
        'Weighted modified julian date of combined observations'
    )
    fits_f[0].header.set(
        'MJD', mjd_weightedMean,
        'Weighted modified julian date of combined observations'
    )
    
    ## Also update date field in header
    fits_f[0].header.set(
        'DATE', '{0}'.format(time_obs.fits),
        'Weighted observation date'
    )

    # Put in number of drizzle images
    fits_f[0].header.set('NDRIZIM', len(roots), 
                         'Drizzle, number of images drizzled onto this out')

    # save weight file
    f_dlog.write('- Output weighting image: combo' + _wgt.split('/combo')[1] + '\n')
    wgt_hdu = fits.PrimaryHDU(data=driz.out_wht, header=hdr)
    wgt_hdu.writeto(_wgt, output_verify='ignore', 
                                overwrite=True)
    
    # Save to final fits file.
    f_dlog.write('- Output data image: combo' + _tmpfits.split('/combo')[1] + '\n')
    fits_f[0].writeto(_fits, output_verify=outputVerify)
    util.rmall([_tmpfits])

    for i in range(len(roots)):
        _cdwt = cleanDir + 'weight/cdwt' + roots[i] + '.fits'
        util.rmall([_cdwt])

    f_dlog.close()
                       

def combine_submaps(
        imgsize, cleanDir, roots, outroot, weights,
        shifts, submaps, wave, diffPA,
        fixDAR=True, use_koa_weather=False,
        mask=True, instrument=instruments.default_inst,
    ):
    """
    Assumes the list of roots are pre-sorted based on quality. Images are then
          divided up with every Nth image going into the Nth submap.

    mask: (def=True) Set to false for maser mosaics since they have only
          one image at each positions. Masking produces artifacts that
          Starfinder can't deal with.
    """

    extend = []
    for i in range(1 ,submaps+1):
        extend.append('_{}'.format(i))
    #extend = ['_1', '_2', '_3']
    _out = [outroot + end for end in extend]
    _fits = [o + '.fits' for o in _out]
    _tmp = [o + '_tmp.fits' for o in _out]
    _wgt = [o + '_sig.fits' for o in _out]
    _log = [o + '_driz.log' for o in _out]
    _max = [o + '.max' for o in _out]
    output_hdrs = [{} for o in _out]

    util.rmall(_fits + _tmp + _wgt + _log + _max)

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

    # Make one drizzle object per submap
    driz = []
    kernel = 'lanczos3'
    for i in range(submaps):
        driz.append(drizzle.resample.Drizzle(kernel = kernel,
                        out_shape = (imgsize, imgsize),
                        fillval = 0
                        ))

    for log in f_log:
        log.write('- Base dir: ' + cleanDir + '\n')
        
    for i in range(len(roots)):
        # Cleaned image
        _c = cleanDir + 'c' + roots[i] + '.fits'

        # Cleaned but distorted image
        _cd = cleanDir + 'distort/cd' + roots[i] + '.fits'
        cdwt = cleanDir + 'weight/cdwt' + roots[i] + '.fits'

        # Multiply each distorted image by it's weight
        util.rmall([cdwt])

        fits_cd = fits.open(_cd)
        cdwt_img = fits_cd[0].data * weights[i]   # Keep image variable for drizzling later
        fits_cd[0].data = cdwt_img
        fits_cd[0].header[instrument.hdr_keys['itime']] *= weights[i]
        fits_cd.writeto(cdwt, output_verify=outputVerify)
        
        # Padd size (must be same as in combine_size)
        imgsize_orig = np.max(cdwt_img.shape)
        padd_1percent = int(np.floor(imgsize_orig * 0.01))
        padd_bigger_img = (imgsize - imgsize_orig) / 2.0
        
        # Get pixel shifts
        xsh = shifts[i][1]
        ysh = shifts[i][2]
        
        # Determine which submap we should be drizzling to.
        sub = int(i % submaps)
        fits_im = _tmp[sub]
        wgt = _wgt[sub]
        log = f_log[sub]

        log.write(time.ctime() + '\n')
        log.write('- {} is image {} to be drizzled'.format(roots[i], i) + '\n')
        
        # For the first image of each submap, read in the header, otherwise use
        # the loaded in header from the previously drizzled image
        img_in_submap = int(i/submaps)
        if img_in_submap == 0:
            hdr = fits.getheader(_c, ignore_missing_end=True)
        else:
            hdr = fits.getheader(fits_im, ignore_missing_end=True)
            
        # Read in PA of each file to feed into drizzle for rotation
        hdr_current_img = fits.getheader(_c, ignore_missing_end=True)
        # Each submap will build its header on the first image in submap
        if bool(output_hdrs[sub]) == False:
            output_hdrs[sub] = hdr_current_img
        phi = instrument.get_position_angle(hdr_current_img)
        if (diffPA == 1):
             _, _, cd_mat = rotate_wcs(hdr_current_img, phi)

        # Calculate saturation level for submaps
        # by weighting each individual satLevel and summing.
        # Read in each satLevel from individual .max files
        max_indiv = cleanDir + 'c' + roots[i] + '.max'
        satfile = open(max_indiv)
        satLvl = float(satfile.read())
        satLvl_wt = satLvl * weights[i]
        satLvl_tot[sub] += satLvl_wt
        
        # Add up the weights that go into each submap
        weightsTot[sub] += weights[i]
        
        satLvl_sub[sub] = satLvl_tot[sub] / weightsTot[sub]
        
        if (fixDAR == True):
            darRoot = cdwt.replace('.fits', 'geo')
            print('submap: ',cdwt)
            (_xgeoim, _ygeoim) = dar.darPlusDistortion(
                                   cdwt, darRoot,
                                   xgeoim=distXgeoim,
                                   ygeoim=distYgeoim,
                                   instrument=instrument,
                                   use_koa_weather=use_koa_weather)
        else:
            _xgeoim = distXgeoim
            _ygeoim = distYgeoim

        log.write('- Input data image: clean' + cdwt.split('/clean')[1] + '\n')
        log.write('- X-shift distortion image: clean' + _xgeoim.split('/clean')[1] + '\n')
        log.write('- Y-shift distortion image: clean' + _ygeoim.split('/clean')[1] + '\n')

        # Get exposure time
        itime_keyword = 'ITIME'
        exp_time = hdr_current_img[itime_keyword]

        # Read in MJD of current file from FITS header
        mjd = instrument.get_mjd(hdr)
        mjd_weightedSums[sub] += weights[i] * mjd
        
        # Drizzle this file ontop of all previous ones.

        # weight the image by multiplying by mask
        # this is what is said to be done by in_mask in the iraf version 
        # (https://ftp.eso.org/scisoft/scisoft4/sources/iraf/extern/eis/doc/drizzle.hlp.html)
        if (mask == True):
            _mask = cleanDir + 'masks/mask' + roots[i] + '.fits'
            mask_img = fits.getdata(_mask)
            wgt_in = np.ones(np.shape(mask_img))*mask_img
            log.write('- Mask image: clean' + _mask.split('/clean')[1] + '\n')
        else:
            wgt_in = np.ones(np.shape(cdwt_img))
            
        # We tell it the input its distorted/shfited and we want to undistort it
        wcs_in = wcs.WCS(hdr_current_img)
        wcs_out = wcs.WCS(hdr_current_img)

        log.write('- Shifting image. xshift = {0:8.2f}, yshift = {1:8.2f} \n'.format(xsh, ysh))
        wcs_in.wcs.crpix = [cdwt_img.shape[0]/2.0, cdwt_img.shape[1]/2.0]
        wcs_out.wcs.crpix = [(imgsize / 2.0) + xsh, (imgsize / 2.0) + ysh]
        
        xgeoim = fits.getdata(_xgeoim).astype('float32')
        ygeoim = fits.getdata(_ygeoim).astype('float32')
    
        xdist = wcs.DistortionLookupTable( xgeoim, [0, 0], [0, 0], [1, 1])
        ydist = wcs.DistortionLookupTable( ygeoim, [0, 0], [0, 0], [1, 1])
    
        wcs_in.cpdis1 = xdist
        wcs_in.cpdis2 = ydist
    
        pixmap = drizzle.utils.calc_pixmap(wcs_in, wcs_out)

        # Catches case when exposure time is a fraction of a second
        if exp_time > 0 and exp_time < 1:
            wht_scale = exp_time
            pixfrac = 1.0
            driz[sub].add_image(cdwt_img, pixmap = pixmap, 
                                weight_map = wgt_in,
                                exptime = 1,
                                xmax = int(imgsize),
                                ymax = int(imgsize),
                                wht_scale = wht_scale,
                                pixfrac = pixfrac,
                                in_units = 'counts')
        else:
            wht_scale = 1.0
            pixfrac = 1.0
            driz[sub].add_image(cdwt_img, pixmap = pixmap, 
                                weight_map = wgt_in,
                                exptime = exp_time,
                                xmax = int(imgsize),
                                ymax = int(imgsize),
                                wht_scale = wht_scale,
                                pixfrac = pixfrac,
                                in_units = 'counts')
        log.write('- Drizzling onto full output image. Kernel: ' + kernel + '\n')
        
        #swtich from output cps to counts by multiplying by total counts
        out_img = driz[sub].out_img * driz[sub]._texptime
        img_hdu = fits.PrimaryHDU(data=out_img, header=hdr)

        # make header
        if img_in_submap == 0:
            # set CRPIX by the first image
            shift_crpix1 = wcs_out.wcs.crpix[0] + (imgsize - np.shape(cdwt_img)[0])/2 + xsh
            shift_crpix2 = wcs_out.wcs.crpix[1] + (imgsize - np.shape(cdwt_img)[0])/2 + ysh
            img_hdu.header.set('CRPIX1', shift_crpix1 + wcs_in.cpdis1.get_offset(shift_crpix1, shift_crpix2))
            img_hdu.header.set('CRPIX2', shift_crpix2 + wcs_in.cpdis2.get_offset(shift_crpix1, shift_crpix2))
        img_hdu.header.set('D{0:03d}VER'.format(img_in_submap + 1), 'DRIZZLE VERSION {}'.format(drizzle.__version__))
        img_hdu.header.set('D{0:03d}DATA'.format(img_in_submap + 1), 'clean' + cdwt.split('/clean')[1], 'Drizzle, input data image')
        img_hdu.header.set('D{0:03d}DEXP'.format(img_in_submap + 1), exp_time, 'Drizzle, input image exposure time (s)')
        img_hdu.header.set('D{0:03d}OUDA'.format(img_in_submap + 1), 'combo' + fits_im.split('/combo')[1], 'Drizzle, output data image')
        img_hdu.header.set('D{0:03d}OUWE'.format(img_in_submap + 1), 'combo' + wgt.split('/combo')[1], 'Drizzle, output weighting image')
        img_hdu.header.set('D{0:03d}OUCO'.format(img_in_submap + 1), '', 'Drizzle, output context image')
        if (mask == True):
            img_hdu.header.set('D{0:03d}MASK'.format(img_in_submap + 1), 'clean' + _mask.split('/clean')[1], 'Drizzle, input mask')
        else:
            img_hdu.header.set('D{0:03d}MASK'.format(img_in_submap + 1), '', 'Drizzle, input mask')
        img_hdu.header.set('D{0:03d}WTSC'.format(img_in_submap + 1), wht_scale, 'Drizzle, weighting factor for input image')
        img_hdu.header.set('D{0:03d}KERN'.format(img_in_submap + 1), kernel, 'Drizzle, form of weight distribution kernel')
        img_hdu.header.set('D{0:03d}PIXF'.format(img_in_submap + 1), pixfrac, 'Drizzle, linear size of drop')
        img_hdu.header.set('D{0:03d}COEF'.format(img_in_submap + 1), '', 'Drizzle, coefficients file name')
        img_hdu.header.set('D{0:03d}XGIM'.format(img_in_submap + 1), 'clean' + _xgeoim.split('/clean')[1], 'Drizzle, X distortion image name')
        img_hdu.header.set('D{0:03d}YGIM'.format(img_in_submap + 1), 'clean' + _ygeoim.split('/clean')[1], 'Drizzle, Y distortion image name')
        img_hdu.header.set('D{0:03d}SCAL'.format(img_in_submap + 1), 1, 'Drizzle, scale (pixel size) of output image')
        if (diffPA == 1):
            img_hdu.header.set('D{0:03d}ROT'.format(img_in_submap + 1), phi, 'Drizzle, rotation angle, degrees anticlockwise')
        else:
            img_hdu.header.set('D{0:03d}ROT'.format(img_in_submap + 1), 0, 'Drizzle, rotation angle, degrees anticlockwise')
        img_hdu.header.set('D{0:03d}XSH'.format(img_in_submap + 1), xsh, 'Drizzle, X shift applied')
        img_hdu.header.set('D{0:03d}YSH'.format(img_in_submap + 1), ysh, 'Drizzle, Y shift applied')
        img_hdu.header.set('D{0:03d}SFTU'.format(img_in_submap + 1), 'pixels', 'Drizzle, units used for shifts (output or input)')
        img_hdu.header.set('D{0:03d}SFTF'.format(img_in_submap + 1), 'pixels', 'Drizzle, frame in which shifts were applied') #this might be wrong
        img_hdu.header.set('D{0:03d}EXKY'.format(img_in_submap + 1), itime_keyword, 'Drizzle, exposure keyword name in input image')
        img_hdu.header.set('D{0:03d}INUN'.format(img_in_submap + 1), 'counts', 'Drizzle, units of input image - counts or cps')
        img_hdu.header.set('D{0:03d}OUUN'.format(img_in_submap + 1), 'counts', 'Drizzle, units of output image - counts or cps')
        img_hdu.header.set('D{0:03d}FVAL'.format(img_in_submap + 1), '0', 'Drizzle, fill value for zero weight output pixel')
        
        img_hdu.writeto(fits_im, output_verify='ignore', 
                                    overwrite=True)
        
        log.write('- Output data image: combo' + fits_im.split('/combo')[1] + '\n')
    
    # Calculate weighted MJDs for each submap
    mjd_weightedMeans = mjd_weightedSums / weightsTot
    submaps_time_obs = Time(mjd_weightedMeans, format='mjd')
        
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
                                          maxiters=5)
        sci_mean = tmp_stats[0]
        sci_stddev = tmp_stats[2]

        # Find and fix really bad pixels
        idx = np.where(fits_f[0].data < (sci_mean - 10*sci_stddev))
        fits_f[0].data[idx] = 0.0
        
        f_log[s].write('- Removed {} bad pixels'.format(len(idx)) + '\n')

        # Normalize properly
        fits_f[0].data = fits_f[0].data / weightsTot[s]

        # Fix the ITIME header keyword so that it matches (weighted).
        itime = fits_f[0].header.get('ITIME')
        itime /= weightsTot[s]
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

        # Fix the DATASEC header keyword, if it exists.
        if 'DATASEC' in fits_f[0].header:
            fits_f[0].header['DATASEC'] = '[1:{0:d},1:{0:d}]'.format(imgsize)
        
        # Store weighted MJDs in header
        fits_f[0].header.set(
            'MJD-OBS',
            mjd_weightedMeans[s],
            'Weighted modified julian date of combined observations',
        )
        fits_f[0].header.set(
            'MJD',
            mjd_weightedMeans[s],
            'Weighted modified julian date of combined observations',
        )
        
        ## Also update date field in header
        fits_f[0].header.set(
            'DATE',
            '{0}'.format(submaps_time_obs[s].fits),
            'Weighted observation date',
        )

        # Deletes ref pixels and strehl pixel ref on submaps
        # CRITICAL for making starfinder run
        del fits_f[0].header['XREF']
        del fits_f[0].header['YREF']
        del fits_f[0].header['XSTREHL']
        del fits_f[0].header['YSTREHL']

        f_log[s].write('Writing final images \n')
        
        # Write out final submap fits file
        f_log[s].write('- Output data image: combo' + _fits[s].split('/combo')[1] + '\n')
        fits_f[0].writeto(_fits[s], output_verify=outputVerify)

        f_log[s].write('- Output weighting image: combo' + _wgt[s].split('/combo')[1] + '\n')
        wgt_hdu = fits.PrimaryHDU(data=driz[s].out_wht, header=hdr)
        wgt_hdu.writeto(_wgt[s], output_verify='ignore', 
                                    overwrite=True)
    
    for f in f_log:
        f.close()
        
    util.rmall(_tmp)
    for i in range(len(roots)):
        _cdwt = cleanDir + 'weight/cdwt' + roots[i] + '.fits'
        util.rmall([_cdwt])
    #util.rmall([cdwt])


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
    phis = np.zeros(len(clean_files))

    for cc in range(len(clean_files)):
        hdr = fits.getheader(clean_files[cc], ignore_missing_end=True)
        phi = instrument.get_position_angle(hdr)
        phis[cc] = phi

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
            phis[cc] = phi
            rot_img(roots[cc], phi, cleanDir)

    return (diffPA, phis)

def sort_frames(roots, strehls, fwhm, weights, shiftsTab):
    sidx = np.argsort(fwhm)

    # Make sorted lists.
    strehls = strehls[sidx]
    fwhm = fwhm[sidx]
    weights = weights[sidx]
    roots = [roots[i] for i in sidx]
    shiftsX = shiftsTab['xshift']
    shiftsX = shiftsX[sidx]
    shiftsY = shiftsTab['yshift']
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

def xregister_correlation_fourier(I, R, Nx, Ny):
    """
    Pythonified version of algorithm from iraf xregister
    (see algorithms https://astro.uni-bonn.de/~sysstw/lfa_html/iraf/images.xregister.html)
    """
    sum_I = np.sum(I) / (Nx * Ny)
    sum_R = np.sum(R) / (Nx * Ny)
    
    sumsqI = np.sqrt(np.sum((I - sum_I) ** 2))
    sumsqR = np.sqrt(np.sum((R - sum_R) ** 2))
    
    FFTI = np.fft.fft2((I - sum_I) / sumsqI)
    FFTR = np.fft.fft2((R - sum_R) / sumsqR)
    
    correlation = np.fft.ifft2(FFTR * np.conj(FFTI))

    return correlation

def elliptical_gaussian_2d(params, x, y):
    amplitude, x0, y0, sigma_x, sigma_y, theta, offset = params
    x_prime = (x - x0) * np.cos(theta) + (y - y0) * np.sin(theta)
    y_prime = -(x - x0) * np.sin(theta) + (y - y0) * np.cos(theta)
    return amplitude * np.exp(-((x_prime / sigma_x) ** 2 + (y_prime / sigma_y) ** 2) / 2) + offset

def residuals(params, x, y, image):
    # Residual function for least_squares
    return elliptical_gaussian_2d(params, x, y) - image

def combine_register(outroot, refImage, diffPA, plot_correlation = False, instrument = instruments.default_inst):
    """
    Register all the images through a cross-correlation method.
    """

    shiftFile = outroot + '.shifts'

    print('combine: registering images')
    if (diffPA == 1):
        input = '@' + outroot + '.lis_r'
    else:
        input = '@' + outroot + '.lis'

    # Fetch the reference image, which all shifts will be calculated against.
    # Also get the coordinates of the cooStar in the reference image (coo_file_ref)
    ref_img = fits.getdata(refImage)
    
    # coo file starts with c even if image was rotated
    refImage_filename = refImage.split('/')[-1]
    if refImage_filename[0] == 'r':
        coo_file_ref = Table.read(refImage.split(refImage_filename)[0] + 'c' + refImage_filename[1:-5] + '.rcoo',
                                  format='ascii', header_start=None)
    else:
        coo_file_ref = Table.read(refImage[:-5] + '.coo', format='ascii', header_start=None)

    # Fetch the file names and cooStar position for all images.
    fileNames = Table.read(input[1:], format='ascii.no_header')
    fileNames = np.array(fileNames)
    fileNames = np.array(fileNames, dtype='S')
    
    coords = Table.read(outroot + '.coo', format='ascii', header_start=None)

    # Setup a shifts table for storing the final calculated shifts.
    shiftsTable_empty = np.zeros((len(fileNames), 3), dtype=float)
    shiftsTable = Table(shiftsTable_empty, dtype=('S50', float, float))
    
    hdrRef = fits.getheader(refImage, ignore_missing_end=True)

    # Calculate the size of the large-crop box (crop-val) in pixels.
    plate_scale = instrument.get_plate_scale(hdrRef) #arcsec/pixels
    crop_val = 1 / plate_scale # 1 arcsec/(arcsec/pix), units = pixels
    crop_val = int(np.round(crop_val))

    # For speed, we pre-trim the reference image edges by 5%.
    # This helps remove edge effects. 
    five_percent = int(np.shape(ref_img)[0]*0.05)
    ref_img_noedge = ref_img[five_percent:-five_percent, five_percent:-five_percent]
    
    # Loop through the images and calculate the shift for each (w.r.t. the reference image).
    for ii in range(len(fileNames)):
        fileName = fileNames[ii].decode("utf-8")
        
        # Load up the image, and cooStar coordinates.
        shift_img = fits.getdata(fileName)
        fileName_filename = fileName.split('/')[-1]
        if fileName_filename[0] == 'r':
            coo_name = fileName.split(fileName_filename)[0] + 'c' + fileName_filename[1:-5] + '.rcoo'
        else:
            coo_name = fileName[:-5] + '.coo'
        coo_file = Table.read(coo_name, format='ascii', header_start=None)

        # Crop off 5% of the edges of ii-th image
        # to avoid edge effects contaminating the cross correlation.
        shift_img_noedge = shift_img[five_percent:-five_percent, five_percent:-five_percent]
        
        # Figure out the overlap window between this ii-th image and the ref image.
        xmin_ref = 0
        xmax_ref = ref_img_noedge.shape[1]
        ymin_ref = 0
        ymax_ref = ref_img_noedge.shape[0]

        xmin_ii = 0
        xmax_ii = shift_img_noedge.shape[1]
        ymin_ii = 0
        ymax_ii = shift_img_noedge.shape[0]
        
        # Apply a shift to the image to center up on the cooStar.
        xshift = int(coo_file_ref['col1'][0] - coo_file['col1'][0])
        yshift = int(coo_file_ref['col2'][0] - coo_file['col2'][0])
        #global_shift_img = shift(shift_img_noedge, (xshift, yshift), mode = 'constant')  # shifted image.

        xmin_ii += xshift
        xmax_ii += xshift
        ymin_ii += yshift
        ymax_ii += yshift

        # Calculate the vertices of the overlap rectangle. Note these are in the coordinate system
        # of the 5% cropped reference image (ref_img_noedge). 
        overlap_xmin = max(xmin_ref, xmin_ii)
        overlap_xmax = min(xmax_ref, xmax_ii)
        overlap_ymin = max(ymin_ref, ymin_ii)
        overlap_ymax = min(ymax_ref, ymax_ii)

        # For sanity and ease of use later on, make sure the overlap rectangle is divisible by 2.
        if (overlap_xmax - overlap_xmin) % 2 != 0:
            overlap_xmax -= 1
        if (overlap_ymax - overlap_ymin) % 2 != 0:
            overlap_ymax -= 1

        # Double-check that we have overlap in both directions.
        if overlap_xmin >= overlap_xmax or overlap_ymin >= overlap_ymax:
            raise RuntimeError(f'Two images do not overlap: {fileName_filename} and {refImage_filename}')

        # Rectangle vertices in the original coordinates of shift_img_noedge:
        overlap_xmin_ii = overlap_xmin - xshift
        overlap_xmax_ii = overlap_xmax - xshift
        overlap_ymin_ii = overlap_ymin - yshift
        overlap_ymax_ii = overlap_ymax - yshift

        # Return cropped images of both the reference and ii-th image that contains only
        # the overlapping area. This is what should be used as input into the cross-correlation.
        # Note that this action should also center up the two images onto the same coordinate system.
        ref_img_cropped = ref_img_noedge[overlap_ymin:overlap_ymax, overlap_xmin:overlap_xmax]
        ii_img_cropped = shift_img_noedge[overlap_ymin_ii:overlap_ymax_ii, overlap_xmin_ii:overlap_xmax_ii]

        from skimage.registration import phase_cross_correlation
        if plot_correlation:
            plot_corr_filename = fileName.replace('.fits', '_corr.png')
        else:
            plot_corr_filename = None

        shift_pcc = subpixel_windowed_pcc(ref_img_cropped, ii_img_cropped, window_radius=100,
                                          upsample_factor=50, alpha=0.5,
                                          plot_corr_map=plot_correlation, plot_outfile=plot_corr_filename)
        total_y_shift_pcc = shift_pcc[0] + yshift
        total_x_shift_pcc = shift_pcc[1] + xshift
        pcc_half_size = np.array(ref_img_cropped.shape) / 2.0

        print('*** Comparing correlation methods:')
        print(f'\t {"":<10s} {"x":<10s} {"y":<10s}')
        print(f'\t {"coo":<10s} {xshift:10.2f} {yshift:10.2f}')
        # print(f'\t {"manual":<10s} {total_x_shift:10.2f} {total_y_shift:8.3f}')
        print(f'\t {"pcc":<10s} {total_x_shift_pcc:10.2f} {total_y_shift_pcc:10.2f}')

        if plot_correlation:
            import pylab as plt
            import matplotlib.colors as mcolors
            
            norm = mcolors.SymLogNorm(vmin=-10, vmax=1e4, linthresh=100) 
            
            # plot the image cutouts used in cross-correlation.
            plt.close(1)
            plt.figure(1, figsize=(10, 5.2))
            plt.clf()
            plt.subplots_adjust(left=0.1, right=0.95, bottom=0.10, top=0.85)
            ax1 = plt.subplot(1, 2, 1)
            plt.imshow(ref_img_cropped, cmap='gist_heat', norm=norm)
            plt.axis('equal')
            plt.xlim(0, ref_img_cropped.shape[1])
            plt.xlim(0, ref_img_cropped.shape[0])
            plt.title(f'Cropped Ref:\n{refImage_filename}')
            
            plt.subplot(1, 2, 2)
            plt.imshow(ii_img_cropped, cmap='gist_heat', norm=norm)
            plt.axis('equal')
            plt.xlim(0, ref_img_cropped.shape[1])
            plt.xlim(0, ref_img_cropped.shape[0])
            plt.title(f'Cropped ii Img:\n{fileName_filename}')
            save_name1 = fileName.replace('.fits', '_crop.png')
            plt.savefig(save_name1)

        # Save the final shifts to the shift table.
        shiftsTable[ii][0] = fileName.split('/')[-1]
        shiftsTable[ii][1] = total_x_shift_pcc
        shiftsTable[ii][2] = total_y_shift_pcc


    # Save out the shifts table. 
    util.rmall([shiftFile])
    shiftsTable.rename_column('col0', 'file')
    shiftsTable.rename_column('col1', 'xshift')
    shiftsTable.rename_column('col2', 'yshift')
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

def combine_size(shiftsTable, refImage, outroot, subroot, submaps, phis):
    """Determine the final size of the fully combined image, accounting for
    both shifts and rotations. Use the shifts stored in the shiftsTable and
    rotation angles in phis.
    
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
    @type submaps: int
    @param phis: Array of rotation angles in degrees for each image
    @type phis: array-like
    """
    x_allShifts = np.array(shiftsTable['xshift'])
    y_allShifts = np.array(shiftsTable['yshift'])
    phis_rad = np.deg2rad(np.array(phis))
    
    # Get original image size
    orig_img = fits.getdata(refImage)
    orig_size = (orig_img.shape)[0]
    half_size = orig_size / 2.0
    

    # Calculate the four corners of the original image relative to its center
    corners = np.array([
        [-half_size, -half_size],
        [half_size, -half_size],
        [half_size, half_size],
        [-half_size, half_size]
    ])  # Shape: (4, 2)
    
    # Compute cos and sin for all angles
    cos_phi = np.cos(phis_rad)  # Shape: (n_images,)
    sin_phi = np.sin(phis_rad)  # Shape: (n_images,)
    
    # Rotate corners for all images at once
    corners_x = corners[:, 0][:, np.newaxis]  # Shape: (4, 1)
    corners_y = corners[:, 1][:, np.newaxis]  # Shape: (4, 1)

    # Unrotate shifts (inverse rotation)
    unrotated_xsh = x_allShifts * cos_phi + y_allShifts * sin_phi
    unrotated_ysh = -x_allShifts * sin_phi + y_allShifts * cos_phi
    
    # Add unrotated shifts to corners
    shifted_x = corners_x + unrotated_xsh
    shifted_y = corners_y + unrotated_ysh
    
    rotated_x = shifted_x * cos_phi - shifted_y * sin_phi  # Shape: (4, n_images)
    rotated_y = shifted_x * sin_phi + shifted_y * cos_phi  # Shape: (4, n_images)
    
    # Find the range needed: max positive and max negative extent
    xmax = np.max(rotated_x)
    xmin = np.min(rotated_x)
    ymax = np.max(rotated_y)
    ymin = np.min(rotated_y)

    # Maximum offset needed from center in any direction
    maxoffset = max(abs(xmax), abs(xmin), abs(ymax), abs(ymin))
    
    # Add padding
    padd = int(np.floor(orig_size * 0.01))

    
    # Read in the reference star's position in the ref image and translate
    # it into the coordinates of the final main and sub maps.
    hdr = fits.getheader(refImage, ignore_missing_end=True)
    xrefSrc = float(hdr['XREF'])
    yrefSrc = float(hdr['YREF'])
    xrefSrc = xrefSrc + ((maxoffset + padd) - half_size)
    yrefSrc = yrefSrc + ((maxoffset + padd) - half_size)
    
    cooMain = [outroot + '.coo']
    cooSubs = ['%s_%d.coo' % (subroot, i) for i in range(submaps+1)]
    cooAll = cooMain + cooSubs
    util.rmall(cooAll)
    
    for coo in cooAll:
        _allCoo = open(coo, 'w')
        _allCoo.write('%9.3f %9.3f\n' % (xrefSrc, yrefSrc))
        _allCoo.close()

    
    xysize = int((maxoffset + padd) * 2.0) #int(float(orig_size) + ((maxoffset + padd) * 2.0))
    
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
    #from pyraf import iraf as ir

    # Setup the drizzle parameters we will use
    drizzle.outweig = ''
    drizzle.in_mask = ''
    drizzle.wt_scl = 1
    drizzle.outnx = imgsize
    drizzle.outny = imgsize
    drizzle.pixfrac = 1
    drizzle.kernel = 'lanczos3'
    drizzle.scale = 1
    drizzle.shft_un = 'input'
    drizzle.shft_fr = 'output'
    drizzle.align = 'center'
    drizzle.expkey = 'ITIME'
    drizzle.in_un = 'counts'
    drizzle.out_un = 'counts'

def clean_drizzle(xgeoim, ygeoim, _bp, _cd, _wgt, _dlog,
        fixDAR=True, instrument=instruments.default_inst,
        use_koa_weather=False):

    # Get the distortion maps for this instrument.

    bp_file = fits.open(_bp)
    hdr = bp_file[0].header
    bp_img = bp_file[0].data
    distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr)
    itime_keyword = 'ITIME'
    exp_time = hdr[itime_keyword]

    # Input image size
    imgsizeX = float(hdr['NAXIS1'])
    imgsizeY = float(hdr['NAXIS2'])
    if (imgsizeX >= imgsizeY):
        imgsize = imgsizeX
    else:
        imgsize = imgsizeY
    outnx = imgsize
    outny = imgsize

    if (fixDAR == True):
        darRoot = _cd.replace('.fits', 'geo')

        (_xgeoim, _ygeoim) = dar.darPlusDistortion(
                               _bp, darRoot, xgeoim, ygeoim,
                               instrument=instrument,
                               use_koa_weather=use_koa_weather)
    else:
        _xgeoim = distXgeoim
        _ygeoim = distYgeoim

    wgt_in = np.ones((int(outnx),int(outny)))
    wcs_in = wcs.WCS(hdr)
    wcs_out = wcs.WCS(hdr)

    xgeoim = fits.getdata(_xgeoim).astype('float32')
    ygeoim = fits.getdata(_ygeoim).astype('float32')

    xdist = wcs.DistortionLookupTable( xgeoim, [0, 0], [0, 0], [1, 1])
    ydist = wcs.DistortionLookupTable( ygeoim, [0, 0], [0, 0], [1, 1])

    wcs_in.cpdis1 = xdist
    wcs_in.cpdis2 = ydist

    pixmap = drizzle.utils.calc_pixmap(wcs_in, wcs_out)

    kernel = 'lanczos3'
    driz = drizzle.resample.Drizzle(kernel = kernel,
                                    out_shape = np.shape(bp_img),
                                    fillval = 0
                                    )

    # Catches case when exposure time is a fraction of a second
    if exp_time > 0 and exp_time < 1:
        wht_scale = exp_time
        pixfrac = 1.0
        driz.add_image(bp_img, pixmap = pixmap, 
                            exptime = 1,
                            xmax = int(outnx),
                            ymax = int(outny),
                            wht_scale = wht_scale,
                            pixfrac = pixfrac,
                            in_units = 'counts')
    else:
        wht_scale = 1.0
        pixfrac = 1.0
        driz.add_image(bp_img, pixmap = pixmap, 
                            exptime = exp_time,
                            xmax = int(outnx),
                            ymax = int(outny),
                            wht_scale = wht_scale,
                            pixfrac = pixfrac,
                            in_units = 'counts')

    #swtich from output cps to counts by multiplying by total counts
    out_img = driz.out_img * driz._texptime
    img_hdu = fits.PrimaryHDU(data=out_img, header=hdr)

    # make header
    img_hdu.header.set('CRPIX1', wcs_in.wcs.crpix[0] + wcs_in.cpdis1.get_offset(wcs_in.wcs.crpix[0], wcs_in.wcs.crpix[1]))
    img_hdu.header.set('CRPIX2', wcs_in.wcs.crpix[1] + wcs_in.cpdis2.get_offset(wcs_in.wcs.crpix[0], wcs_in.wcs.crpix[1]))
    img_hdu.header.set('NDRIZIM', 1, 'Drizzle, number of images drizzled onto this out')
    img_hdu.header.set('D001VER', 'DRIZZLE VERSION {}'.format(drizzle.__version__))
    img_hdu.header.set('D001DATA', _bp, 'Drizzle, input data image')
    img_hdu.header.set('D001DEXP', driz._texptime, 'Drizzle, input image exposure time (s)')
    img_hdu.header.set('D001OUDA', _cd, 'Drizzle, output data image')
    img_hdu.header.set('D001OUWE', _wgt, 'Drizzle, output weighting image')
    img_hdu.header.set('D001OUCO', '', 'Drizzle, output context image')
    img_hdu.header.set('D001MASK', '', 'Drizzle, input weighting image')
    img_hdu.header.set('D001WTSC', wht_scale, 'Drizzle, weighting factor for input image')
    img_hdu.header.set('D001KERN', kernel, 'Drizzle, form of weight distribution kernel')
    img_hdu.header.set('D001PIXF', pixfrac, 'Drizzle, linear size of drop')
    img_hdu.header.set('D001COEF', '', 'Drizzle, coefficients file name')
    img_hdu.header.set('D001XGIM', _xgeoim, 'Drizzle, X distortion image name')
    img_hdu.header.set('D001YGIM', _ygeoim, 'Drizzle, Y distortion image name')
    img_hdu.header.set('D001SCAL', 1, 'Drizzle, scale (pixel size) of output image')
    img_hdu.header.set('D001ROT', 0, 'Drizzle, rotation angle, degrees anticlockwise')
    img_hdu.header.set('D001XSH', 0, 'Drizzle, X shift applied')
    img_hdu.header.set('D001YSH', 0, 'Drizzle, Y shift applied')
    img_hdu.header.set('D001SFTU', 'pixels', 'Drizzle, units used for shifts (output or input)')
    img_hdu.header.set('D001SFTF', 'pixels', 'Drizzle, frame in which shifts were applied') #this might be wrong
    img_hdu.header.set('D001EXKY', itime_keyword, 'Drizzle, exposure keyword name in input image')
    img_hdu.header.set('D001INUN', 'counts', 'Drizzle, units of input image - counts or cps')
    img_hdu.header.set('D001OUUN', 'counts', 'Drizzle, units of output image - counts or cps')
    img_hdu.header.set('D001FVAL', '0', 'Drizzle, fill value for zero weight output pixel')

    img_hdu.writeto(_cd, output_verify='ignore', 
                                overwrite=True)
    #hdulist=fits.open(_cd)
    #img_data = hdulist['SCI'].data
    #img_header = hdulist['SCI'].header

    #wgt_data = hdulist['WHT'].data
    #wgt_header = hdulist['WHT'].header
    wgt_hdu = fits.PrimaryHDU(data=driz.out_wht, header=hdr)
    wgt_hdu.writeto(_wgt, output_verify='ignore', 
                                overwrite=True)


def cosmicray_median(ccd, input_mask, error_image=None, thresh=5, mbox=5, gbox=0,
                     rbox=0, fratio = 5, star_thresh=4, thresh_in_star = 10):
    """
    Modified from ccdproc
    
    Identify cosmic rays through median technique. The median technique
    identifies cosmic rays by identifying pixels by subtracting a median image
    from the initial data array.

    Parameters
    ----------
    ccd : `numpy.ndarray` or `numpy.ma.MaskedArray`
        Data to have cosmic ray cleaned.

    thresh : float, optional
        Threshold for detecting cosmic rays.
        Default is ``5``.

    error_image : `numpy.ndarray`, float or None, optional
        Error level. If None, the task will use the standard
        deviation of the data. If an ndarray, it should have the same shape
        as data.
        Default is ``None``.

    mbox : int, optional
        Median box for detecting cosmic rays.
        Default is ``11``.

    gbox : int, optional
        Box size to grow cosmic rays. If zero, no growing will be done.
        Default is ``0``.

    rbox : int, optional
        Median box for calculating replacement values. If zero, no pixels will
        be replaced.
        Default is ``0``.

    fratio : int, optional
        Flux ratio determination for cosmic rays in stars by comparing residual
        to nearby median.
        Default is ``5``.

    star_thresh : int, optional
        Threshold for determining if pixels belong to a star as compared to 
        their 25th percentile background.
        Default is ``4``.

    thresh_in_star : int, optional
        Threshold for cosmicrays in star.
        Default is ``10".

    Notes
    -----
    Similar implementation to crmedian in iraf.imred.crutil.crmedian.

    Returns
    -------
    nccd : `~astropy.nddata.CCDData` or `numpy.ndarray`
        An object of the same type as ccd is returned. If it is a
        `~astropy.nddata.CCDData`, the mask attribute will also be updated with
        areas identified with cosmic rays masked.

    nccd : `numpy.ndarray`
        If an `numpy.ndarray` is provided as ccd, a boolean ndarray with the
        cosmic rays identified will also be returned.

    Examples
    --------
    1) Given an numpy.ndarray object, the syntax for running
       cosmicray_median would be:

       >>> newdata, mask = cosmicray_median(data, error_image=error,
       ...                                  thresh=5, mbox=11,
       ...                                  rbox=11, gbox=5)   # doctest: +SKIP

       where error is an array that is the same shape as data but
       includes the pixel error. This would return a data array, newdata,
       with the bad pixels replaced by the local median from a box of 11
       pixels; and it would return a mask indicating the bad pixels.

    2) Given an `~astropy.nddata.CCDData` object with an uncertainty frame, the syntax
       for running cosmicray_median would be:

       >>> newccd = cosmicray_median(ccd, thresh=5, mbox=11,
       ...                           rbox=11, gbox=5)   # doctest: +SKIP

       The newccd object will have bad pixels in its data array replace and the
       mask of the object will be created if it did not previously exist or be
       updated with the detected cosmic rays.
    """
    data = ccd

    if error_image is None:
        error_image = data.std()
    else:
        if not isinstance(error_image, (float, np.ndarray, np.float32)):
            raise TypeError('error_image is not a float or ndarray.')

    # create the median image
    marr = ndimage.median_filter(data, size=(mbox, mbox))


    # Only look at the data array
    if isinstance(data, np.ma.MaskedArray):
        data = data.data

    # Find the residual image
    # Compare residual to the error in the image
    rarr = (data - marr) /error_image

    # Compare residual to the mean nearby
    # in order to exclude stars
    rarr2 = (data - marr) /marr

    local_bg = ndimage.percentile_filter(data, 25, size=mbox*3)  # Larger box for background
    star_indicator = (marr - local_bg) / error_image
                         
    # identify all sources
    # Different criteria for stars (star indicator > star thresh) and for background (star indicator < star thresh)
    # Larger criteria to find a cosmic ray on a star
    crarr = ((rarr > thresh_in_star) & (rarr2 > fratio) & (star_indicator > star_thresh)) | ((rarr > thresh) & (star_indicator <= star_thresh))

    # Remove reference pixels
    ref_pixels = np.where(input_mask == 2)
    crarr[ref_pixels] = 0

    # grow the pixels
    if gbox > 0:
        crarr = ndimage.maximum_filter(crarr, gbox)

    # replace bad pixels in the image
    if rbox > 0:
        # Assert shape of data is a square
        assert np.shape(data)[0] == np.shape(data)[1]
        ndata = loop_through_crs(data, crarr, dim = np.shape(data)[0])

    return ndata, crarr

def loop_through_crs(data, crarr, search_box = 5, interp_box = 10, dim = 1024):
    # Box to search for cosmic rays must be smaller than the area you're interpretting over
    assert search_box < interp_box
    
    masked_data = np.ma.masked_array(data, (crarr == 1))
    new_data = copy.deepcopy(data)

    crarr_xs = np.where(crarr == True)[0]
    crarr_ys = np.where(crarr == True)[1]

    mark_done = np.ones((dim, dim))
    mark_done[crarr == 1] = 0

    filled_vals = masked_data.filled(np.nan)
    
    for crarr_x, crarr_y in zip(crarr_xs, crarr_ys):
        if mark_done[crarr_x, crarr_y] == 1:
            continue
        
        # Set edges of search box and round to edges of image 
        min_search_x = crarr_x - search_box
        max_search_x = crarr_x + search_box
        min_search_y = crarr_y - search_box
        max_search_y = crarr_y + search_box

        min_search_x, max_search_x = round_to_edge(min_search_x, max_search_x, 0, dim)
        min_search_y, max_search_y = round_to_edge(min_search_y, max_search_y, 0, dim)

        # Set values in search box (which will be filled in) as done
        mark_done[min_search_x:max_search_x, min_search_y:max_search_y] = 1

        # Set edges of interp box and round to edges of image 
        min_interp_x = crarr_x - interp_box
        max_interp_x = crarr_x + interp_box
        min_interp_y = crarr_y - interp_box
        max_interp_y = crarr_y + interp_box

        min_interp_x, max_interp_x = round_to_edge(min_interp_x, max_interp_x, 0, dim)
        min_interp_y, max_interp_y = round_to_edge(min_interp_y, max_interp_y, 0, dim)

        # Cut out cosmic rays
        patch_x = np.arange(min_interp_x, max_interp_x)
        patch_y = np.arange(min_interp_y, max_interp_y)
        masked_data_patch = masked_data[min_interp_x:max_interp_x, min_interp_y:max_interp_y]
        x1, y1 = np.meshgrid(patch_x, patch_y, indexing='ij')
        x = x1[~masked_data_patch.mask]
        y = y1[~masked_data_patch.mask]
        filtered_data = masked_data_patch[~masked_data_patch.mask]

        # Prepare patch related to where we're looking for cosmic rays
        search_patch_x = np.arange(min_search_x, max_search_x)
        search_patch_y = np.arange(min_search_y, max_search_y)
        x1_search_patch, y1_search_patch = np.meshgrid(search_patch_x, search_patch_y, indexing='ij')
        masked_data_search_patch = masked_data[min_search_x:max_search_x, min_search_y:max_search_y]

        # Interpolate over masked data and find interpolated values
        # over search area (where we've looked for cosmic rays)
        linear_interp = griddata((x, y), filtered_data.ravel(), (x1_search_patch,y1_search_patch), fill_value = 0)#np.nan)

        # Set masked out cosmic ray values to interpolated values
        new_data_patch = new_data[min_search_x:max_search_x, min_search_y:max_search_y]
        new_data_patch[masked_data_search_patch.mask] = linear_interp[masked_data_search_patch.mask]
        new_data[min_search_x:max_search_x, min_search_y:max_search_y] = new_data_patch

    return new_data


def round_to_edge(min_box_val, max_box_val, min_val, max_val):
    if min_box_val < min_val:
        min_box_val = min_val
    if max_box_val > max_val:
        max_box_val = max_val

    return min_box_val, max_box_val
    
def clean_cosmicrays(_ff, _mask, wave, _input_mask, thresh=5, mbox=5, rbox=10, fratio = 5, gbox = 0, star_thresh=4, thresh_in_star = 10):
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
    ff_header = fits.getheader(_ff)

    input_mask = fits.getdata(_input_mask)
    
    tmp_stats = stats.sigma_clipped_stats(ff_img,
                                          sigma_upper=2, sigma_lower=5,
                                          maxiters=5)
    mean = tmp_stats[0]
    stddev = tmp_stats[2]
    
    newdata, crmask = cosmicray_median(ff_img, input_mask, error_image = stddev, thresh=thresh, mbox=mbox, gbox=gbox, rbox=rbox, fratio=fratio,
                                      star_thresh = star_thresh, thresh_in_star = thresh_in_star)
    crmask = crmask.astype(int)

    ff_header.set('CRCOR', 'removed={}, thresh={}, mbox={}, gbox={}, rbox={}, fratio={}, star_thresh={}, thresh_in_star={}'.format(np.sum(crmask), thresh, mbox, gbox, rbox, fratio, star_thresh, thresh_in_star))

    # Save to a temporary file.
    fits.writeto(_mask, crmask, output_verify=outputVerify)
    fits.writeto(_ff, newdata, header=ff_header, output_verify=outputVerify, overwrite = True)

    return stddev

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
                                          maxiters=5)
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
        update_fits=True,cent_box=12,
        offset_method='aotsxy',pcuxyRef=None):
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

    check_loc : bool, default=True
        If True the reference source is recentered for this frame.
        Use False if the offsets are large enough to move the reference source
        off of the image.
    update_fits : bool, default=True
        Update the fits files with the reference pixel values
    cent_box : float, default: 12
        Box size to center the source
    offset_method : str, default='aotsxy'
        Method to calculate offsets from reference image.
        Options are 'aotsxy' or 'radec'.
        In images where 'aotsxy' keywords aren't reliable, 'radec' calculated
        offsets may work better.
    """
    hdr = fits.getheader(_ce, ignore_missing_end=True)
    radec = instrument.get_radec(hdr)
    aotsxy = kai_util.getAotsxy(hdr)
    
    if offset_method == 'pcu':
        pcuxy = [float(hdr['PCSFX']), float(hdr['PCSFY'])]  #New version may be PCUX and PCUY
        pcu_scale = instrument.get_pcu_scale(hdr)

    # Determine the image's PA and plate scale
    phi = instrument.get_position_angle(hdr)
    scale = instrument.get_plate_scale(hdr)

    # Determine the instrument angle w.r.t. the AO bench.
    inst_angle = instrument.get_instrument_angle(hdr)

    # Calculate the pixel offsets from the reference image
    if offset_method == 'radec':
        d_xy = kai_util.radec2pix(radec, phi, scale, radecRef)
    elif offset_method == 'aotsxy':
        d_xy = kai_util.aotsxy2pix(aotsxy, scale, aotsxyRef,
                                   inst_angle=inst_angle)
    elif offset_method == 'pcu':
            d_xy = kai_util.pcuxy2pix(pcuxy, phi, pcu_scale, pcuxyRef)
    else:
        d_xy = kai_util.aotsxy2pix(aotsxy, scale, aotsxyRef,
                                   inst_angle=inst_angle)

    # In the new image, find the REF and STRL coords
    xref = refSrc[0] + d_xy[0]
    yref = refSrc[1] + d_xy[1]
    xstr = strSrc[0] + d_xy[0]
    ystr = strSrc[1] + d_xy[1]
    print('clean_makecoo: xref, yref start = {0:.2f} {1:.2f}'.format(xref, yref))

    # re-center stars to get exact coordinates
    if check_loc:
        image_data = fits.getdata(_ce)
        for _ in range(5):
            x0 = int(np.round(xref - (cent_box - 1)/2))
            y0 = int(np.round(yref - (cent_box - 1)/2))
            cutout = image_data[
                y0: y0 + cent_box,
                x0: x0 + cent_box
            ]
            dy, dx = ndimage.center_of_mass(cutout)
            xref = x0 + dx
            yref = y0 + dy

        for _ in range(5):
            x0 = int(np.round(xstr - (cent_box - 1)/2))
            y0 = int(np.round(ystr - (cent_box - 1)/2))
            cutout = image_data[
                y0: y0 + cent_box,
                x0: x0 + cent_box
            ]
            dy, dx = ndimage.center_of_mass(cutout)
            xstr = x0 + dx
            ystr = y0 + dy

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

    #file(_cc.replace('.fits', '.coo'), 'w').write('%7.2f  %7.2f\n' % (xref, yref))
    open(_cc.replace('.fits', '.coo'), 'w').write('%7.2f  %7.2f\n' % (xref, yref))

    # Make a temporary rotated coo file, in case there are any data sets
    # with various PAs; needed for xregister; remove later
    xyRef_rot = kai_util.rotate_coo(xref, yref, phi, instrument)
    xref_r = xyRef_rot[0]
    yref_r = xyRef_rot[1]

    xyStr_rot = kai_util.rotate_coo(xstr, ystr, phi, instrument)
    xstr_r = xyStr_rot[0]
    ystr_r = xyStr_rot[1]

    #file(_cc.replace('.fits', '.rcoo'), 'w').write('%7.2f  %7.2f\n' % (xref_r, yref_r))
    open(_cc.replace('.fits', '.rcoo'), 'w').write('%7.2f  %7.2f\n' % (xref_r, yref_r))

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
            n_img = fits.getdata(_n, ignore_missing_end=True)
            sci_stats = stats.sigma_clipped_stats(n_img,
                                                  sigma_upper=1, sigma_lower=10,
                                                  maxiters=20)
            sci_mean = sci_stats[0]

            sky_img = fits.getdata(_sky, ignore_missing_end=True)
            sky_stats = stats.sigma_clipped_stats(sky_img,
                                                  sigma_upper=5, sigma_lower=15,
                                                  maxiters=5)
            sky_mean = sky_stats[0]
            
            
            fact = sci_mean/sky_mean
            #print 'scaleSky: factor = %5f  sci_mean = %5f  sky_mean = %5f' % \
            #      (fact, sci_mean, sky_mean)
            fits_sky = fits.open(_sky)
            fits_sky[0].data *= fact
            fits_sky.writeto(self.skyName, output_verify=outputVerify)
        else:
            img = fits.getdata(_sky)
            fits.writeto(self.skyName, img, output_verify=outputVerify)

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
        sciAng = fits.getheader(_n)['ROTPPOSN']

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

        sky = self.skyDir + self.skies[skyidx] + ".fits"

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
                                              maxiters=4)
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
           trim=False, weight=None, fwhm_max=0, submaps=0,
           clean_dirs=None, combo_dir=None,
           fixDAR=True, maskSubmap=False,
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
        If set to default=0 and trim=True, then we use FWHM < 1.25 * FWHM.min().
    strehl_trim_min : float or None, default = None
        The minimum strehl value allowed when trimming is turned on. When a value
        is specified, both FWHM and strehl trimming are performed.
        Recommended *ONLY* when AO performance is bad leading to a tight core, but
        very poor PSF otherwise.
        If None, then will not trim on strehl. 
    submaps : int, default=0
        Set to the number of submaps to be made (def=0).
    fixDAR : boolean, default = True
        Whether or not to calculate and apply DAR correction coefficients.
    use_koa_weather : boolean, default = False
        If calculating DAR correction, this keyword specifies if the atmosphere
        conditions should be downloaded from the KOA weather data. If False,
        atmosphere conditions are downloaded from the MKWC CFHT data.
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
    waveDir = util.getcwd()
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
        outroot += outSuffix

    # Set up combo directory. This is the final output directory.
    comboDir = rootDir + 'combo/'
    
    if combo_dir is not None:
        comboDir = util.trimdir(os.path.abspath(combo_dir) + '/')
    
    util.mkdir(comboDir)

    
    # Make strings out of all the filename roots.
    allroots = instrument.make_filenames(files, prefix='')
    allroots = [aa.replace('.fits', '') for aa in allroots]

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
    diffPA, phis = combine_rotation(cleanDir, roots, instrument=instrument)

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
    BROKEN Register images for a mosaic. This only calculates the exact
    shifts between each image... it doesn't do the combining.

    @param outroot: The root for the output image. The resulting
    shifts will be written into a file called <outroot>.shifts
    @type outroot: string
    @param refImage: The name of the reference image.
    @type refImage: string
    """
    from pyraf import iraf as ir

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
    x_allShifts = shiftsTable['xshift']
    y_allShifts = shiftsTable['yshift']

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
    

def subpixel_windowed_pcc(ref_img, mov_img, window_radius=100, upsample_factor=10, alpha=0.5,
                          plot_corr_map=False, plot_outfile=None):
    """
    Finds sub-pixel shifts restricted to a window_radius box.
    """
    shape = np.array(ref_img.shape)
    
    #### 1. Compute Cross-Correlation (Coarse)
    src_freq = fftn(ref_img)
    target_freq = fftn(mov_img)
    cross_corr = ifftn(src_freq * target_freq.conj())
    
    #### 2. Apply Tukey Window to restrict the search area
    cross_corr_centered = fftshift(cross_corr)
    win_list = [windows.tukey(dim, alpha=alpha) for dim in shape]
    
    # (Simplified windowing: zeroing out anything beyond the window_radius)
    mask = np.zeros(shape)
    slices = tuple(slice(max(0, d//2 - window_radius), min(d, d//2 + window_radius)) for d in shape)
    mask[slices] = 1.0
    
    # Create the tapered mask
    full_window = np.ones(shape)
    for i, w in enumerate(win_list):
        slice_obj = [None] * len(shape)
        slice_obj[i] = slice(None)
        full_window *= w[tuple(slice_obj)]
    
    # Apply both the hard boundary and the taper
    constrained_corr = ifftshift(cross_corr_centered * full_window * mask)

    #### 3. Locate Coarse Peak
    max_idx = np.unravel_index(np.argmax(np.abs(constrained_corr)), shape)
    midpoints = shape // 2
    # Convert to shift (handling wrap-around)
    coarse_shift = np.array([i if i <= d // 2 else i - d for i, d in zip(max_idx, shape)])

    #### 4. Sub-pixel Refinement (Upsampled DFT)
    if upsample_factor > 1:
        # We "zoom in" around the coarse_shift using a Matrix Multiply DFT
        # This is the same math skimage uses for their sub-pixel refinement
        # Define the sub-pixel grid around the coarse peak
        upsampled_region_size = np.ceil(upsample_factor * 1.5).astype(int)
        dft_shift = coarse_shift * upsample_factor
        
        # Run upsampled DFT at the specific coarse location
        # (Simplified for brevity: we call the peak-finding logic on the small patch)
        # In practice, you'd use skimage.registration._upsampled_dft
        refined_shift = _upsampled_refinement(src_freq, target_freq, upsample_factor, coarse_shift)
    else:
        refined_shift = coarse_shift

    if plot_corr_map and plot_outfile is not None:
        # plot the correlation map.
        plt.figure(2)
        plt.clf()
        plt.imshow(np.abs(fftshift(constrained_corr)), cmap="viridis")
        # Plot center points (the 0.5 is just for plotting purposes.)
        plt.scatter(shape[1]/2 + 0.5, shape[0]/2 + 0.5,
                    color="green", marker="*", s=40, label="Coo Shift")
        plt.scatter(shape[1]/2 + refined_shift[1] + 0.5, shape[0]/2 + refined_shift[0] + 0.5,
                    color="red", marker="*", s=40, label="PCC Shift")
        plt.legend()
        plt.xlim(shape[1]/2 - window_radius/3, shape[1]/2 + window_radius/3)
        plt.ylim(shape[0]/2 - window_radius/3, shape[0]/2 + window_radius/3)

        plt.savefig(plot_outfile)
        
        
    return refined_shift.astype(float)


def _upsampled_refinement(src_freq, target_freq, upsample_factor, initial_shift):
    """
    Internal helper to perform the Matrix Multiply DFT centered on initial_shift.
    """
    # This mimics the core logic of skimage's subpixel refinement
    # It essentially 'zooms' into the correlation map at the specific pixel
    from skimage.registration._phase_cross_correlation import _upsampled_dft
    
    shape = src_freq.shape
    upsample_factor = float(upsample_factor)
    
    # Center the upsampled DFT on the coarse shift
    sample_region_offset = initial_shift * upsample_factor
    
    # Compute the cross-power spectrum in the zoomed region
    image_product = src_freq * target_freq.conj()
    
    # Upsample by factor around the initial peak
    upsampled_corr = _upsampled_dft(image_product, 
                                    upsampled_region_size=1.5, 
                                    upsample_factor=upsample_factor, 
                                    axis_offsets=initial_shift)
    
    # Find peak in the zoomed-in array
    max_idx = np.unravel_index(np.argmax(np.abs(upsampled_corr)), upsampled_corr.shape)
    
    # Adjust shift by the fractional amount
    fractional_shift = (np.array(max_idx) - 0.75 * upsample_factor) / upsample_factor
    
    return initial_shift + fractional_shift
