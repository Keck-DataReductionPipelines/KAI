import os, sys, shutil
from kai.reduce import util, lin_correction
from astropy.io import fits
from astropy import stats
from kai import instruments
import numpy as np
from astropy import stats
import astropy
from datetime import datetime
from pkg_resources import parse_version
import warnings
import pdb

module_dir = os.path.dirname(__file__)


def makedark(files, output,
             raw_dir=None,
             reduce_dir=None,
             instrument=instruments.default_inst):
    """
    Make dark image for imaging data. Makes a calib/ directory
    and stores all output there. All output and temporary files
    will be created in a darks/ subdirectory.

    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    output : str
        Output file name. Include the .fits extension.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    reduce_dir : str, optional
        Directory such as <epoch>/reduce/ with contents including
        the calib/, calib/darks/, etc. directories live.
        Files will be output into reduce_dir + calib/darks/.
        If epoch_dir is None, then use the current working directory.
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    rawDir, redDir = get_raw_reduce_directories(raw_dir, reduce_dir)
    curDir = redDir + 'calib/'
    darkDir = util.trimdir(curDir + 'darks/')

    util.mkdir(curDir)
    util.mkdir(darkDir)

    _out = darkDir + output
    _outlis = _out.replace('.fits', '.lis')
    util.rmall([_out, _outlis])

    darks_orig = instrument.make_filenames(files, rootDir=rawDir)

    # Write out the sources of the dark files
    data_sources_file = open(redDir + 'data_sources.txt', 'a')
    data_sources_file.write('---\n# Dark Files for {0} \n'.format(output))

    for cur_file in darks_orig:
        out_line = '{0} ({1})\n'.format(cur_file, datetime.now())
        data_sources_file.write(out_line)

    data_sources_file.close()

    # Create a file with the list of images that go into the dark
    f_on = open(_outlis, 'w')
    f_on.write('\n'.join(darks_orig) + '\n')
    f_on.close()

    img_stack = []
    hdr_stack = []

    print(f'makedark: Making dark {_out} with the following images:')
    for ii in range(len(darks_orig)):
        print(f'\t{darks_orig[ii]}')
        img, hdr = fits.getdata(darks_orig[ii], header=True)
        img_stack.append(img)
        hdr_stack.append(hdr)
    img_stack = np.array(img_stack)

    # Stack the darks with sigma clipping, median combine.
    dk_avg, dk_med, dk_std = astropy.stats.sigma_clipped_stats(img_stack,
                                                               cenfunc='median',
                                                               sigma_lower=3,
                                                               sigma_upper=3,
                                                               axis=0)

    # Save to an output file.
    fits.writeto(_out, dk_med, header=hdr_stack[0])

    return


def makeflat(onFiles, offFiles, output,
             dark_frame=None, normalizeFirst=False,
             raw_dir=None, reduce_dir=None,
             instrument=instruments.default_inst):
    """
    Make flat field image for imaging data. Makes a calib/ directory
    and stores all output there. All output and temporary files
    will be created in a flats/ subdirectory.

    If only twilight flats were taken (as in 05jullgs), use these flats as
    the onFiles, and use 0,0 for offFiles. So the reduce.py file should look
    something like this: onFiles = range(22, 26+1) and offFiles = range(0,0)
    The flat will then be made by doing a median combine using just the
    twilight flats.

    Parameters
    ----------
    onFiles : list of int
        Integer list of lamps ON files. Does not require padded zeros.
    offFiles : list of int
        Integer list of lamps OFF files. Does not require padded zeros.
    output : str
        Output file name. Include the .fits extension.
    dark_frame : str, default=None
        File name for the dark frame in order to carry out dark correction.
        If not provided, dark frame is not subtracted and a warning is thrown.
        Assumes dark file is located under ./calib/darks/
    normalizeFirst : bool, default=False
        If the individual flats should be normalized first,
        such as in the case of twilight flats.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    reduce_dir : str, optional
        Directory such as <epoch>/reduce/ with contents including
        the calib/, calib/darks/, etc. directories live.
        Files will be output into reduce_dir + calib/darks/.
        If epoch_dir is None, then use the current working directory.
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    rawDir, redDir = get_raw_reduce_directories(raw_dir, reduce_dir)
    curDir = redDir + 'calib/'
    flatDir = util.trimdir(curDir + 'flats/')

    util.mkdir(curDir)
    util.mkdir(flatDir)

    _on = flatDir + 'lampsOn.fits'
    _off = flatDir + 'lampsOff.fits'
    _out = flatDir + output
    _onlis = flatDir + 'on.lis'
    _offlis = flatDir + 'off.lis'

    util.rmall([_on, _off, _out, _onlis, _offlis])

    lampson = instrument.make_filenames(onFiles, rootDir=rawDir)
    lampsoff = instrument.make_filenames(offFiles, rootDir=rawDir)

    lampson_copied = instrument.make_filenames(onFiles, rootDir=flatDir)
    lampsoff_copied = instrument.make_filenames(offFiles, rootDir=flatDir)

    # Copy files - we will modify these in place.
    for ii in range(len(lampson_copied)):
        if os.path.exists(lampson_copied[ii]): os.remove(lampson_copied[ii])
        shutil.copy(lampson[ii], lampson_copied[ii])

    for ii in range(len(lampsoff_copied)):
        if os.path.exists(lampsoff_copied[ii]): os.remove(lampsoff_copied[ii])
        shutil.copy(lampsoff[ii], lampsoff_copied[ii])

    # Write out the sources of the flat files
    data_sources_file = open(redDir + 'data_sources.txt', 'a')

    data_sources_file.write(
        '---\n# Flat Files for {0}, Lamps On\n'.format(output))
    for cur_file in lampson:
        out_line = '{0} ({1})\n'.format(cur_file, datetime.now())
        data_sources_file.write(out_line)

    data_sources_file.write(
        '---\n# Flat Files for {0}, Lamps Off\n'.format(output))
    for cur_file in lampsoff:
        out_line = '{0} ({1})\n'.format(cur_file, datetime.now())
        data_sources_file.write(out_line)

    data_sources_file.close()

    # If dark frame is provided, load it up.
    if dark_frame is not None:
        dark_file = redDir + '/calib/darks/' + dark_frame

        # Read in dark frame data
        dark_data = fits.getdata(dark_file, ignore_missing_end=True)
    else:
        warning_message = 'Dark frame not provided for makeflat().'
        warning_message += '\nUsing flat frames without dark subtraction.'
        warnings.warn(warning_message)

    # Go through each flat file and subtract dark, linearity correct (inline, overwrite copied file).
    print(f'makeflat: Making flat lamps on with the following images:')
    for i in range(len(lampson_copied)):
        print(f'\t{lampson_copied[i]}')

        # Dark subtraction
        if dark_frame is not None:
            cur_flat = fits.open(lampson_copied[i], mode='update', ignore_missing_end=True, output_verify='ignore')
            cur_flat[0].data -= dark_data
            cur_flat.flush(output_verify='ignore')
            cur_flat.close(output_verify='ignore')

        # Linearity correction
        lin_correction.lin_correction(lampson_copied[i], instrument=instrument)

    print(f'makeflat: Making flat lamps on with the following images:')
    for i in range(len(lampsoff_copied)):
        print(f'\t{lampson_copied[i]}')

        # Dark subtraction
        if dark_frame is not None:
            cur_flat = fits.open(lampsoff_copied[i], mode='update', ignore_missing_end=True, output_verify='ignore')
            cur_flat[0].data -= dark_data
            cur_flat.flush(output_verify='ignore')
            cur_flat.close(output_verify='ignore')

        # Linearity correction
        lin_correction.lin_correction(lampsoff_copied[i], instrument=instrument)


    # Gather up and stack the lamps-off files. Final median-combined image is 'off_med'
    if len(lampsoff_copied) > 0:
        stack_off = []
        for ii in range(len(lampsoff_copied)):
            img = fits.getdata(lampsoff_copied[ii])
            stack_off.append(img)
        stack_off = np.array(stack_off)

        # Stack the lamps off flats with sigma clipping, median combine.
        off_avg, off_med, off_std = astropy.stats.sigma_clipped_stats(stack_off,
                                                               cenfunc='median',
                                                               sigma_lower=3,
                                                               sigma_upper=3,
                                                               axis=0)

        del stack_off, off_avg, off_std
    else:
        off_med = None

    # Load up the lamps on images
    stack_on = []
    hdr = None

    for ii in range(len(lampson_copied)):
        # Load images (and first header)
        if ii == 0:
            img, hdr = fits.getdata(lampson_copied[ii], header=True)
        else:
            img = fits.getdata(lampson_copied[ii])

        # Subtract the off-lamps (if we have them)
        if off_med is not None:
            img -= off_med

        # Normalize before combining if specified.
        if normalizeFirst:
            # Determine central area for normalization, avoid edges. Leave a 10% gap on each side.
            norm_gap_size = (0.1 * np.array(img.shape)).astype('int')
            norm_lo = norm_gap_size
            norm_hi = img.shape - norm_gap_size
            norm_value = np.median( img[norm_lo[0]:norm_hi[0], norm_lo[1]:norm_hi[1]] )

            img /= norm_value

        stack_on.append(img)

    stack_on = np.array(stack_on)

    # Stack the lamps off flats with sigma clipping, median combine.
    on_avg, on_med, on_std = astropy.stats.sigma_clipped_stats(stack_on,
                                                               cenfunc='median',
                                                               sigma_lower=3,
                                                               sigma_upper=3,
                                                               axis=0)

    # Normalize the final flat.
    norm_gap_size = (0.1 * np.array(on_med.shape)).astype('int')
    norm_lo = norm_gap_size
    norm_hi = on_med.shape - norm_gap_size
    norm_value = np.median(on_med[norm_lo[0]:norm_hi[0], norm_lo[1]:norm_hi[1]])

    on_med /= norm_value

    # Save to an output file.
    fits.writeto(_out, on_med, header=hdr)

    # Save some records of which files went into the stack.
    f_on = open(_onlis, 'w')
    f_on.write('\n'.join(lampson_copied) + '\n')
    f_on.close()

    if off_med is not None:
        f_off = open(_offlis, 'w')
        f_off.write('\n'.join(lampsoff_copied) + '\n')
        f_off.close()

    return


def makemask(dark, flat, output, reduce_dir=None, instrument=instruments.default_inst):
    """
    Make bad pixel mask for imaging data. Makes a calib/ directory
    and stores all output there. All output and temporary files
    will be created in a masks/ subdirectory.
    
    Parameters
    ----------
    dark : str
        The filename of the dark file (must be in the
        calib/darks/ directory). This is used to
        construct a hot pixel mask. Use a long (t>20sec) exposure dark.
    flat : str
        The filename of a flat file  (must be in the
        calib/flats/ directory). This is used to
        construct a dead pixel mask. The flat should be normalized.
    output : str
        The output file name. This will be created in the masks/
        subdirectory.
    reduce_dir : str, optional
        Directory such as <epoch>/reduce/ with contents including
        the calib/, calib/darks/, etc. directories live.
        Files will be output into reduce_dir + calib/darks/.
        If epoch_dir is None, then use the current working directory.
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    rawDir, redDir = get_raw_reduce_directories(None, reduce_dir)
    calDir = redDir + 'calib/'
    maskDir = util.trimdir(calDir + 'masks/')
    flatDir = util.trimdir(calDir + 'flats/')
    darkDir = util.trimdir(calDir + 'darks/')

    util.mkdir(calDir)
    util.mkdir(maskDir)

    _out = maskDir + output
    _dark = darkDir + dark
    _flat = flatDir + flat
    _inst_mask = module_dir + '/masks/' +  instrument.get_bad_pixel_mask_name()

    util.rmall([_out])

    ##########
    # Make hot pixel mask
    ##########
    whatDir = redDir + dark
    print(whatDir)
    print(f'makemask: Making mask {_out} with the following:')
    print(f'\t{_dark}')
    print(f'\t{_flat}')

    # Get the sigma-clipped mean and stddev on the dark
    img_dk = fits.getdata(_dark)
    if parse_version(astropy.__version__) < parse_version('3.0'):
        dark_stats = stats.sigma_clipped_stats(img_dk,
                                               sigma=3,
                                               iters=10)
    else:
        dark_stats = stats.sigma_clipped_stats(img_dk,
                                               sigma=3,
                                               maxiters=10)
        
    dark_mean = dark_stats[0]
    dark_stddev = dark_stats[2]

    # Clip out the very hot pixels.
    hi = dark_mean + (10.0 * dark_stddev)
    hot = img_dk > hi

    ##########
    # Make dead pixel mask
    ##########
    img_fl = fits.getdata(_flat)
    if parse_version(astropy.__version__) < parse_version('3.0'):
        flat_stats = stats.sigma_clipped_stats(img_dk,
                                               sigma=3,
                                               iters=10)
    else:
        flat_stats = stats.sigma_clipped_stats(img_fl,
                                           sigma=3,
                                           maxiters=10)
    flat_mean = flat_stats[0]
    flat_stddev = flat_stats[2]

    # Clip out the dead pixels
    lo = 0.5
    hi = flat_mean + (15.0 * flat_stddev)

    dead = np.logical_or(img_fl > hi, img_fl < lo)

    # We also need the original instrument mask (with cracks and such)
    inst_mask = fits.getdata(_inst_mask)

    # Combine into a final supermask. Use the flat file just as a template
    # to get the header from.
    ofile = fits.open(_flat)

    if ((hot.shape)[0] == (inst_mask.shape)[0]):
        mask = hot + dead + inst_mask
    else:
        mask = hot + dead
    mask = (mask != 0)
    unmask = (mask == 0)
    ofile[0].data[unmask] = 0
    ofile[0].data[mask] = 1
    ofile[0].writeto(_out, output_verify='silentfix')

    return


def make_instrument_mask(dark, flat, outDir, instrument=instruments.default_inst):
    """Make the static bad pixel mask for the instrument. This only needs to be
    run once. This creates a file called nirc2mask.fits or osiris_img_mask.fits
    which is subsequently used throughout the pipeline. The dark should be a long
    integration dark.
    
    Parameters
    ----------
    dark : str
        The full absolute path to a medianed dark file. This is
        used to construct a hot pixel mask (4 sigma detection thresh).
    flat : str
        The full absolute path to a medianed flat file. This is
        used to construct a dead pixel mask.
    outDir : str
        full path to output directory with '/' at the end.
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    _out = outDir + instrument.get_bad_pixel_mask_name()
    _dark = dark
    _flat = flat

    util.rmall([_out])

    ##########
    # Make hot pixel mask
    ##########
    # Get the sigma-clipped mean and stddev on the dark
    img_dk = fits.getdata(_dark)
    if parse_version(astropy.__version__) < parse_version('3.0'):
        dark_stats = stats.sigma_clipped_stats(img_dk,
                                               sigma=3,
                                               iters=10)
    else:
        dark_stats = stats.sigma_clipped_stats(img_dk,
                                            sigma=3,
                                               maxiters=10)
    dark_mean = dark_stats[0]
    dark_stddev = dark_stats[2]

    # Clip out the very hot pixels.
    hi = dark_mean + (15.0 * dark_stddev)
    hot = img_dk > hi
    print(('Found %d hot pixels' % (hot.sum())))

    ##########
    # Make dead pixel mask
    ##########
    img_fl = fits.getdata(_flat)
    if parse_version(astropy.__version__) < parse_version('3.0'):
        flat_stats = stats.sigma_clipped_stats(img_dk,
                                               sigma=3,
                                               iters=10)
    else:
    
        flat_stats = stats.sigma_clipped_stats(img_fl,
                                               sigma=3,
                                               maxiters=10)
    flat_mean = flat_stats[0]
    flat_stddev = flat_stats[2]

    # Clip out the dead pixels
    lo = 0.5
    hi = flat_mean + (15.0 * flat_stddev)

    dead = np.logical_or(img_fl > hi, img_fl < lo)
    print(('Found %d dead pixels' % (dead.sum())))

    # Combine into a final supermask
    new_file = fits.open(_flat)

    mask = hot + dead
    mask = (mask != 0)
    unmask = (mask == 0)
    new_file[0].data[unmask] = 0
    new_file[0].data[mask] = 1
    new_file[0].writeto(_out, output_verify='silentfix')


def analyzeDarkCalib(firstFrame, raw_dir=None, reduce_dir=None, skipcombo=False):
    """
    Reduce data from the dark_calib script that should be run once
    a summer in order to test the dark current and readnoise.

    This should be run in the reduce/calib/ directory for a particular
    run.

    Parameters
    ------
    firstFrame : float
        Number of the first frame.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    reduce_dir : str, optional
        Directory such as <epoch>/reduce/ with contents including
        the calib/, calib/darks/, etc. directories live.
        Files will be output into reduce_dir + calib/darks/.
        If epoch_dir is None, then use the current working directory.
    """
    rawDir, redDir = get_raw_reduce_directories(raw_dir, reduce_dir)
    curDir = redDir + 'calib/'
    darkDir = util.trimdir(curDir + 'darks/')

    util.mkdir(curDir)
    util.mkdir(darkDir)

    def printStats(frame, tint, sampmode, reads):
        files = list(range(frame, frame+3))

        fileName = 'dark_%ds_1ca_%d_%dsm.fits' % (tint, sampmode, reads)

        if (skipcombo == False):
            makedark(files, fileName)

        # Get the sigma-clipped mean and stddev on the dark
        img_dk = fits.getdata(darkDir + fileName)
        if parse_version(astropy.__version__) < parse_version('3.0'):
            dark_stats = stats.sigma_clipped_stats(img_dk,
                                                   sigma=3,
                                                   iters=10)
        else:        
            dark_stats = stats.sigma_clipped_stats(img_dk,
                                                   sigma=3,
                                                   maxiters=10)

        darkMean = dark_stats[0]
        darkStdv = dark_stats[2]

        return darkMean, darkStdv

    frame = firstFrame

    lenDarks = 11

    tints = np.zeros(lenDarks) + 12
    tints[-3] = 10
    tints[-2] = 50
    tints[-1] = 100

    reads = np.zeros(lenDarks)
    reads[0] = 1
    reads[1] = 2
    reads[2] = 4
    reads[3] = 8
    reads[4] = 16
    reads[5] = 32
    reads[6] = 64
    reads[7] = 92
    reads[-3:] = 16

    samps = np.zeros(lenDarks) + 3
    samps[0] = 2

    dMeans = np.zeros(lenDarks, dtype=float)
    dStdvs = np.zeros(lenDarks, dtype=float)

    for ii in range(lenDarks):
        (dMeans[ii], dStdvs[ii]) = printStats(frame, tints[ii],samps[ii], reads[ii])
        dStdvs[ii] *= np.sqrt(3)
        frame += 3

    # Calculate the readnoise
    rdnoise = dStdvs * 4.0 * np.sqrt(reads) / (np.sqrt(2.0))
    print(('READNOISE per read: ', rdnoise))


    ##########
    # Print Stuff Out
    ##########
    outFile = darkDir + 'analyzeDarkCalib.out'
    util.rmall([outFile])
    _out = open(outFile,'w')
    hdr = '%8s  %5s  &9s  %9s  %4s  %6s'
    print('Sampmode  Reads  Noise(DN)  Noise(e-)  Tint  Coadds')
    print('--------  -----  ---------  ---------  ----  ------')

    _out.write('Sampmode  Reads  Noise(DN)  Noise(e-)  Tint  Coadds\n')
    _out.write('--------  -----  ---------  ---------  ----  ------\n')

    for ii in range(lenDarks):
        print(('%8d  %5d  %9.1f  %9.1f  %4d  1' % \
               (samps[ii], reads[ii], dStdvs[ii], dStdvs[ii] * 4.0, tints[ii])))
        
    for ii in range(lenDarks):
        _out.write('%8d  %5d  %9.1f  %9.1f  %4d  1\n' % \
                   (samps[ii], reads[ii], dStdvs[ii], dStdvs[ii] * 4.0, tints[ii]))

    _out.close()

    return

def get_raw_reduce_directories(raw_dir, reduce_dir):
    # Determine directory locations
    if reduce_dir is None:
        redDir = os.getcwd() + '/'  # Reduce directory.
    else:
        redDir = util.trimdir(os.path.abspath(reduce_dir) + '/')
    # Set location of raw data
    if raw_dir is not None:
        rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')
    else:
        rawDir = util.trimdir(os.path.abspath(redDir + '../raw') + '/')
    return rawDir, redDir
  
  
