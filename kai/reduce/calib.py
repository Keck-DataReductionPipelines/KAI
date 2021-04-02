import os, sys
from . import util
from astropy.io import fits
from astropy import stats
from pyraf import iraf as ir
from kai import instruments
import numpy as np
from astropy import stats
import astropy
from datetime import datetime
from pkg_resources import parse_version

module_dir = os.path.dirname(__file__)

def makedark(files, output,
             raw_dir=None,
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
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    redDir = os.getcwd() + '/'  # Reduce directory.
    curDir = redDir + 'calib/'
    darkDir = util.trimdir(curDir + 'darks/')
    
    # Set location of raw data
    rawDir = util.trimdir(os.path.abspath(redDir + '../raw') + '/')
    
    # Check if user has specified a specific raw directory
    if raw_dir is not None:
        rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')
    
    util.mkdir(curDir)
    util.mkdir(darkDir)
    
    _out = darkDir + output
    _outlis = darkDir + 'dark.lis'
    util.rmall([_out, _outlis])

    darks = instrument.make_filenames(files, rootDir=rawDir)
    
    # Write out the sources of the dark files
    data_sources_file = open(redDir + 'data_sources.txt', 'a')
    data_sources_file.write(
        '---\n# Dark Files for {0} \n'.format(output))
    
    for cur_file in darks:
        out_line = '{0} ({1})\n'.format(cur_file, datetime.now())
        data_sources_file.write(out_line)
    
    data_sources_file.close()
    
    f_on = open(_outlis, 'w')
    f_on.write('\n'.join(darks) + '\n')
    f_on.close()

    ir.unlearn('imcombine')
    ir.imcombine.combine = 'median'
    ir.imcombine.reject = 'sigclip'
    ir.imcombine.nlow = 1
    ir.imcombine.nhigh = 1
    ir.imcombine('@' + _outlis, _out)


def makeflat(onFiles, offFiles, output, normalizeFirst=False,
             raw_dir=None,
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
    normalizeFirst : bool, default=False
        If the individual flats should be normalized first,
        such as in the case of twilight flats.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    redDir = os.getcwd() + '/'
    curDir = redDir + 'calib/'
    flatDir = util.trimdir(curDir + 'flats/')
    
    # Set location of raw data
    rawDir = util.trimdir(os.path.abspath(redDir + '../raw') + '/')
    
    # Check if user has specified a specific raw directory
    if raw_dir is not None:
        rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')
    
    util.mkdir(curDir)
    util.mkdir(flatDir)

    _on = flatDir + 'lampsOn.fits'
    _off = flatDir + 'lampsOff.fits'
    _norm = flatDir + 'flatNotNorm.fits'
    _out = flatDir + output
    _onlis = flatDir + 'on.lis'
    _offlis = flatDir + 'off.lis'
    _onNormLis = flatDir + 'onNorm.lis'

    util.rmall([_on, _off, _norm, _out, _onlis, _offlis, _onNormLis])

    lampson = instrument.make_filenames(onFiles, rootDir=rawDir)
    lampsoff = instrument.make_filenames(offFiles, rootDir=rawDir)
    lampsonNorm = instrument.make_filenames(onFiles, rootDir=flatDir + 'norm')
    util.rmall(lampsonNorm)
    
    # Write out the sources of the dark files
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
    
    if (len(offFiles) != 0):
        f_on = open(_onlis, 'w')
        f_on.write('\n'.join(lampson) + '\n')
        f_on.close()
        f_on = open(_offlis, 'w')
        f_on.write('\n'.join(lampsoff) + '\n')
        f_on.close()
        f_onn = open(_onNormLis, 'w')
        f_onn.write('\n'.join(lampsonNorm) + '\n')
        f_onn.close()

        # Combine to make a lamps on and lamps off
        ir.unlearn('imcombine')
        ir.imcombine.combine = 'median'
        ir.imcombine.reject = 'sigclip'
        ir.imcombine.nlow = 1
        ir.imcombine.nhigh = 1
        ir.imcombine('@' + _offlis, _off)

        # Check if we should normalize individual flats first
        # such as in the case of twilight flats.
        if normalizeFirst:
            f_on = open(_offlis, 'w')
            f_on.write('\n'.join(lampsoff) + '\n')
            f_on.close()

            # Subtract "off" from individual frames
            ir.imarith('@'+_onlis, '-', _off, '@'+_onNormLis)

            # Scale them and combine
            ir.imcombine.scale = 'median'
            ir.imcombine('@' + _onNormLis, _norm)
        else:
            # Combine all "on" frames
            ir.imcombine('@' + _onlis, _on)

            # Now do lampsOn - lampsOff
            ir.imarith(_on, '-', _off, _norm)


        # Normalize the final flat
        ir.module.load('noao', doprint=0, hush=1)
        ir.module.load('imred', doprint=0, hush=1)
        ir.module.load('generic', doprint=0, hush=1)
        orig_img = fits.getdata(_norm)
        orig_size = (orig_img.shape)[0]
        if (orig_size >= 1024):
            flatRegion = '[100:900,513:950]'
        else:
            flatRegion = ''
        ir.normflat(_norm, _out, sample=flatRegion)

    else:
        f_on = open(_onlis, 'w')
        f_on.write('\n'.join(lampson) + '\n')
        f_on.close()

        # Combine twilight flats
        ir.unlearn('imcombine')
        ir.imcombine.combine = 'median'
        ir.imcombine.reject = 'sigclip'
        ir.imcombine.nlow = 1
        ir.imcombine.nhigh = 1
        if normalizeFirst:
            # Scale them
            ir.imcombine.scale = 'median'
        ir.imcombine('@' + _onlis, _norm)

        # Normalize the flat
        ir.module.load('noao', doprint=0, hush=1)
        ir.module.load('imred', doprint=0, hush=1)
        ir.module.load('generic', doprint=0, hush=1)
        flatRegion = '[100:900,513:950]'
        ir.normflat(_norm, _out, sample=flatRegion)

def makemask(dark, flat, output, instrument=instruments.default_inst):
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
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    redDir = os.getcwd() + '/'
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

def analyzeDarkCalib(firstFrame, skipcombo=False):
    """
    Reduce data from the dark_calib script that should be run once
    a summer in order to test the dark current and readnoise.

    This should be run in the reduce/calib/ directory for a particular
    run.
    """
    redDir = os.getcwd() + '/'  # Reduce directory.
    curDir = redDir + 'calib/'
    darkDir = util.trimdir(curDir + 'darks/')
    rawDir = util.trimdir(os.path.abspath(redDir + '../raw') + '/')

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
