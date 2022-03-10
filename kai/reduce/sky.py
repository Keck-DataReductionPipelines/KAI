from astropy.io import fits
from astropy.table import Table
from astropy import stats
import os, sys, shutil
from kai.reduce import util
#import util
import numpy as np
from pyraf import iraf as ir
from kai import instruments
from datetime import datetime
import pdb
import astropy
from pkg_resources import parse_version

def makesky(files, nite,
            wave, skyscale=True,
            raw_dir=None,
            instrument=instruments.default_inst):
    """
    Make short wavelength (not L-band or longer) skies.
    
    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    nite : str
        Name for night of observation (e.g.: "nite1"), used as suffix
        inside the reduce sub-directories.
    wave : str
        Name for the observation passband (e.g.: "kp")
    skyscale : bool, default=True
        Whether or not to scale the sky files to the common median.
        Turn on for scaling skies before subtraction.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    # Make new directory for the current passband and switch into it
    util.mkdir(wave)
    os.chdir(wave)
    
    # Determine directory locatons
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')
    skyDir = waveDir + 'sky_' + nite + '/'
    
    # Set location of raw data
    rawDir = rootDir + 'raw/'
    
    # Check if user has specified a specific raw directory
    if raw_dir is not None:
        rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')
    
    util.mkdir(skyDir)
    print('sky dir: ',skyDir)
    print('wave dir: ',waveDir)

    skylist = skyDir + 'skies_to_combine.lis'
    output = skyDir + 'sky_' + wave + '.fits'

    util.rmall([skylist, output])

    nn = instrument.make_filenames(files, rootDir=skyDir)
    nsc = instrument.make_filenames(files, rootDir=skyDir, prefix='scale')
    skies = instrument.make_filenames(files, rootDir=rawDir)

    for ii in range(len(nn)):
        if os.path.exists(nn[ii]): os.remove(nn[ii])
        if os.path.exists(nsc[ii]): os.remove(nsc[ii])
        shutil.copy(skies[ii], nn[ii])
    
    # Write out the sources of the sky files
    data_sources_file = open(redDir + 'data_sources.txt', 'a')
    data_sources_file.write('---\n# Sky Files ({0})\n'.format(wave))
    
    for cur_file in skies:
        out_line = '{0} ({1})\n'.format(cur_file, datetime.now())
        data_sources_file.write(out_line)
    
    data_sources_file.close()
    
    # scale skies to common median
    if skyscale:
        _skylog = skyDir + 'sky_scale.log'
        util.rmall([_skylog])
        f_skylog = open(_skylog, 'w')

        sky_mean = np.zeros([len(skies)], dtype=float)

        for i in range(len(skies)):
            # Get the sigma-clipped mean and stddev on the dark
            img_sky = fits.getdata(nn[i], ignore_missing_end=True)
            if parse_version(astropy.__version__) < parse_version('3.0'):
                sky_stats = stats.sigma_clipped_stats(img_sky,
                                                       sigma=3,
                                                       iters=4)
            else:
                sky_stats = stats.sigma_clipped_stats(img_sky,
                                                      sigma=10,
                                                      maxiters=4)

            sky_mean[i] = sky_stats[0]

        sky_all = sky_mean.mean()
        sky_scale = sky_all/sky_mean

        for i in range(len(skies)):
            _nn = fits.open(nn[i], ignore_missing_end=True)
            _nn[0].data = _nn[0].data * sky_scale[i]
            _nn[0].writeto(nsc[i])

            skyf = nn[i].split('/')
            print(('%s   skymean=%10.2f   skyscale=%10.2f' % 
                  (skyf[len(skyf)-1], sky_mean[i],sky_scale[i])))
            f_skylog.write('%s   %10.2f  %10.2f\n' % 
                           (nn[i], sky_mean[i], sky_scale[i]))

        # Make list for combinng
        f_on = open(skylist, 'w')
        f_on.write('\n'.join(nsc) + '\n')
        f_on.close()

        #skylist = skyDir + 'scale????.fits'
        f_skylog.close()
    else:
        # Make list for combinng
        f_on = open(skylist, 'w')
        f_on.write('\n'.join(nn) + '\n')
        f_on.close()

        #skylist = skyDir + 'n????.fits' 

    if os.path.exists(output): os.remove(output)
    ir.unlearn('imcombine')
    ir.imcombine.combine = 'median'
    ir.imcombine.reject = 'none'
    ir.imcombine.nlow = 1
    ir.imcombine.nhigh = 1

    ir.imcombine('@' + skylist, output)
    
    # Change back to original directory
    os.chdir('../')


def makesky_lp(files, nite, wave, number=3, rejectHsigma=None,
               raw_dir=None,
               instrument=instruments.default_inst):
    """
    Make L' skies by carefully treating the ROTPPOSN angle
    of the K-mirror. Uses 3 skies combined (set by number keyword).
    
    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    nite : str
        Name for night of observation (e.g.: "nite1"), used as suffix
        inside the reduce sub-directories.
    wave : str
        Name for the observation passband (e.g.: "lp")
    number : int, default=3
        Number of skies to be combined
    rejectHsigma : int, default:None
        Flag to pass for rejectHsigma for IRAF imcombine
        By default, no flags are passed
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    
    # Make new directory for the current passband and switch into it
    util.mkdir(wave)
    os.chdir(wave)
    
    # Determine directory locatons
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')
    skyDir = waveDir + 'sky_' + nite + '/'
    
    rawDir = rootDir + 'raw/'
    
    # Check if user has specified a specific raw directory
    if raw_dir is not None:
        rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')
    
    util.mkdir(skyDir)
    print('sky dir: ',skyDir)
    print('wave dir: ',waveDir)
    
    raw = instrument.make_filenames(files, rootDir=rawDir)
    skies = instrument.make_filenames(files, rootDir=skyDir)
    
    # Write out the sources of the sky files
    data_sources_file = open(redDir + 'data_sources.txt', 'a')
    data_sources_file.write('---\n# Sky Files ({0})\n'.format(wave))
    
    for cur_file in skies:
        out_line = '{0} ({1})\n'.format(cur_file, datetime.now())
        data_sources_file.write(out_line)
    
    data_sources_file.close()
    
    
    _rawlis = skyDir + 'raw.lis'
    _nlis = skyDir + 'n.lis'
    _skyRot = skyDir + 'skyRot.txt'
    _txt = skyDir + 'rotpposn.txt'
    _out = skyDir + 'sky'
    _log = _out + '.log'
    util.rmall([_rawlis, _nlis, _skyRot, _txt, _out, _log])
    util.rmall([sky + '.fits' for sky in skies])

    open(_rawlis, 'w').write('\n'.join(raw)+'\n')
    open(_nlis, 'w').write('\n'.join(skies)+'\n')

    print('makesky_lp: Getting raw files')
    ir.imcopy('@' + _rawlis, '@' + _nlis, verbose='no')
    ir.hselect('@' + _nlis, "$I,ROTPPOSN", 'yes', Stdout=_skyRot) 

    # Read in the list of files and rotation angles
    files, angles = read_sky_rot_file(_skyRot)

    # Fix angles to be between -180 and 180
    angles[angles > 180] -= 360.0
    angles[angles < -180] += 360.0
    
    sidx = np.argsort(angles)
    
    # Make sorted numarrays
    angles = angles[sidx]
    files = files[sidx]

    f_log = open(_log, 'w')
    f_txt = open(_txt, 'w')
    
    # Skip the first and last since we are going to 
    # average every NN files.
    print('makesky_lp: Combining to make skies.')
    startIdx = number / 2
    stopIdx = len(sidx) - (number / 2)
    for i in range(startIdx, stopIdx):
        sky = 'sky%.1f' % (angles[i])
        skyFits = skyDir + sky + '.fits'
        util.rmall([skyFits])

        # Take NN images
        start = i - (number/2)
        stop = start + number
        list = [file for file in files[start:stop]]
        short = [file for file in files[start:stop]]
        angleTmp = angles[start:stop]

        # Make short names
        for j in range(len(list)):
            tmp = (short[j]).rsplit('/', 1)
            short[j] = tmp[len(tmp)-1]

        print('%s: %s' % (sky, " ".join(short)))
        f_log.write('%s:' % sky)
        for j in range(len(short)):
            f_log.write(' %s' % short[j])
        for j in range(len(angleTmp)):
            f_log.write(' %6.1f' % angleTmp[j])
        f_log.write('\n')

        ir.unlearn('imcombine')
        ir.imcombine.combine = 'median'

        if (rejectHsigma == None):
            ir.imcombine.reject = 'none'
            ir.imcombine.nlow = 1
            ir.imcombine.nhigh = 1
        else:
            ir.imcombine.reject = 'sigclip'
            ir.imcombine.lsigma = 100
            ir.imcombine.hsigma = rejectHsigma
            ir.imcombine.zero = 'median'

        ir.imcombine.logfile = ''
        ir.imcombine(','.join(list), skyFits)

        ir.hedit(skyFits, 'SKYCOMB', 
                 '%s: %s' % (sky, ' '.join(short)), 
                 add='yes', show='no', verify='no')

        f_txt.write('%13s %8.3f\n' % (sky, angles[i]))
	
    f_txt.close()
    f_log.close()
    
    # Change back to original directory
    os.chdir('../')


def makesky_lp2(files, nite, wave):
    """Make L' skies by carefully treating the ROTPPOSN angle
    of the K-mirror. Uses only 2 skies combined."""

    # Start out in something like '06maylgs1/reduce/kp/'
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')
    skyDir = waveDir + 'sky_' + nite + '/'
    rawDir = rootDir + 'raw/'

    util.mkdir(skyDir)

    raw = instrument.make_filenames(files, rootDir=rawDir)
    skies = instrument.make_filenames(files, rootDir=skyDir)
    
    _rawlis = skyDir + 'raw.lis'
    _nlis = skyDir + 'n.lis'
    _skyRot = skyDir + 'skyRot.txt'
    _txt = skyDir + 'rotpposn.txt'
    _out = skyDir + 'sky'
    _log = _out + '.log'
    util.rmall([_rawlis, _nlis, _skyRot, _txt, _out, _log])
    util.rmall([sky + '.fits' for sky in skies])

    open(_rawlis, 'w').write('\n'.join(raw)+'\n')
    open(_nlis, 'w').write('\n'.join(skies)+'\n')

    print('makesky_lp: Getting raw files')
    ir.imcopy('@' + _rawlis, '@' + _nlis, verbose='no')
    ir.hselect('@' + _nlis, "$I,ROTPPOSN", 'yes', Stdout=_skyRot) 

    # Read in the list of files and rotation angles
    files, angles = read_sky_rot_file(_skyRot)

    # Fix angles to be between -180 and 180
    angles[angles > 180] -= 360.0
    angles[angles < -180] += 360.0
    
    sidx = np.argsort(angles)
    
    # Make sorted numarrays
    angles = angles[sidx]
    files = files[sidx]

    f_log = open(_log, 'w')
    f_txt = open(_txt, 'w')

    # Skip the first and last since we are going to 
    # average every 3 files.
    print('makesky_lp: Combining to make skies.')
    for i in range(1, len(sidx)):
        angav = (angles[i] + angles[i-1])/2.
        sky = 'sky%.1f' % (angav)
        skyFits = skyDir + sky + '.fits'
        util.rmall([skyFits])

        # Average 2 images
        list = [file for file in files[i-1:i+1]]
        short = [file for file in files[i-1:i+1]]

        # Make short names
        for j in range(len(list)):
            tmp = (short[j]).rsplit('/', 1)
            short[j] = tmp[len(tmp)-1]

        print('%s: %s %s' % (sky, short[0], short[1]))
        f_log.write('%s: %s %s  %6.1f %6.1f\n' %
                    (sky, short[0], short[1], 
                     angles[i-1], angles[i]))

        ir.unlearn('imcombine')
        ir.imcombine.combine = 'average'
        ir.imcombine.reject = 'none'
        ir.imcombine.nlow = 1
        ir.imcombine.nhigh = 1
        ir.imcombine.logfile = ''
        ir.imcombine(list[1]+','+list[0], skyFits)

        ir.hedit(skyFits, 'SKYCOMB', 
                 '%s: %s %s' % (sky, short[0], short[1]), 
                 add='yes', show='no', verify='no')

        f_txt.write('%13s %8.3f\n' % (sky, angav))
	
    f_txt.close()
    f_log.close()

    #ir.imdelete('@' + _nlis)

def makesky_fromsci(files, nite, wave):
    """Make short wavelength (not L-band or longer) skies."""

    # Start out in something like '06maylgs1/reduce/kp/'
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')
    skyDir = waveDir + 'sky_' + nite + '/'
    rawDir = rootDir + 'raw/'

    util.mkdir(skyDir)
    print('sky dir: ',skyDir)
    print('wave dir: ',waveDir)

    skylist = skyDir + 'skies_to_combine.lis'
    output = skyDir + 'sky_' + wave + '.fits'

    util.rmall([skylist, output])

    nn = instrument.make_filenames(files, rootDir=skyDir)
    nsc = instrument.make_filenames(files, rootDir=skyDir, prefix='scale')
    skies = instrument.make_filenames(files, rootDir=rawDir)

    for ii in range(len(nn)):
        ir.imdelete(nn[ii])
        ir.imdelete(nsc[ii])
        ir.imcopy(skies[ii], nn[ii], verbose="no")

    # Make list for combinng. Reset the skyDir to an IRAF variable.
    ir.set(skydir=skyDir)
    f_on = open(skylist, 'w')
    for ii in range(len(nn)):
        nn_new = nn[ii].replace(skyDir, "skydir$")
        f_on.write(nn_new + '\n')
    f_on.close()

    # Calculate some sky statistics, but reject high (star-like) pixels
    sky_mean = np.zeros([len(skies)], dtype=float)
    sky_std = np.zeros([len(skies)], dtype=float)

    for ii in range(len(nn)):
        img_sky = fits.getdata(nn[i], ignore_missing_end=True)
        if parse_version(astropy.__version__) < parse_version('3.0'):
            sky_stats = stats.sigma_clipped_stats(img_sky,
                                                  sigma_lower=10, sigma_upper=3,
                                                  iters=10)
        else:
            sky_stats = stats.sigma_clipped_stats(img_sky,
                                                  sigma_lower=10, sigma_upper=3,
                                                  maxiters=10)
        sky_mean[ii] = sky_stats[0]
        sky_std[ii] = sky_stats[2]

    sky_mean_all = sky_mean.mean()
    sky_std_all = sky_std.mean()

    # Upper threshold above which we will ignore pixels when combining.
    hthreshold = sky_mean_all + 3.0 * sky_std_all

    ir.imdelete(output)
    ir.unlearn('imcombine')
    ir.imcombine.combine = 'median'
    ir.imcombine.reject = 'sigclip'
    ir.imcombine.mclip = 'yes'
    ir.imcombine.hsigma = 2
    ir.imcombine.lsigma = 10
    ir.imcombine.hthreshold = hthreshold

    ir.imcombine('@' + skylist, output)

def makesky_lp_fromsci(files, nite, wave, number=3, rejectHsigma=None):
    """Make L' skies by carefully treating the ROTPPOSN angle
    of the K-mirror. Uses 3 skies combined (set by number keyword)."""

    # Start out in something like '06maylgs1/reduce/kp/'
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')
    skyDir = waveDir + 'sky_' + nite + '/'
    rawDir = rootDir + 'raw/'

    util.mkdir(skyDir)

    raw = instrument.make_filenames(files, rootDir=rawDir)
    skies = instrument.make_filenames(files, rootDir=skyDir)

    flatDir = redDir + 'calib/flats/'
    flat = flatDir + 'flat_' + wave + '.fits'
    if not os.access(flat, os.F_OK): 
        flat = flatDir + 'flat.fits'
    
    _rawlis = skyDir + 'raw.lis'
    _nlis = skyDir + 'n.lis'
    _skyRot = skyDir + 'skyRot.txt'
    _txt = skyDir + 'rotpposn.txt'
    _out = skyDir + 'sky'
    _log = _out + '.log'
    util.rmall([_rawlis, _nlis, _skyRot, _txt, _out, _log])
    util.rmall([sky + '.fits' for sky in skies])

    open(_rawlis, 'w').write('\n'.join(raw)+'\n')
    open(_nlis, 'w').write('\n'.join(skies)+'\n')

    print('makesky_lp: Getting raw files')
    ir.imarith('@'+_rawlis, '/', flat, '@'+_nlis)
    #ir.imcopy('@' + _rawlis, '@' + _nlis, verbose='no')
    ir.hselect('@' + _nlis, "$I,ROTPPOSN", 'yes', Stdout=_skyRot) 

    # Read in the list of files and rotation angles
    files, angles = read_sky_rot_file(_skyRot)

    # Fix angles to be between -180 and 180
    angles[angles > 180] -= 360.0
    angles[angles < -180] += 360.0
    
    sidx = np.argsort(angles)
    
    # Make sorted numarrays
    angles = angles[sidx]
    files = files[sidx]

    f_log = open(_log, 'w')
    f_txt = open(_txt, 'w')

    # Skip the first and last since we are going to 
    # average every NN files.
    print('makesky_lp: Combining to make skies.')
    startIdx = number / 2
    stopIdx = len(sidx) - (number / 2)
    for i in range(startIdx, stopIdx):
        sky = 'sky%.1f' % (angles[i])
        skyFitsTmp = skyDir + sky + '_tmp.fits'
        skyFits = skyDir + sky + '.fits'
        util.rmall([skyFitsTmp, skyFits])

        # Take NN images
        start = i - (number/2)
        stop = start + number
        list = [file for file in files[start:stop]]
        short = [file for file in files[start:stop]]
        angleTmp = angles[start:stop]

        # Make short names
        for j in range(len(list)):
            tmp = (short[j]).rsplit('/', 1)
            short[j] = tmp[len(tmp)-1]

        print('%s: %s' % (sky, " ".join(short)))
        f_log.write('%s:' % sky)
        for j in range(len(short)):
            f_log.write(' %s' % short[j])
        for j in range(len(angleTmp)):
            f_log.write(' %6.1f' % angleTmp[j])
        f_log.write('\n')

        ir.unlearn('imcombine')
        ir.imcombine.combine = 'median'

        if (rejectHsigma == None):
            ir.imcombine.reject = 'none'
            ir.imcombine.nlow = 1
            ir.imcombine.nhigh = 1
        else:
            ir.imcombine.reject = 'sigclip'
            ir.imcombine.lsigma = 100
            ir.imcombine.hsigma = rejectHsigma
            ir.imcombine.zero = 'median'

        ir.imcombine.logfile = ''
        ir.imcombine(','.join(list), skyFitsTmp)

        ir.imarith(skyFitsTmp, '*', flat, skyFits)

        ir.hedit(skyFits, 'SKYCOMB', 
                 '%s: %s' % (sky, ' '.join(short)), 
                 add='yes', show='no', verify='no')

        f_txt.write('%13s %8.3f\n' % (sky, angles[i]))
	
    f_txt.close()
    f_log.close()


def read_sky_rot_file(sky_rot_file):
    """Read in the list of files and rotation angles."""
    
    rotTab = Table.read(sky_rot_file, format='ascii', header_start=None)
    cols = list(rotTab.columns.keys())
    files = rotTab[cols[0]]
    angles = rotTab[cols[1]]

    return files, angles
    
