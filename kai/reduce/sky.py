from astropy.io import fits
from astropy.table import Table
from astropy import stats
import os, shutil
from kai.reduce import util, lin_correction
import numpy as np
from kai import instruments
from datetime import datetime
import pdb
import astropy
import warnings
from pkg_resources import parse_version


def makesky(files, nite, wave,
            dark_frame=None, skyscale=True,
            raw_dir=None, reduce_dir=None,
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
    dark_frame : str, default=None
        File name for the dark frame in order to carry out dark correction.
        If not provided, dark frame is not subtracted and a warning is thrown.
        Assumes dark file is located under ./calib/darks/
    skyscale : bool, default=True
        Whether or not to scale the sky files to the common median.
        Turn on for scaling skies before subtraction.
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

    waveDir = redDir + wave + '/'
    skyDir = waveDir + 'sky_' + nite + '/'

    # Make new directory for the current passband and switch into it
    util.mkdir(wave)
    util.mkdir(skyDir)
    print('raw dir: ', rawDir)
    print('sky dir: ', skyDir)
    print('wave dir: ', waveDir)

    skylist = skyDir + 'skies_to_combine.lis'
    output = skyDir + 'sky_' + wave + '.fits'

    util.rmall([skylist, output])

    nn = instrument.make_filenames(files, rootDir=skyDir)  # copy of raw sky
    nsc = instrument.make_filenames(files, rootDir=skyDir, prefix='scale')  # scaled sky
    skies = instrument.make_filenames(files, rootDir=rawDir)  # original raw sky

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

    # If dark frame is provided, carry out dark correction
    if dark_frame is not None:
        dark_file = redDir + '/calib/darks/' + dark_frame

        # Read in dark frame data
        dark_data = fits.getdata(dark_file, ignore_missing_end=True)

        # Go through each sky file
        for i in range(len(skies)):
            with fits.open(nn[i], mode='readonly', output_verify='ignore', ignore_missing_end=True) as cur_sky:
                sky_data = cur_sky[0].data
                sky_header = cur_sky[0].header
            sky_data = sky_data - dark_data
            sky_hdu = fits.PrimaryHDU(data=sky_data, header=sky_header)
            sky_hdu.writeto(nn[i], output_verify='ignore', overwrite=True)
    else:
        warning_message = 'Dark frame not provided for makesky().'
        warning_message += '\nUsing sky frames without dark subtraction.'

        warnings.warn(warning_message)

    # Perform linearity correction
    for i in range(len(skies)):
        lin_correction.lin_correction(nn[i], instrument=instrument)

    # List of skies to combine (might be changed after scaling).
    skies_to_combine = nn

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
        sky_scale = sky_all / sky_mean

        for i in range(len(skies)):
            _nn = fits.open(nn[i], ignore_missing_end=True)
            _nn[0].data = _nn[0].data * sky_scale[i]
            _nn[0].writeto(nsc[i])

            skyf = nn[i].split('/')
            print(('%s   skymean=%10.2f   skyscale=%10.2f' %
                   (skyf[len(skyf) - 1], sky_mean[i], sky_scale[i])))
            f_skylog.write('%s   %10.2f  %10.2f\n' %
                           (nn[i], sky_mean[i], sky_scale[i]))

        # Make list for combinng
        f_on = open(skylist, 'w')
        f_on.write('\n'.join(nsc) + '\n')
        f_on.close()

        # skylist = skyDir + 'scale????.fits'
        f_skylog.close()
        skies_to_combine = nsc
    else:
        # Make list for combinng
        f_on = open(skylist, 'w')
        f_on.write('\n'.join(nn) + '\n')
        f_on.close()

        # skylist = skyDir + 'n????.fits'

    # Combine the skies
    img_stack = []  # stack of image data
    hdr_stack = []  # stack of image headers
    for ii in range(len(skies_to_combine)):
        img, hdr = fits.getdata(skies_to_combine[ii], header=True)
        img_stack.append(img)
        hdr_stack.append(hdr)
    img_stack = np.array(img_stack)

    # Stack the darks with sigma clipping, median combine.
    sk_avg, sk_med, sk_std = astropy.stats.sigma_clipped_stats(img_stack,
                                                               cenfunc='median',
                                                               sigma_lower=3,
                                                               sigma_upper=3,
                                                               axis=0)

    # Save to an output file. Use first header.
    fits.writeto(output, sk_med, header=hdr_stack[0], overwrite=True)

    return


def makesky_lp(files, nite, wave,
               dark_frame=None, number=3, rejectHsigma=5,
               raw_dir=None, reduce_dir=None,
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
    dark_frame : str, default=None
        File name for the dark frame in order to carry out dark correction.
        If not provided, dark frame is not subtracted and a warning is thrown.
        Assumes dark file is located under ./calib/darks/
    number : int, default=3
        Number of skies to be combined
    rejectHsigma : int, default:None
        Apply a sigma_high threshold for sigma clipping. This works to
        remove stars or galaxies when using science images to make the sky.
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

    waveDir = redDir + wave + '/'
    skyDir = waveDir + 'sky_' + nite + '/'

    # Make new directory for the current passband and switch into it
    util.mkdir(waveDir)
    os.chdir(waveDir)

    util.mkdir(skyDir)
    print('raw dir: ', rawDir)
    print('sky dir: ', skyDir)
    print('wave dir: ', waveDir)

    # Copy over the raw sky files into the sky directory
    skies = instrument.make_filenames(files, rootDir=skyDir)
    raw = instrument.make_filenames(files, rootDir=rawDir)

    for ii in range(len(skies)):
        if os.path.exists(skies[ii]): os.remove(skies[ii])
        shutil.copy(raw[ii], skies[ii])

    # Write out the sources of the sky files
    data_sources_file = open(redDir + 'data_sources.txt', 'a')
    data_sources_file.write('---\n# Sky Files ({0})\n'.format(wave))

    for cur_file in raw:
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

    open(_rawlis, 'w').write('\n'.join(raw) + '\n')
    open(_nlis, 'w').write('\n'.join(skies) + '\n')

    print('makesky_lp: Getting raw files')

    write_sky_rot_file(_rawlis, _nlis, _skyRot)

    # If dark frame is provided, carry out dark correction
    if dark_frame is not None:
        dark_file = redDir + '/calib/darks/' + dark_frame

        # Read in dark frame data
        dark_data = fits.getdata(dark_file, ignore_missing_end=True)

        # Go through each sky file
        for i in range(len(skies)):
            with fits.open(skies[i], mode='readonly', output_verify='ignore', ignore_missing_end=True) as cur_sky:
                sky_data = cur_sky[0].data
                sky_header = cur_sky[0].header
            sky_data = sky_data - dark_data
            sky_hdu = fits.PrimaryHDU(data=sky_data, header=sky_header)
            sky_hdu.writeto(skies[i], output_verify='ignore', overwrite=True)
    else:
        warning_message = 'Dark frame not provided for makesky_lp().'
        warning_message += '\nUsing sky frames without dark subtraction.'

        warnings.warn(warning_message)

    # Perform linearity correction
    for i in range(len(skies)):
        lin_correction.lin_correction(skies[i], instrument=instrument)

    # Read in the list of files and rotation angles
    files, angles = read_sky_rot_file(_skyRot)

    # Fix angles to be between -180 and 180
    angles[angles > 180] -= 360.0
    angles[angles < -180] += 360.0

    sidx = np.argsort(angles)

    # Make sorted numarrays
    angles = angles[sidx]
    files = files[sidx]

    # Open some logging files
    f_log = open(_log, 'w')
    f_txt = open(_txt, 'w')

    # Skip the first and last since we are going to 
    # average every NN files.
    print('makesky_lp: Combining to make skies.')
    startIdx = number / 2
    stopIdx = len(sidx) - (number / 2)
    for i in np.arange(startIdx, stopIdx):
        # Take NN images
        start = int(i - (number / 2))
        stop = int(start + number)

        list = [file for file in files[start:stop]]
        short = [file for file in files[start:stop]]
        angleTmp = angles[start:stop]
        angle_avg = angleTmp.mean()

        sky = 'sky%.1f' % (angle_avg)
        skyFits = skyDir + sky + '.fits'
        util.rmall([skyFits])

        # Make short names
        for j in range(len(list)):
            tmp = (short[j]).rsplit('/', 1)
            short[j] = tmp[len(tmp) - 1]

        # Load up the FITS files.
        img_stack = []
        hdr_stack = []
        for j in range(len(list)):
            img, hdr = fits.getdata(list[j], header=True)
            img_stack.append(img)
            hdr_stack.append(hdr)
        img_stack = np.array(img_stack)

        # Log the short names of the sky frames in this stack.
        print('%s: %s' % (sky, " ".join(short)))
        f_log.write('%s:' % sky)
        for j in range(len(short)):
            f_log.write(' %s' % short[j])
        for j in range(len(angleTmp)):
            f_log.write(' %6.1f' % angleTmp[j])
        f_log.write('\n')

        # Decide if we are using mean or median.
        if number < 3:
            cenfunc = 'mean'
        else:
            cenfunc = 'median'

        # Stack the images.
        sk_avg, sk_med, sk_std = astropy.stats.sigma_clipped_stats(img_stack,
                                                                   cenfunc=cenfunc,
                                                                   sigma_lower=100,
                                                                   sigma_upper=rejectHsigma,
                                                                   axis=0)
        if number < 3:
            sk_img = sk_avg
        else:
            sk_img = sk_med

        sk_hdr = hdr_stack[0]
        sk_hdr['SKYCOMB'] = '%s: %s' % (sky, ' '.join(short))
        fits.writeto(skyFits, sk_img, header=sk_hdr, overwrite=True)

        # Save output sky angle to txt file.
        f_txt.write('%13s %8.3f\n' % (sky, angle_avg))

    f_txt.close()
    f_log.close()

    # Change back to original directory
    os.chdir('../')

    return


def makesky_fromsci(files, nite, wave,
                    raw_dir=None, reduce_dir=None,
                    instrument=instruments.default_inst):
    """
    Make short wavelength (not L-band or longer) skies from the science exposures
    themselves. We use strong clipping to get rid of stars. This should only be used
    on sparser science fields.

    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    nite : str
        Name for night of observation (e.g.: "nite1"), used as suffix
        inside the reduce sub-directories.
    wave : str
        Name for the observation passband (e.g.: "kp")
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

    util.mkdir(wave)
    os.chdir(wave)

    waveDir = redDir + wave + '/'
    skyDir = waveDir + 'sky_' + nite + '/'

    # Make new directory for the current passband and switch into it
    util.mkdir(wave)
    os.chdir(wave)

    util.mkdir(skyDir)
    print('raw dir: ', rawDir)
    print('sky dir: ', skyDir)
    print('wave dir: ', waveDir)

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

    # Perform linearity correction
    for i in range(len(skies)):
        lin_correction.lin_correction(nn[i], instrument=instrument)

    # Calculate some sky statistics, but reject high (star-like) pixels
    sky_mean = np.zeros([len(skies)], dtype=float)
    sky_std = np.zeros([len(skies)], dtype=float)

    for ii in range(len(nn)):
        img_sky = fits.getdata(nn[ii], ignore_missing_end=True)
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

    # Combine the skies
    img_stack = []  # stack of image data
    hdr_stack = []  # stack of image headers
    for ii in range(len(skies)):
        img, hdr = fits.getdata(skies[ii], header=True)
        img_stack.append(img)
        hdr_stack.append(hdr)
    img_stack = np.array(img_stack)

    # Mask all pixels greater than the high threshold.
    img_stack_msk = np.ma.masked_greater(img_stack, hthreshold)

    # Stack the darks with sigma clipping, median combine.
    sk_avg, sk_med, sk_std = astropy.stats.sigma_clipped_stats(img_stack_msk,
                                                               cenfunc='median',
                                                               sigma_lower=10,
                                                               sigma_upper=2,
                                                               axis=0)

    # Save to an output file. Use first header.
    fits.writeto(output, sk_med, header=hdr_stack[0], overwrite=True)

    # Change back to original directory
    os.chdir('../')

    return


def read_sky_rot_file(sky_rot_file):
    """Read in the list of files and rotation angles."""
    rotTab = Table.read(sky_rot_file, format='ascii', header_start=None)
    cols = list(rotTab.columns.keys())
    files = rotTab[cols[0]]
    angles = rotTab[cols[1]]

    return files, angles


def write_sky_rot_file(rawlis, nlis, skyRot):
    """Write the list of files and rotation angles in Lp.
    Created to avoid using pyraf which causes issues when dealing
    with non-fits conforming headers."""
    # 1. Copy file
    shutil.copyfile(rawlis, nlis)
    # 2. Loop through files to obtain ROTPPOSN header and append to list
    with open(nlis) as file:
        contents = file.read().split('\n')
    rot = []
    for fit in contents:
        try:
            with fits.open(fit, mode='readonly', output_verify='ignore',
                           ignore_missing_end=True) as rot_hdu:
                rot_header = rot_hdu[0].header
                rot.append(rot_header['ROTPPOSN'])
        except:
            continue
    # 3. Write in skyRot the name of the file and the header values
    with open(skyRot, 'at') as edit:
        for jj, ii in enumerate(rot):
            edit.write(contents[jj] + '\t' + str(ii) + '\n')

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

