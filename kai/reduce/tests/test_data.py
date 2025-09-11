import unittest
import pytest
import shutil
import os
from kai import instruments
from kai.reduce import data
from kai.reduce import bfixpix
from kai.reduce import util
from astropy.io import fits
import math
import numpy as np
from astropy.table import Table
from astropy import stats
from copy import deepcopy


def test_clean():
    nirc2 = instruments.NIRC2()
    
    mod_path = os.path.dirname(os.path.abspath(data.__file__))
    
    epoch_dir = mod_path + '/../data/test_epoch/17may21/'
    reduce_dir = epoch_dir + 'reduce/'
    clean_dir = epoch_dir + 'clean/'
    raw_dir = epoch_dir + 'raw/'
    
    dark_files = range(214, 223 + 1)
    dark_out = 'dark_30.0s_1ca.fits'
    dark_out_full_path = reduce_dir + 'calib/darks/' + dark_out  # only use for testing

    os.chdir(reduce_dir)
    
    sci_files = ['{0:03d}'.format(ii) for ii in range(7, 19+1)]
    refSrc = [556, 564]
    data.clean(sci_files, 'ob170095', 'kp', refSrc, refSrc, field='ob170095',
               raw_dir=raw_dir, clean_dir=clean_dir, instrument=nirc2)

    assert os.path.exists(clean_dir + 'ob170095_kp/c0007.fits') == True
    assert os.path.exists(clean_dir + 'ob170095_kp/c0019.fits') == True
    
    return


def test_combine():
    nirc2 = instruments.NIRC2()

    mod_path = os.path.dirname(os.path.abspath(data.__file__))
    
    epoch_dir = mod_path + '/../data/test_epoch/17may21/'
    reduce_dir = epoch_dir + 'reduce/'
    combo_dir = epoch_dir + 'combo/'

    os.chdir(reduce_dir)

    sci_files = ['{0:03d}'.format(ii) for ii in range(7, 19+1)]
    
    data.combine(sci_files, 'kp', '17may21', trim=1,
             weight='strehl', field='ob170095', submaps=3)

    assert os.path.exists(combo_dir + 'mag17may21_ob170095_kp.fits') == True
    assert os.path.exists(combo_dir + 'm17may21_ob170095_kp_3.fits') == True

    return


def test_rot_img():
    mod_path = os.path.dirname(os.path.abspath(data.__file__))

    epoch_dir = mod_path + '/../data/test_epoch/17may21/'
    clean_dir = epoch_dir + 'clean/ob170095_kp/'

    img_root = '0007'
    phi = 30.0   # arbitrary.

    rot_img_filename = clean_dir + 'r' + img_root + '.fits'

    data.rot_img(img_root, phi, cleanDir=clean_dir)

    rot_img_filename = clean_dir + 'r' + img_root + '.fits'

    img_rot, img_rot_hdr = fits.getdata(rot_img_filename, header = True)

    assert img_rot.shape == (1024, 1024)

    # Check various updated headers (1e-6 somewhat arbitrary)
    assert math.isclose(img_rot_hdr['CRPIX1'], 512.4051640853286, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['CRPIX2'], 512.3818982131779, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['LTV1'], 324.5949932623673, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['LTV2'], -187.40500673763256, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['LTM1_1'], 0.8660254037844387, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['LTM1_2'], 0.49999999999999994, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['LTM2_1'], -0.49999999999999994, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['LTM2_2'], 0.8660254037844387, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['CD1_1'], -2.4001372467297e-06, abs_tol = 1e-12)
    assert math.isclose(img_rot_hdr['CD1_2'], -1.3716803435255e-06, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['CD2_1'], -1.3716803435255e-06, abs_tol = 1e-6)
    assert math.isclose(img_rot_hdr['CD2_2'], 2.40013724672977e-06, abs_tol = 1e-6)

    return

def test_run_clean_cosmicrays():
    mod_path = os.path.dirname(os.path.abspath(data.__file__))

    epoch_dir = mod_path + '/../data/test_epoch/17may21/'
    reduce_dir = epoch_dir + 'reduce/kp/sci_ob170095/'

    #for ii in range(7, 20):
    img_root = '0007' #'{0:04d}'.format(ii)

    to_clean_img_filename = reduce_dir + 'ff' + img_root + '.fits'
    clean_img_filename = reduce_dir + 'ff' + img_root + '_f.fits'
    
    crmask_filename = reduce_dir + 'crmask' + img_root + '.fits'
    _statmask = reduce_dir + 'stat_mask' + img_root + '.fits'
    _ff_s = reduce_dir + 'ff' + img_root + '_s.fits'

    _supermask = epoch_dir + 'reduce/calib/masks/supermask.fits'

    util.rmall([crmask_filename])

    ### Fix bad pixels ###
    # Produces _ff_f file
    bfixpix.bfixpix(to_clean_img_filename, _statmask)
    util.rmall([_ff_s])

    data.clean_cosmicrays(_ff=clean_img_filename, _mask=crmask_filename, wave='kp', _input_mask = _supermask)

    crmask = fits.getdata(reduce_dir + 'crmask{}.fits'.format(img_root))
    img = fits.getdata(reduce_dir + 'ff{}_f.fits'.format(img_root))

    total_cosmic_rays = np.sum(crmask)

    assert total_cosmic_rays == 40

    replaced_vals = img[np.where(crmask == True)]
    reference_replaced_vals = np.array([ 93.22003923, 109.46820358, 104.71930535, 108.7328641 ,
                                       114.35664048, 128.26219344, 120.95855976, 121.0454987 ,
                                       121.88991271, 109.69561369, 111.76535249, 103.19907258,
                                       119.03382031, 112.67812526, 105.74626037, 103.82823455,
                                       100.79214175, 112.74748025, 110.15684462, 120.08868667,
                                       117.81130218, 123.66382118, 120.7714691 , 107.01930555,
                                       104.82423279, 102.46210371,  97.71296729, 104.68712229,
                                       102.99442105, 139.74373983, 106.26805862, 121.44722921,
                                       120.31932389, 121.95712798, 102.49784852, 119.8155971 ,
                                       103.76733508, 126.02076741, 127.17467951, 126.54006193])
    replaced_diffs = np.abs(replaced_vals - reference_replaced_vals)
    assert np.all(replaced_diffs < 2)

    return

def test_run_combine_register():
    mod_path = os.path.dirname(os.path.abspath(data.__file__))

    epoch_dir = mod_path + '/../data/test_epoch/17may21/'
    cleanDir = epoch_dir + 'clean/ob170095_kp/'
    comboDir = epoch_dir + 'combo/'
    
    outroot = '17may21_ob170095'
    wave = 'kp'
    instrument = instruments.default_inst

    files = ['n{0:04d}'.format(ii) for ii in range(7, 19+1)]
    # Make strings out of all the filename roots.
    roots = instrument.make_filenames(files, prefix='')
    #roots = [aa.replace('.fits', '') for aa in roots]
    rootsc = [aa.replace('n', 'c') for aa in roots]
    
    _out = comboDir + 'mag' + outroot + '_' + wave
    
    strehls, fwhm = data.loadStrehl(cleanDir, rootsc)
    # skip trimming since no trimming in this test

    roots = [aa.replace('.fits', '') for aa in roots]
    roots = [aa.replace('n', '') for aa in roots]

    # Determine the reference image
    # refImage_index = 0    # Use the first image from night
    refImage_index = np.argmin(fwhm)    # Use the lowest FWHM frame
    refImage = cleanDir + 'c' + roots[refImage_index] + '.fits'
    print('combine: reference image - %s' % refImage)

    # See if all images are at same PA, if not, rotate all to PA = 0
    # temporarily. This needs to be done to get correct shifts.
    diffPA = data.combine_rotation(cleanDir, roots, instrument=instrument)

    # Make a table of coordinates for the reference source.
    # These serve as initial estimates for the shifts.
    data.combine_coo(_out + '.coo', cleanDir, roots, diffPA, refImage_index)

    # Keep record of files that went into this combine
    data.combine_lis(_out + '.lis', cleanDir, roots, diffPA)

    # Register images to get shifts.
    shiftsTab = data.combine_register(_out, refImage, diffPA)

    reference_shifts = np.array([('c0007.fits', -9.84464849e+00,  1.49251340e+01),
                               ('c0008.fits',  2.78008820e+01,  2.65070686e+01),
                               ('c0009.fits',  2.76459768e+01,  2.64722525e+01),
                               ('c0010.fits',  2.78401224e+01,  2.61418346e+01),
                               ('c0011.fits',  2.79354609e+01,  2.62435264e+01),
                               ('c0012.fits',  2.73448354e+01,  2.67489579e+01),
                               ('c0013.fits', -5.43577967e-01, -5.16297914e-01),
                               ('c0014.fits', -6.69372700e-01, -3.53712142e-01),
                               ('c0015.fits', -6.05979081e-01, -2.79023924e-01),
                               ('c0016.fits', -4.10188242e-03,  1.03180342e-03),
                               ('c0017.fits', -3.78230707e-01,  2.72357711e-01),
                               ('c0018.fits',  6.12677200e+01,  2.98161363e+01),
                               ('c0019.fits',  6.14623497e+01,  2.96408862e+01)],
                              dtype=[('col0', '<U10'), ('col1', '<f8'), ('col2', '<f8')])

    assert np.all(np.abs(np.array(shiftsTab['col1']) - np.array(reference_shifts['col1'])) < 1e-2)
    assert np.all(np.abs(np.array(shiftsTab['col2']) - np.array(reference_shifts['col2'])) < 1e-2)


    return


def test_run_clean_drizzle():
    # To be run in the reduce directory of the test_epoch
    mod_path = os.path.dirname(os.path.abspath(data.__file__))

    epoch_dir = mod_path + '/../data/test_epoch/17may21_noiraf/'
    rawDir = epoch_dir + 'raw/'
    
    instrument = instruments.NIRC2()
    files = ['{0:04d}'.format(ii) for ii in range(7, 19+1)]
    firstFile = instrument.make_filenames([files[0]], rootDir=rawDir)[0]
    hdr1 = fits.getheader(firstFile, ignore_missing_end=True)

    
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
    #data.setup_drizzle(imgsize)
    
    os.chdir(epoch_dir + 'reduce/kp/sci_ob170095/')
    for f in [files[0]]:
        _bp = instrument.make_filenames([f], prefix='bp')[0]
        _ce = instrument.make_filenames([f], prefix='ce')[0]
        _wgt = instrument.make_filenames([f], prefix='wgt')[0]
        _ff = instrument.make_filenames([f], prefix='ff')[0]
        _ff_f = _ff.replace('_.fits', '_f.fits')
        #_ff_f = _ff.replace('_noiraf.fits', '_f_noiraf.fits')
        _dlog_tmp = instrument.make_filenames([f], prefix='driz')[0]
        _dlog = _dlog_tmp.replace('.fits', '.log')

        ### Background Subtraction ###
        if os.path.exists(_bp) == False:
            bkg = data.clean_bkgsubtract(_ff_f, _bp)
        else:
            os.remove(_bp)
            bkg = data.clean_bkgsubtract(_ff_f, _bp)
        
        ### Drizzle individual file ###
        data.clean_drizzle(distXgeoim, distYgeoim, _bp, _ce, _wgt, _dlog,
                      fixDAR=True, instrument=instrument,
                      use_koa_weather=False)

    drizzled_file = fits.open('ce0007.fits')
    weight_file = fits.open('wgt0007.fits')

    # check shape is correct
    assert np.array_equal(np.shape(drizzled_file[0].data), np.array([1024, 1024]))
    assert np.array_equal(np.shape(weight_file[0].data), np.array([1024, 1024]))

    # makes sure there's no nans in the drizzle file
    assert np.sum(np.isnan(drizzled_file[0].data)) == 0

    # Check header vals
    hdr = drizzled_file[0].header
    assert hdr['D001DATA'] == 'bp0007.fits'
    assert hdr['D001DEXP'] == 30.0
    assert hdr['D001OUDA'] == 'ce0007.fits'
    assert hdr['D001OUWE'] == 'wgt0007.fits'
    assert hdr['D001OUCO'] == ''
    assert hdr['D001MASK'] == ''
    assert hdr['D001WTSC'] == 1.0
    assert hdr['D001KERN'] == 'lanczos3'
    assert hdr['D001PIXF'] == 1.0
    assert hdr['D001COEF'] == ''
    assert hdr['D001XGIM'] == 'ce0007geo_x.fits'
    assert hdr['D001YGIM'] == 'ce0007geo_y.fits'
    assert hdr['D001SCAL'] == 1
    assert hdr['D001ROT'] == 0
    assert hdr['D001XSH'] == 0
    assert hdr['D001YSH'] == 0
    assert hdr['D001SFTU'] == 'pixels'
    assert hdr['D001SFTF'] == 'pixels'
    assert hdr['D001EXKY'] == 'ITIME'
    assert hdr['D001INUN'] == 'counts'
    assert hdr['D001OUUN'] == 'counts'
    assert hdr['D001FVAL'] == '0'

    os.chdir('../..')

    return
        

def test_run_combine_drizzle():
    # To be run in the reduce directory of the test_epoch
    mod_path = os.path.dirname(os.path.abspath(data.__file__))
    epoch_dir = mod_path + '/../data/test_epoch/17may21/'
    cleanDir = epoch_dir + 'clean/ob170095_kp/'
    combo_dir = epoch_dir + 'combo/'

    _out = epoch_dir + '/combo/mag17may21_ob170095_kp'
    _sub = epoch_dir + '/combo/m17may21_ob170095_kp'
    refImage = epoch_dir + 'clean/ob170095_kp/c0016.fits'
    diffPA = 0
    submaps = 3
    wave = 'kp'
    instrument = instruments.NIRC2()

    roots = ['0007', '0008', '0009', '0010', '0011', '0012', '0013', 
             '0014', '0015', '0016', '0017', '0018', '0019']
    strehls = np.array([0.314, 0.369, 0.348, 0.376, 0.39, 0.373, 
                        0.347, 0.371, 0.357, 0.388, 0.368, 0.326, 0.333])
    fwhm = np.array([65.85, 58.66, 60.65, 57.72, 56.37, 57.86, 
                     57.84, 55.89, 58.17, 55.78, 56.03, 60.95, 59.87])
    weights = np.array([0.06738197424892704, 0.07918454935622317, 0.07467811158798282, 
                        0.08068669527896996, 0.08369098712446352, 0.08004291845493562, 
                        0.07446351931330471,  0.0796137339055794, 0.07660944206008583, 
                        0.08326180257510729, 0.07896995708154506, 0.06995708154506437, 0.07145922746781116])
    

    shiftsTab = Table.read(combo_dir + 'mag17may21_ob170095_kp.shifts', format='ascii', data_start=1)
    roots, strehls, fwhm, weights, shiftsTab = data.sort_frames(roots, strehls, fwhm, weights, shiftsTab)
    xysize = 1170

    data.combine_drizzle(xysize, cleanDir, roots, _out, weights, shiftsTab,
                wave, diffPA, fixDAR=True, mask=True, instrument=instrument,
                use_koa_weather=False)

    drizzled_img = fits.open(combo_dir + 'mag17may21_ob170095_kp.fits')
    hdr = drizzled_img[0].header

    # Note, shape determined by max shift
    assert(np.shape(drizzled_img[0].data) == (1170, 1170))
    
    mean, med, std = stats.sigma_clipped_stats(drizzled_img[0].data, sigma_upper = 4, sigma_lower = 4, maxiters=5)
    assert(np.isclose(mean, 14.911523, rtol = 1e-2))

    # makes sure there's no nans in the drizzle file
    assert np.sum(np.isnan(drizzled_img[0].data)) == 0

    max_sat = open(combo_dir + 'mag17may21_ob170095_kp.max')
    max_sat = float(max_sat.read())
    assert(np.isclose(max_sat, 11776.6318, rtol = 1))

    # Check first and last stacked headers
    assert hdr['D001DATA'] == 'clean/ob170095_kp/weight/cdwt0016.fits'
    assert hdr['D001DEXP'] == 30.0
    assert hdr['D001OUDA'] == 'combo/mag17may21_ob170095_kp_tmp.fits'
    assert hdr['D001OUWE'] == 'combo/mag17may21_ob170095_kp_sig.fits'
    assert hdr['D001OUCO'] == ''
    assert hdr['D001MASK'] == 'clean/ob170095_kp/masks/mask0016.fits'
    assert hdr['D001WTSC'] == 1.0
    assert hdr['D001KERN'] == 'lanczos3'
    assert hdr['D001PIXF'] == 1.0
    assert hdr['D001COEF'] == ''
    assert hdr['D001XGIM'] == 'clean/ob170095_kp/weight/cdwt0016geo_x.fits'
    assert hdr['D001YGIM'] == 'clean/ob170095_kp/weight/cdwt0016geo_y.fits'
    assert hdr['D001SCAL'] == 1
    assert hdr['D001ROT'] == 0
    assert np.isclose(hdr['D001XSH'], -0.00410188245410836, rtol = 1e-4)
    assert np.isclose(hdr['D001YSH'], 0.001031803046146251, rtol = 1e-4)
    assert hdr['D001SFTU'] == 'pixels'
    assert hdr['D001SFTF'] == 'pixels'
    assert hdr['D001EXKY'] == 'ITIME'
    assert hdr['D001INUN'] == 'counts'
    assert hdr['D001OUUN'] == 'counts'
    assert hdr['D001FVAL'] == '0'
    assert hdr['D013DATA'] == 'clean/ob170095_kp/weight/cdwt0007.fits'
    assert hdr['D013DEXP'] == 30.0
    assert hdr['D013OUDA'] == 'combo/mag17may21_ob170095_kp_tmp.fits'
    assert hdr['D013OUWE'] == 'combo/mag17may21_ob170095_kp_sig.fits'
    assert hdr['D013OUCO'] == ''
    assert hdr['D013MASK'] == 'clean/ob170095_kp/masks/mask0007.fits'
    assert hdr['D013WTSC'] == 1.0
    assert hdr['D013KERN'] == 'lanczos3'
    assert hdr['D013PIXF'] == 1.0
    assert hdr['D013COEF'] == ''
    assert hdr['D013XGIM'] == 'clean/ob170095_kp/weight/cdwt0007geo_x.fits'
    assert hdr['D013YGIM'] == 'clean/ob170095_kp/weight/cdwt0007geo_y.fits'
    assert hdr['D013SCAL'] == 1
    assert hdr['D013ROT'] == 0
    assert np.isclose(hdr['D013XSH'], -9.8446484895112, rtol = 1e-4)
    assert np.isclose(hdr['D013YSH'], 14.925133998969272, rtol = 1e-4)
    assert hdr['D013SFTU'] == 'pixels'
    assert hdr['D013SFTF'] == 'pixels'
    assert hdr['D013EXKY'] == 'ITIME'
    assert hdr['D013INUN'] == 'counts'
    assert hdr['D013OUUN'] == 'counts'
    assert hdr['D013FVAL'] == '0'

    return

@pytest.mark.skip
def test_run_combine_drizzle_with_rotation(self):
    # To be run in the reduce directory of the test_epoch
    mod_path = os.path.dirname(os.path.abspath(data.__file__))
    epoch_dir = mod_path + '/../data/test_epoch/17may21_noiraf/'
    cleanDir = epoch_dir + 'clean/ob170095_kp/'
    combo_dir = epoch_dir + 'combo/'

    _out = epoch_dir + '/combo/mag17may21_ob170095_kp'
    _sub = epoch_dir + '/combo/m17may21_ob170095_kp'
    refImage = epoch_dir + 'clean/ob170095_kp/c0016.fits'
    diffPA = 1
    submaps = 3
    wave = 'kp'
    instrument = instruments.NIRC2()

    roots = ['0007', '0008', '0009', '0010', '0011', '0012', '0013', 
             '0014', '0015', '0016', '0017', '0018', '0019']
    strehls = np.array([0.314, 0.369, 0.348, 0.376, 0.39, 0.373, 
                        0.347, 0.371, 0.357, 0.388, 0.368, 0.326, 0.333])
    fwhm = np.array([65.85, 58.66, 60.65, 57.72, 56.37, 57.86, 
                     57.84, 55.89, 58.17, 55.78, 56.03, 60.95, 59.87])
    weights = np.array([0.06738197424892704, 0.07918454935622317, 0.07467811158798282, 
                        0.08068669527896996, 0.08369098712446352, 0.08004291845493562, 
                        0.07446351931330471,  0.0796137339055794, 0.07660944206008583, 
                        0.08326180257510729, 0.07896995708154506, 0.06995708154506437, 0.07145922746781116])
    

    shiftsTab = Table.read(combo_dir + 'mag17may21_ob170095_kp.shifts', format='ascii', data_start=1)
    
    roots, strehls, fwhm, weights, shiftsTab = data.sort_frames(roots, strehls, fwhm, weights, shiftsTab)
    xysize = 1170

    data.combine_drizzle(xysize, cleanDir, roots, _out, weights, shiftsTab,
                wave, diffPA, fixDAR=True, mask=True, instrument=instrument,
                use_koa_weather=False)

    return

def test_run_submaps():
    # To be run in the reduce directory of the test_epoch
    mod_path = os.path.dirname(os.path.abspath(data.__file__))
    epoch_dir = mod_path + '/../data/test_epoch/17may21/'
    cleanDir = epoch_dir + 'clean/ob170095_kp/'
    combo_dir = epoch_dir + 'combo/'

    _out = epoch_dir + '/combo/mag17may21_ob170095_kp'
    _sub = epoch_dir + '/combo/m17may21_ob170095_kp'
    refImage = epoch_dir + 'clean/ob170095_kp/c0016.fits'
    diffPA = 0
    submaps = 3
    wave = 'kp'
    instrument = instruments.NIRC2()

    roots = ['0007', '0008', '0009', '0010', '0011', '0012', '0013', 
             '0014', '0015', '0016', '0017', '0018', '0019']
    strehls = np.array([0.314, 0.369, 0.348, 0.376, 0.39, 0.373, 
                        0.347, 0.371, 0.357, 0.388, 0.368, 0.326, 0.333])
    fwhm = np.array([65.85, 58.66, 60.65, 57.72, 56.37, 57.86, 
                     57.84, 55.89, 58.17, 55.78, 56.03, 60.95, 59.87])
    weights = np.array([0.06738197424892704, 0.07918454935622317, 0.07467811158798282, 
                        0.08068669527896996, 0.08369098712446352, 0.08004291845493562, 
                        0.07446351931330471,  0.0796137339055794, 0.07660944206008583, 
                        0.08326180257510729, 0.07896995708154506, 0.06995708154506437, 0.07145922746781116])
    

    shiftsTab = Table.read(combo_dir + 'mag17may21_ob170095_kp.shifts', format='ascii', data_start=1)
    
    roots, strehls, fwhm, weights, shiftsTab = data.sort_frames(roots, strehls, fwhm, weights, shiftsTab)
    xysize = 1170

    data.combine_submaps(xysize, cleanDir, roots, _sub, weights,
                    shiftsTab, submaps, wave, diffPA, fixDAR=True,
                    mask=True, instrument=instrument,
                    use_koa_weather=False)

    drizzled_img = fits.open(combo_dir + 'm17may21_ob170095_kp_1.fits')
    hdr = drizzled_img[0].header

    # Note, shape determined by max shift
    assert(np.shape(drizzled_img[0].data) == (1170, 1170))
    
    mean, med, std = stats.sigma_clipped_stats(drizzled_img[0].data, sigma_upper = 4, sigma_lower = 4, maxiters=5)
    assert(np.isclose(mean, 14.930861, rtol = 1e-2))

    # makes sure there's no nans in the drizzle file
    assert np.sum(np.isnan(drizzled_img[0].data)) == 0

    max_sat = open(combo_dir + 'm17may21_ob170095_kp_1.max')
    max_sat = float(max_sat.read())
    assert(np.isclose(max_sat, 11772.3401, rtol = 1))

    # Check first and last stacked headers in the first submap
    assert hdr['D001DATA'] == 'clean/ob170095_kp/weight/cdwt0016.fits'
    assert hdr['D001DEXP'] == 30.0
    assert hdr['D001OUDA'] == 'combo/m17may21_ob170095_kp_1_tmp.fits'
    assert hdr['D001OUWE'] == 'combo/m17may21_ob170095_kp_1_sig.fits'
    assert hdr['D001OUCO'] == ''
    assert hdr['D001MASK'] == 'clean/ob170095_kp/masks/mask0016.fits'
    assert hdr['D001WTSC'] == 1.0
    assert hdr['D001KERN'] == 'lanczos3'
    assert hdr['D001PIXF'] == 1.0
    assert hdr['D001COEF'] == ''
    assert hdr['D001XGIM'] == 'clean/ob170095_kp/weight/cdwt0016geo_x.fits'
    assert hdr['D001YGIM'] == 'clean/ob170095_kp/weight/cdwt0016geo_y.fits'
    assert hdr['D001SCAL'] == 1
    assert hdr['D001ROT'] == 0
    assert np.isclose(hdr['D001XSH'], -0.00410188245410836, rtol = 1e-4)
    assert np.isclose(hdr['D001YSH'], 0.001031803046146251, rtol = 1e-4)
    assert hdr['D001SFTU'] == 'pixels'
    assert hdr['D001SFTF'] == 'pixels'
    assert hdr['D001EXKY'] == 'ITIME'
    assert hdr['D001INUN'] == 'counts'
    assert hdr['D001OUUN'] == 'counts'
    assert hdr['D001FVAL'] == '0'
    assert hdr['D005DATA'] == 'clean/ob170095_kp/weight/cdwt0007.fits'
    assert hdr['D005DEXP'] == 30.0
    assert hdr['D005OUDA'] == 'combo/m17may21_ob170095_kp_1_tmp.fits'
    assert hdr['D005OUWE'] == 'combo/m17may21_ob170095_kp_1_sig.fits'
    assert hdr['D005OUCO'] == ''
    assert hdr['D005MASK'] == 'clean/ob170095_kp/masks/mask0007.fits'
    assert hdr['D005WTSC'] == 1.0
    assert hdr['D005KERN'] == 'lanczos3'
    assert hdr['D005PIXF'] == 1.0
    assert hdr['D005COEF'] == ''
    assert hdr['D005XGIM'] == 'clean/ob170095_kp/weight/cdwt0007geo_x.fits'
    assert hdr['D005YGIM'] == 'clean/ob170095_kp/weight/cdwt0007geo_y.fits'
    assert hdr['D005SCAL'] == 1
    assert hdr['D005ROT'] == 0
    assert np.isclose(hdr['D005XSH'], -9.8446484895112, rtol = 1e-4)
    assert np.isclose(hdr['D005YSH'], 14.925133998969272, rtol = 1e-4)
    assert hdr['D005SFTU'] == 'pixels'
    assert hdr['D005SFTF'] == 'pixels'
    assert hdr['D005EXKY'] == 'ITIME'
    assert hdr['D005INUN'] == 'counts'
    assert hdr['D005OUUN'] == 'counts'
    assert hdr['D005FVAL'] == '0'

    return
