import unittest
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

class TestClean(unittest.TestCase):
     def test_clean_iraf():#self):
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/'
        raw_dir = epoch_dir + 'raw/'

        dark_files = range(214, 223 + 1)
        dark_out = 'dark_30.0s_1ca.fits'
        dark_out_full_path = reduce_dir + 'calib/darks/' + dark_out  # only use for testing

        sci_files = ['{0:04d}_iraf'.format(ii) for ii in range(7, 19+1)]
        refSrc = [556, 564]
        data.clean_iraf(sci_files, 'ob170095', 'kp', refSrc, refSrc, field='ob170095',
                   raw_dir='../raw/', clean_dir='../clean/', instrument=nirc2)

        return
         
     def test_clean(self):
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/'
        raw_dir = epoch_dir + 'raw/'

        dark_files = range(214, 223 + 1)
        dark_out = 'dark_30.0s_1ca.fits'
        dark_out_full_path = reduce_dir + 'calib/darks/' + dark_out  # only use for testing

        sci_files = ['{0:03d}'.format(ii) for ii in range(7, 19+1)]
        refSrc = [556, 564]
        data.clean(sci_files, 'ob170095', 'kp', refSrc, refSrc, field='ob170095',
                   raw_dir='../raw/', clean_dir='../clean/', instrument=nirc2)

        return


class TestCombine(unittest.TestCase):
    def test_combine_iraf():
        nirc2 = instruments.NIRC2()

        sci_files = ['{0:03d}_iraf'.format(ii) for ii in range(7, 19+1)]
        
        data.combine_iraf(sci_files, 'kp', '17may21', trim=1,
                 weight='strehl', field='ob170095', submaps=3)

        return

    def test_combine_noiraf(self):
        nirc2 = instruments.NIRC2()

        sci_files = ['{0:03d}'.format(ii) for ii in range(7, 19+1)]
        
        data.combine(sci_files, 'kp', '17may21', trim=1,
                 weight='strehl', field='ob170095', submaps=3)

        return


class TestRotImg(unittest.TestCase):
    def test_run_rot_img_noiraf(self):
        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        clean_dir = epoch_dir + 'clean/ob170095_kp/'

        img_root = '0007'
        phi = 30.0   # arbitrary.

        rot_img_filename = clean_dir + 'r' + img_root + '.fits'

        data.rot_img(img_root, phi, cleanDir=clean_dir)
        shutil.move(rot_img_filename, rot_img_filename.replace('.fits', '_noiraf.fits'))

        return

    def test_run_rot_img_iraf(self):
        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        clean_dir = epoch_dir + 'clean/ob170095_kp/'

        img_root = '0007'
        phi = 30.0   # arbitrary.

        rot_img_filename = clean_dir + 'r' + img_root + '.fits'

        data.rot_img_iraf(img_root, phi, cleanDir=clean_dir)
        shutil.move(rot_img_filename, rot_img_filename.replace('.fits', '_iraf.fits'))

        return

    def test_compare_rot_img_iraf_noiraf(self):
        """Compare IRAF vs. no IRAF image rotation routine."""

        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        clean_dir = epoch_dir + 'clean/ob170095_kp/'

        img_root = '0007'
        phi = 30.0   # arbitrary.

        rot_img_filename = clean_dir + 'r' + img_root + '.fits'

        img_rot_iraf, img_rot_iraf_hdr = fits.getdata(rot_img_filename.replace('.fits', '_iraf.fits'), header = True)
        img_rot_noiraf, img_rot_noiraf_hdr = fits.getdata(rot_img_filename.replace('.fits', '_noiraf.fits'), header = True)

        assert img_rot_iraf.shape == img_rot_noiraf.shape

        # Checks headers are the same
        hdr_keys = list(dict(img_rot_iraf_hdr).keys())
        skip_keys = ['EXTEND', 'DATE', 'IRAF-TLM', 'COMMENT', '']
        for key in hdr_keys:
            if key in skip_keys:
                continue
            if type(img_rot_iraf_hdr[key]) == str:
                assert img_rot_iraf_hdr[key] == img_rot_noiraf_hdr[key]
            # Reference pixel
            elif key == 'CRPIX1' or key == 'CRPIX2':
                assert math.isclose(img_rot_iraf_hdr[key], img_rot_noiraf_hdr[key], abs_tol = 0.1)
            # CCD to image coords - some sort of shift
            elif key == 'LTV1' or key == 'LTV2':
                assert math.isclose(img_rot_iraf_hdr[key], img_rot_noiraf_hdr[key], abs_tol = 0.5)
            else:
                assert math.isclose(img_rot_iraf_hdr[key], img_rot_noiraf_hdr[key], abs_tol = 0.01)
        
        return

class TestCosmicRayRemoval(unittest.TestCase):
    def test_run_clean_cosmicrays_iraf(self):
        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/kp/sci_ob170095/'

        img_root = '0007'
    
        to_clean_img_filename = reduce_dir + 'ff' + img_root + '.fits'
        clean_img_filename = reduce_dir + 'ff' + img_root + '_f.fits'
        
        crmask_filename = reduce_dir + 'crmask' + img_root + '_iraf.fits'
        _statmask = reduce_dir + 'stat_mask' + img_root + '.fits'
        _ff_s = reduce_dir + 'ff' + img_root + '_s.fits'

        util.rmall([crmask_filename])

        ### Fix bad pixels ###
        # Produces _ff_f file
        bfixpix.bfixpix(to_clean_img_filename, _statmask)
        util.rmall([_ff_s])

        data.clean_cosmicrays_iraf(_ff=clean_img_filename, _mask=crmask_filename, wave='kp')
        shutil.move(clean_img_filename, clean_img_filename.replace('.fits', '_iraf.fits'))

        return

    def test_run_clean_cosmicrays_noiraf(self):
        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/kp/sci_ob170095/'

        #for ii in range(7, 20):
        img_root = '0007' #'{0:04d}'.format(ii)
    
        to_clean_img_filename = reduce_dir + 'ff' + img_root + '.fits'
        clean_img_filename = reduce_dir + 'ff' + img_root + '_f.fits'
        
        crmask_filename = reduce_dir + 'crmask' + img_root + '_noiraf.fits'
        _statmask = reduce_dir + 'stat_mask' + img_root + '.fits'
        _ff_s = reduce_dir + 'ff' + img_root + '_s.fits'

        util.rmall([crmask_filename])

        ### Fix bad pixels ###
        # Produces _ff_f file
        bfixpix.bfixpix(to_clean_img_filename, _statmask)
        util.rmall([_ff_s])

        data.clean_cosmicrays(_ff=clean_img_filename, _mask=crmask_filename, wave='kp', thresh=5, mbox=5, rbox=10, fratio = 0.4, gbox=0)
        shutil.move(clean_img_filename, clean_img_filename.replace('.fits', '_noiraf.fits'))

        return

    def test_compare_clean_cosmicrays_iraf_noiraf(self):
        """Compare IRAF vs. no IRAF cosmicray removal routine."""

        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        """
        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/kp/sci_ob170095/'
        
        img_root = '0007'
        
        iraf = fits.getdata(reduce_dir + 'crmask{}_iraf.fits'.format(img_root))
        noiraf = fits.getdata(reduce_dir + 'crmask{}_noiraf.fits'.format(img_root))

        noiraf_img = fits.getdata(reduce_dir + 'ff{}_f_noiraf.fits'.format(img_root))
        iraf_img = fits.getdata(reduce_dir + 'ff{}_f_iraf.fits'.format(img_root))
        """

        iraf_epoch_dir = mod_path + '/../data/test_epoch/17may21_iraf/'
        iraf_reduce_dir = iraf_epoch_dir + 'reduce/kp/sci_ob170095/'

        noiraf_epoch_dir = mod_path + '/../data/test_epoch/17may21_noiraf/'
        noiraf_reduce_dir = noiraf_epoch_dir + 'reduce/kp/sci_ob170095/'
        
        img_root = '0007'
        
        iraf = fits.getdata(iraf_reduce_dir + 'crmask{}.fits'.format(img_root))
        noiraf = fits.getdata(noiraf_reduce_dir + 'crmask{}.fits'.format(img_root))

        noiraf_img = fits.getdata(noiraf_reduce_dir + 'ff{}_f.fits'.format(img_root))
        iraf_img = fits.getdata(iraf_reduce_dir + 'ff{}_f.fits'.format(img_root))

        total_cosmic_rays_iraf = np.sum(iraf)

        in_iraf_only = len(np.where((iraf - noiraf) == 1)[0])
        in_noiraf_only = len(np.where((iraf - noiraf) == -1)[0])
        total_diff = in_iraf_only + in_noiraf_only

        # Total difference in pixels 
        # between the two is less than 35% the total 
        assert total_diff/total_cosmic_rays_iraf < 0.35
        print(total_diff/total_cosmic_rays_iraf)

        # Difference in each individual one less than 25% of total
        # i.e. make sure iraf didn't find way more than no-iraf, just different
        assert in_iraf_only/total_cosmic_rays_iraf < 0.25
        assert in_noiraf_only/total_cosmic_rays_iraf < 0.25

        # Check replaced values are within 40 counts
        # and that there are fewer than 10 pixels with
        # a difference greater than 10 counts
        difference = np.abs(iraf_img - noiraf_img)
        assert np.max(difference) < 55
        assert len(np.where(difference > 10)[0]) < 15
        
        return

class TestCombineRegister(unittest.TestCase):
    def test_run_combine_register_iraf(self):
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
        
        _out = comboDir + 'mag' + outroot + '_' + wave + '_iraf'
        
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
        shiftsTab = data.combine_register_iraf(_out, refImage, diffPA)

        return

    def test_run_combine_register_noiraf(self): #self
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
        
        _out = comboDir + 'mag' + outroot + '_' + wave + '_noiraf'
        
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

        return

    def test_compare_combine_register_iraf_noiraf(self):
        """Compare IRAF vs. no IRAF combine register routine."""

        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        
        combo_dir = mod_path + '/../data/test_epoch/17may21/combo/'
        
        iraf = Table.read(combo_dir + 'mag17may21_ob170095_kp_iraf.shifts', format='ascii', data_start=1)
        noiraf = Table.read(combo_dir + 'mag17may21_ob170095_kp_noiraf.shifts', format='ascii', data_start=1)

        shift_differences_y = []
        shift_differences_x = []
        shift_differences_total = []
        for i in range(len(iraf)):
            shift_differences_y.append(noiraf['col1'][i] - iraf['col1'][i])
            shift_differences_x.append(noiraf['col2'][i] - iraf['col2'][i])
            shift_differences_total.append(np.sqrt((noiraf['col1'][i] - iraf['col1'][i])**2 + (noiraf['col2'][i] - iraf['col2'][i])**2))

        # Mean difference is less than 0.1 pixels
        assert np.mean(shift_differences_total) < 0.1

        # Max difference is less than 0.3 pixels
        assert np.max(shift_differences_total) < 0.3
        
        return

class TestCleanDrizzle(unittest.TestCase):
    """
    Tests create the drizzled files and weighting files and does basic
    tests on no iraf version.
    Note that we do not directly compare the iraf and no iraf cases.
    The weighting file looks much more sharp/defined for the no iraf case
    and more smooth for the iraf case. It's unclear which is more correct.
    It causes differences in the nosie background but it shouldn't change
    the star positions.
    """
    def test_run_clean_drizzle_iraf(self):
        # To be run in the reduce directory of the test_epoch
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        rawDir = epoch_dir + 'raw/'
        
        instrument = instruments.NIRC2()
        files = ['{0:04d}_iraf'.format(ii) for ii in range(7, 19+1)]
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
        data.setup_drizzle_iraf(imgsize)
        
        os.chdir('kp/sci_ob170095/')
        for f in files:
            _bp = instrument.make_filenames([f], prefix='bp')[0]
            _ce = instrument.make_filenames([f], prefix='ce')[0]
            _wgt = instrument.make_filenames([f], prefix='wgt')[0]
            _ff = instrument.make_filenames([f], prefix='ff')[0]
            _ff_f = _ff.replace('_iraf.fits', '_f_iraf.fits')
            _dlog_tmp = instrument.make_filenames([f], prefix='driz')[0]
            _dlog = _dlog_tmp.replace('.fits', '.log')

            ### Background Subtraction ###
            if os.path.exists(_bp) == False:
                bkg = data.clean_bkgsubtract(_ff_f, _bp)
            
            ### Drizzle individual file ###
            data.clean_drizzle_iraf(distXgeoim, distYgeoim, _bp, _ce, _wgt, _dlog,
                          fixDAR=True, instrument=instrument,
                          use_koa_weather=False)
        os.chdir('../..')

        return

    def test_run_clean_drizzle_noiraf(self):
        # To be run in the reduce directory of the test_epoch
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        rawDir = epoch_dir + 'raw/'
        
        instrument = instruments.NIRC2()
        files = ['{0:04d}_noiraf'.format(ii) for ii in range(7, 19+1)]
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
        
        os.chdir('kp/sci_ob170095/')
        for f in [files[0]]:
            _bp = instrument.make_filenames([f], prefix='bp')[0]
            _ce = instrument.make_filenames([f], prefix='ce')[0]
            _wgt = instrument.make_filenames([f], prefix='wgt')[0]
            _ff = instrument.make_filenames([f], prefix='ff')[0]
            _ff_f = _ff.replace('_noiraf.fits', '_f_noiraf.fits')
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

        drizzled_file = fits.open('ce0007_noiraf.fits')
        weight_file = fits.open('wgt0007_noiraf.fits')

        # check shape is correct
        assert np.array_equal(np.shape(drizzled_file[0].data), np.array([1024, 1024]))
        assert np.array_equal(np.shape(weight_file[0].data), np.array([1024, 1024]))

        # makes sure there's no nans in the drizzle file
        assert np.sum(np.isnan(drizzled_file[0].data)) == 0
        
        os.chdir('../..')

        return

class TestCombineDrizzle(unittest.TestCase):
    """
    Tests create the drizzled files and weighting files and compares output
    data file and saturation file.
    
    Note that we do not directly compare the iraf and no iraf cases weighting files.
    The weighting file looks much more sharp/defined for the no iraf case
    and more smooth for the iraf case. It's unclear which is more correct.
    It causes differences in the nosie background but it shouldn't change
    the star positions.
    """
    def test_run_combine_drizzle_iraf():#self):
        # To be run in the reduce directory of the test_epoch
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21_iraf/'
        cleanDir = epoch_dir + 'clean/ob170095_kp/'
        combo_dir = epoch_dir + 'combo/'

        _out = epoch_dir + '/combo/mag17may21_ob170095_kp'
        _sub = epoch_dir + '/combo/m17may21_ob170095_kp'
        refImage = epoch_dir + 'clean/ob170095_kp/c0016.fits'
        diffPA = 0
        submaps = 3
        wave = 'kp'
        instrument = instruments.NIRC2()

        # After trimming
        roots = ['0007', '0008', '0009', '0010', '0011', '0012', 
                 '0013', '0014', '0015', '0016', '0017', '0018', '0019']
        strehls = np.array([0.314, 0.369, 0.348, 0.375, 0.388, 0.373, 
                            0.347, 0.371, 0.356, 0.387, 0.367, 0.325, 0.333])
        fwhm = np.array([65.87, 58.67, 60.66, 57.74, 56.38, 57.87, 57.89, 
                         55.93, 58.22, 55.82, 56.07, 60.98,  59.9])
        weights = np.array([0.0674833440791, 0.0793036750484, 0.0747904577692, 
                            0.0805931656995, 0.0833870621105, 0.0801633354825, 
                            0.0745755426606, 0.0797335052654, 0.0765097786374, 
                            0.0831721470019, 0.0788738448313, 0.0698474102729, 0.0715667311412])

        shiftsTab = Table.read(combo_dir + 'mag17may21_ob170095_kp.shifts', format='ascii', data_start=1)
        roots, strehls, fwhm, weights, shiftsTab = data.sort_frames(roots, strehls, fwhm, weights, shiftsTab)
        xysize = 1166

        data.combine_drizzle_iraf(xysize, cleanDir, roots, _out, weights, shiftsTab,
                    wave, diffPA, fixDAR=True, mask=True, instrument=instrument,
                    use_koa_weather=False)

        return

    def test_run_combine_drizzle_noiraf(self):
        # To be run in the reduce directory of the test_epoch
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21_noiraf/'
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

        return

    def test_compare_combine_drizzle_iraf_noiraf(self):
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        data_dir_iraf = mod_path + '/../data/test_epoch/17may21_iraf/combo/'
        data_dir_noiraf = mod_path + '/../data/test_epoch/17may21_noiraf/combo/'
        
        iraf_drizzled_img = fits.open(data_dir_iraf + 'mag17may21_ob170095_kp.fits')
        noiraf_drizzled_img = fits.open(data_dir_noiraf + 'mag17may21_ob170095_kp.fits')

        #no iraf puts it in the corner of the inflated images 
        #whereas iraf puts it in the center, so we shift them to compare them
        y = 71
        x = 71
        iraf_drizzled_img_cut = iraf_drizzled_img[0].data[y:y+1024, x:x+1024]
        noiraf_drizzled_img_cut = noiraf_drizzled_img[0].data[0:1024, 0:1024]

        # Note total size of the images are different because the shifts
        # are slightly different which is used to calculate the padding.
        
        zero_idxs = np.logical_and((noiraf_drizzled_img_cut == 0), (iraf_drizzled_img_cut == 0))
        assert(np.isclose(np.mean((iraf_drizzled_img_cut/noiraf_drizzled_img_cut)[~zero_idxs]), 1, rtol = 1e-2))

        # makes sure there's no nans in the drizzle file
        assert np.sum(np.isnan(iraf_drizzled_img[0].data)) == 0
        assert np.sum(np.isnan(noiraf_drizzled_img[0].data)) == 0

        # check max saturation levels are close
        max_sat_noiraf = open(data_dir_noiraf + 'mag17may21_ob170095_kp.max')
        max_sat_iraf = open(data_dir_iraf + 'mag17may21_ob170095_kp.max')
        max_sat_noiraf = float(max_sat_noiraf.read())
        max_sat_iraf = float(max_sat_iraf.read())
        assert(np.isclose(max_sat_noiraf, max_sat_iraf, rtol = 1))

        return

    def test_run_combine_drizzle_with_rotation_iraf():#self):
        # To be run in the reduce directory of the test_epoch
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21_iraf/'
        cleanDir = epoch_dir + 'clean/ob170095_kp/'
        combo_dir = epoch_dir + 'combo/'

        _out = epoch_dir + '/combo/mag17may21_ob170095_kp'
        _sub = epoch_dir + '/combo/m17may21_ob170095_kp'
        refImage = epoch_dir + 'clean/ob170095_kp/c0016.fits'
        diffPA = 1
        submaps = 3
        wave = 'kp'
        instrument = instruments.NIRC2()

        # After trimming
        roots = ['0007', '0008', '0009', '0010', '0011', '0012', 
                 '0013', '0014', '0015', '0016', '0017', '0018', '0019']
        strehls = np.array([0.314, 0.369, 0.348, 0.375, 0.388, 0.373, 
                            0.347, 0.371, 0.356, 0.387, 0.367, 0.325, 0.333])
        fwhm = np.array([65.87, 58.67, 60.66, 57.74, 56.38, 57.87, 57.89, 
                         55.93, 58.22, 55.82, 56.07, 60.98,  59.9])
        weights = np.array([0.0674833440791, 0.0793036750484, 0.0747904577692, 
                            0.0805931656995, 0.0833870621105, 0.0801633354825, 
                            0.0745755426606, 0.0797335052654, 0.0765097786374, 
                            0.0831721470019, 0.0788738448313, 0.0698474102729, 0.0715667311412])

        shiftsTab = Table.read(combo_dir + 'mag17may21_ob170095_kp.shifts', format='ascii', data_start=1)
        roots, strehls, fwhm, weights, shiftsTab = data.sort_frames(roots, strehls, fwhm, weights, shiftsTab)
        xysize = 1166

        data.combine_drizzle_iraf(xysize, cleanDir, roots, _out, weights, shiftsTab,
                    wave, diffPA, fixDAR=True, mask=True, instrument=instrument,
                    use_koa_weather=False)

        return

    def test_run_combine_drizzle_with_rotation_noiraf(self):
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

class TestCombineSubmaps(unittest.TestCase):
    """
    Tests create the drizzled files and weighting files and does basic
    tests on no iraf version.
    Note that we do not directly compare the iraf and no iraf cases.
    The weighting file looks much more sharp/defined for the no iraf case
    and more smooth for the iraf case. It's unclear which is more correct.
    It causes differences in the nosie background but it shouldn't change
    the star positions.
    """
    def test_run_submaps_drizzle_iraf():#self):
        # To be run in the reduce directory of the test_epoch
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21_iraf/'
        cleanDir = epoch_dir + 'clean/ob170095_kp/'
        combo_dir = epoch_dir + 'combo/'

        _out = epoch_dir + '/combo/mag17may21_ob170095_kp'
        _sub = epoch_dir + '/combo/m17may21_ob170095_kp'
        refImage = epoch_dir + 'clean/ob170095_kp/c0016.fits'
        diffPA = 0
        submaps = 3
        wave = 'kp'
        instrument = instruments.NIRC2()

        # After trimming
        roots = ['0007', '0008', '0009', '0010', '0011', '0012', 
                 '0013', '0014', '0015', '0016', '0017', '0018', '0019']
        strehls = np.array([0.314, 0.369, 0.348, 0.375, 0.388, 0.373, 
                            0.347, 0.371, 0.356, 0.387, 0.367, 0.325, 0.333])
        fwhm = np.array([65.87, 58.67, 60.66, 57.74, 56.38, 57.87, 57.89, 
                         55.93, 58.22, 55.82, 56.07, 60.98,  59.9])
        weights = np.array([0.0674833440791, 0.0793036750484, 0.0747904577692, 
                            0.0805931656995, 0.0833870621105, 0.0801633354825, 
                            0.0745755426606, 0.0797335052654, 0.0765097786374, 
                            0.0831721470019, 0.0788738448313, 0.0698474102729, 0.0715667311412])

        shiftsTab = Table.read(combo_dir + 'mag17may21_ob170095_kp.shifts', format='ascii', data_start=1)
        roots, strehls, fwhm, weights, shiftsTab = data.sort_frames(roots, strehls, fwhm, weights, shiftsTab)
        xysize = 1166

        data.combine_submaps_iraf(xysize, cleanDir, roots, _sub, weights,
                        shiftsTab, submaps, wave, diffPA, fixDAR=True,
                        mask=True, instrument=instrument,
                        use_koa_weather=False)

        return

    def test_run_submaps_noiraf(self):
        # To be run in the reduce directory of the test_epoch
        mod_path = os.path.dirname(os.path.abspath(data.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21_noiraf/'
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

        return


class TestMosaicRegister(unittest.TestCase):
    def test_run_mosaic_register_iraf(self):
        mod_path = os.path.dirname(os.path.abspath(data.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/kp/sci_ob170095/'

        img_root = '0007'
    
        to_clean_img_filename = reduce_dir + 'ff' + img_root + '.fits'
        clean_img_filename = reduce_dir + 'ff' + img_root + '_f.fits'
        
        crmask_filename = reduce_dir + 'crmask' + img_root + '_iraf.fits'
        _statmask = reduce_dir + 'stat_mask' + img_root + '.fits'
        _ff_s = reduce_dir + 'ff' + img_root + '_s.fits'

        util.rmall([crmask_filename])

        ### Fix bad pixels ###
        # Produces _ff_f file
        bfixpix.bfixpix(to_clean_img_filename, _statmask)
        util.rmall([_ff_s])

        data.clean_cosmicrays_iraf(_ff=clean_img_filename, _mask=crmask_filename, wave='kp')
        shutil.move(clean_img_filename, clean_img_filename.replace('.fits', '_iraf.fits'))

        return

if __name__ == '__main__':
    unittest.main()
