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
        
        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/kp/sci_ob170095/'
        
        img_root = '0007'
        
        iraf = fits.getdata(reduce_dir + 'crmask{}_iraf.fits'.format(img_root))
        noiraf = fits.getdata(reduce_dir + 'crmask{}_noiraf.fits'.format(img_root))

        noiraf_img = fits.getdata(reduce_dir + 'ff{}_f_noiraf.fits'.format(img_root))
        iraf_img = fits.getdata(reduce_dir + 'ff{}_f_iraf.fits'.format(img_root))

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
        data.setup_drizzle(imgsize)
        
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

        # makes sure there's no nans in the dirrzle file
        assert np.sum(np.isnan(drizzled_file[0].data)) == 0
        
        os.chdir('../..')

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
