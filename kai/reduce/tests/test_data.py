import unittest
import shutil
import os
from kai import instruments
from kai.reduce import data
from astropy.io import fits

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

        img_rot_iraf = fits.getdata(rot_img_filename.replace('.fits', '_iraf.fits'))
        img_rot_noiraf = fits.getdata(rot_img_filename.replace('.fits', '_noiraf.fits'))

        assert img_rot_iraf.shape == img_rot_noiraf.shape

        # Add WCS keyword test. Make sure the input and output CD are different and the
        # IRAF CD and NO_IRAF CD are close to the same. 

        return


if __name__ == '__main__':
    unittest.main()
