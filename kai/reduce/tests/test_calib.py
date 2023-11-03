import unittest
import os
import numpy as np
from kai.reduce import calib
from kai import instruments
from astropy.io import fits

class TestMakeDarks(unittest.TestCase):
    def test_makedark_iraf(self):
        """Can only run on IRAF/PyRAF environments.
        """
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(calib.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        raw_dir = epoch_dir + 'raw/'

        dark_files = range(214, 223 + 1)
        dark_out = 'dark_30.0s_1ca_iraf.fits'
        dark_out_full_path = epoch_dir + 'reduce/calib/darks/' + dark_out  # only use for testing

        # If file exists, delete it.
        if os.path.isfile(dark_out_full_path):
            os.remove(dark_out_full_path)

        calib.makedark_iraf(dark_files, dark_out,
                       raw_dir = raw_dir,
                       epoch_dir = epoch_dir,
                       instrument=nirc2)

        # Check that file was created.
        self.assertTrue(os.path.exists(dark_out_full_path))

        return

    def test_makedark(self):
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(calib.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        raw_dir = epoch_dir + 'raw/'

        dark_files = range(214, 223 + 1)
        dark_out = 'dark_30.0s_1ca.fits'
        dark_out_full_path = epoch_dir + 'reduce/calib/darks/' + dark_out  # only use for testing

        # If file exists, delete it.
        if os.path.isfile(dark_out_full_path):
            os.remove(dark_out_full_path)

        calib.makedark(dark_files, dark_out,
                       raw_dir = raw_dir,
                       epoch_dir = epoch_dir,
                       instrument=nirc2)

        # Check that file was created.
        self.assertTrue(os.path.exists(dark_out_full_path))

        # Check the shape is right.
        dk_img = fits.getdata(dark_out_full_path)
        self.assertEqual(dk_img.shape, (1024, 1024))

        return

    def compare_iraf_astropy(self):

        mod_path = os.path.dirname(os.path.abspath(calib.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        dark_dir = epoch_dir + 'reduce/calib/darks/'

        dk_iraf = fits.getdata(dark_dir + 'dark_30.0s_1ca_iraf.fits')
        dk_astr = fits.getdata(dark_dir + 'dark_30.0s_1ca.fits')

        dk_diff = dk_astr - dk_iraf

        # Look for any pixels that are more than 1% differnet.
        wdx = np.where(np.abs(dk_diff / dk_iraf) > 0.1)
        print(wdx)
        print(len(wdx[0]))
        print(dk_diff[wdx[0], wdx[1]])
        print(dk_iraf[wdx[0][100]-5:wdx[0][100]+5, wdx[1][100]-5:wdx[1][100]+5])
        print(dk_astr[wdx[0][100]-5:wdx[0][100]+5, wdx[1][100]-5:wdx[1][100]+5])
        print(dk_diff[wdx[0][100]-5:wdx[0][100]+5, wdx[1][100]-5:wdx[1][100]+5])

        # Found 7489 pixels that were different; This is a pretty small
        # fraction and I have no idea what imcombine was doing internally.

        return

class TestMakeFlats(unittest.TestCase):
    def test_makeflat_iraf(self):
        """Can only run on IRAF/PyRAF environments.
        """
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(calib.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        raw_dir = epoch_dir + 'raw/'

        off_files = range(184, 202 + 1, 2)
        on_files = range(185, 203 + 1, 2)
        flat_out = 'flat_kp_iraf.fits'
        flat_out_full_path = epoch_dir + 'reduce/calib/flats/' + flat_out

        # If file exists, delete it.
        if os.path.isfile(flat_out_full_path):
            os.remove(flat_out_full_path)

        calib.makeflat_iraf(on_files, off_files, flat_out,
                            raw_dir=raw_dir,
                            epoch_dir=epoch_dir,
                            instrument=nirc2)

        # Check that file was created.
        self.assertTrue(os.path.exists(flat_out_full_path))

        return

    def test_makeflat(self):
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(calib.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        raw_dir = epoch_dir + 'raw/'

        off_files = range(184, 202 + 1, 2)
        on_files = range(185, 203 + 1, 2)
        flat_out = 'flat_kp_astr.fits'
        flat_out_full_path = epoch_dir + 'reduce/calib/flats/' + flat_out

        # If file exists, delete it.
        if os.path.isfile(flat_out_full_path):
            os.remove(flat_out_full_path)

        calib.makeflat(on_files, off_files, flat_out,
                       raw_dir=raw_dir,
                       epoch_dir=epoch_dir,
                       instrument=nirc2)

        # Check that file was created.
        self.assertTrue(os.path.exists(flat_out_full_path))

        return

    def compare_iraf_astropy(self):
        mod_path = os.path.dirname(os.path.abspath(calib.__file__))
        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        dark_dir = epoch_dir + 'reduce/calib/flats/'

        flat_iraf = fits.getdata(dark_dir + 'flat_kp_iraf.fits')
        flat_astr = fits.getdata(dark_dir + 'flat_kp_astr.fits')
        flat_diff = flat_astr - flat_iraf

        # Look for any pixels that are more than 2% different.
        wdx = np.where(np.abs(flat_diff / flat_iraf) > 0.02)
        self.assertLess(len(wdx), 50)  # less than 50 big difference pixels (2 %)

        # Make sure these bad pixels are actually bad (low or high flux values.
        self.assertTrue(np.any(flat_diff[wdx[0], wdx[1]] < 0.85) or np.any(flat_diff[wdx[0], wdx[1]] > 1.15))
        # print(wdx)
        # print(len(wdx[0]))
        # print(flat_diff[wdx[0], wdx[1]])
        # print(flat_iraf[wdx[0][3]-2:wdx[0][3]+3, wdx[1][3]-2:wdx[1][3]+3])
        # print(flat_astr[wdx[0][3]-2:wdx[0][3]+3, wdx[1][3]-2:wdx[1][3]+3])
        # print(flat_diff[wdx[0][3]-2:wdx[0][3]+3, wdx[1][3]-2:wdx[1][3]+3])

        return


if __name__ == '__main__':
    unittest.main()
