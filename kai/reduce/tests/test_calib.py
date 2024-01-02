import unittest
import os
import numpy as np
from kai.reduce import calib
from kai import instruments
from astropy.io import fits

class TestMakeDarks(unittest.TestCase):
    def test_makedark(self):
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(calib.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/'
        raw_dir = epoch_dir + 'raw/'

        dark_files = range(214, 223 + 1)
        dark_out = 'dark_30.0s_1ca.fits'
        dark_out_full_path = reduce_dir + 'calib/darks/' + dark_out  # only use for testing

        # If file exists, delete it.
        if os.path.isfile(dark_out_full_path):
            os.remove(dark_out_full_path)

        calib.makedark(dark_files, dark_out,
                       raw_dir=raw_dir,
                       reduce_dir=reduce_dir,
                       instrument=nirc2)

        # Check that file was created.
        self.assertTrue(os.path.exists(dark_out_full_path))

        # Check the shape is right.
        dk_img = fits.getdata(dark_out_full_path)
        self.assertEqual(dk_img.shape, (1024, 1024))

        return


class TestMakeFlats(unittest.TestCase):
    def test_makeflat(self):
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(calib.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/'
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
                       reduce_dir=reduce_dir,
                       instrument=nirc2)

        # Check that file was created.
        self.assertTrue(os.path.exists(flat_out_full_path))

        return


if __name__ == '__main__':
    unittest.main()
