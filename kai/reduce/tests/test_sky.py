import unittest
import os
from kai.reduce import sky
from kai import instruments
from astropy.io import fits


class TestMakeSky(unittest.TestCase):
    def test_makesky_nodark(self):
        nirc2 = instruments.NIRC2()

        mod_path = os.path.dirname(os.path.abspath(sky.__file__))

        epoch_dir = mod_path + '/../data/test_epoch/17may21/'
        reduce_dir = epoch_dir + 'reduce/'
        raw_dir = epoch_dir + 'raw/'

        sky_files = range(143, 151 + 1)
        sky_out = 'sky_kp.fits'
        sky_out_full_path = reduce_dir + 'kp/sky_ob170095/' + sky_out

        # If file exists, delete it.
        if os.path.isfile(sky_out_full_path):
            os.remove(sky_out_full_path)

        sky.makesky(sky_files, 'ob170095', 'kp',
                    dark_frame=None, skyscale=True,
                    raw_dir=raw_dir,
                    reduce_dir=reduce_dir,
                    instrument=nirc2)

        # Check that file was created.
        self.assertTrue(os.path.exists(sky_out_full_path))

        # Check the shape is right.
        sk_img = fits.getdata(sky_out_full_path)
        self.assertEqual(sk_img.shape, (1024, 1024))

        return


if __name__ == '__main__':
    unittest.main()
