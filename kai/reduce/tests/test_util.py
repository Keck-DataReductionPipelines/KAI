import unittest
from kai.reduce import util
import os
from astropy.io import fits
import numpy.testing as nptest

class TestImarith(unittest.TestCase):
    def test_imarith_add_single_file(self):
        mod_path = os.path.dirname(os.path.abspath(util.__file__))
        epoch_dir = os.path.dirname(mod_path + '/../data/test_epoch/17may21/')

        img_file1 = epoch_dir + '/raw/n0216.fits'
        img_file2 = epoch_dir + '/raw/n0217.fits'
        out_file = epoch_dir + '/reduce/calib/tests/test_imarith_add.fits'

        os.makedirs(epoch_dir + '/reduce/calib/tests', exist_ok=True)

        util.imarith(img_file1, '+', img_file2, out_file)

        img1 = fits.getdata(img_file1)
        img2 = fits.getdata(img_file2)
        out = fits.getdata(out_file)

        nptest.assert_almost_equal(img1 + img2, out)

        os.remove(out_file)

        return

    def test_imarith_add_file_list(self):
        mod_path = os.path.dirname(os.path.abspath(util.__file__))
        epoch_dir = os.path.dirname(mod_path + '/../data/test_epoch/17may21/')

        img_files1 = [epoch_dir + '/raw/n0216.fits', epoch_dir + '/raw/n0218.fits']
        img_file2 = epoch_dir + '/raw/n0217.fits'
        out_file = [epoch_dir + '/reduce/calib/tests/test_imarith_add_n0216.fits',
                    epoch_dir + '/reduce/calib/tests/test_imarith_add_n0218.fits']

        os.makedirs(epoch_dir + '/reduce/calib/tests', exist_ok=True)

        util.imarith(img_files1, '+', img_file2, out_file)

        img1 = fits.getdata(img_files1[0])
        img2 = fits.getdata(img_file2)
        out = fits.getdata(out_file[0])

        nptest.assert_almost_equal(img1 + img2, out)

        util.rmall(out_file)

        return

    def test_imarith_add_atfile(self):
        mod_path = os.path.dirname(os.path.abspath(util.__file__))
        epoch_dir = os.path.dirname(mod_path + '/../data/test_epoch/17may21/')
        test_dir = epoch_dir + '/reduce/calib/tests/'
        os.makedirs(test_dir, exist_ok=True)

        img_files1 = [epoch_dir + '/raw/n0216.fits',
                      epoch_dir + '/raw/n0218.fits']
        img_file2 = epoch_dir + '/raw/n0217.fits'
        out_files = [test_dir + 'test_imarith_add_n0216.fits',
                     test_dir + 'test_imarith_add_n0218.fits']

        at_file1_name = test_dir + 'test_imarith_add_atfile1.lis'
        at_file1 = open(at_file1_name, 'w')
        at_file1.write('\n'.join(img_files1))
        at_file1.close()

        at_file3_name = test_dir + 'test_imarith_add_atfile3.lis'
        at_file3 = open(at_file3_name, 'w')
        at_file3.write('\n'.join(out_files))
        at_file3.close()

        util.imarith('@' + at_file1_name,
                     '+', img_file2,
                     '@' + at_file3_name)

        img1 = fits.getdata(img_files1[0])
        img2 = fits.getdata(img_file2)
        out = fits.getdata(out_files[0])

        nptest.assert_almost_equal(img1 + img2, out)

        util.rmall([at_file1_name, at_file3_name] + out_files)

        return

    def test_imarith_subtract_atfile(self):
        mod_path = os.path.dirname(os.path.abspath(util.__file__))
        epoch_dir = os.path.dirname(mod_path + '/../data/test_epoch/17may21/')
        test_dir = epoch_dir + '/reduce/calib/tests/'
        os.makedirs(test_dir, exist_ok=True)

        img_files1 = [epoch_dir + '/raw/n0216.fits',
                      epoch_dir + '/raw/n0218.fits']
        img_file2 = epoch_dir + '/raw/n0217.fits'
        out_files = [test_dir + 'test_imarith_subtract_n0216.fits',
                     test_dir + 'test_imarith_subtract_n0218.fits']

        at_file1_name = test_dir + 'test_imarith_subtract_atfile1.lis'
        at_file1 = open(at_file1_name, 'w')
        at_file1.write('\n'.join(img_files1))
        at_file1.close()

        at_file3_name = test_dir + 'test_imarith_subtract_atfile3.lis'
        at_file3 = open(at_file3_name, 'w')
        at_file3.write('\n'.join(out_files))
        at_file3.close()

        util.imarith('@' + at_file1_name,
                     '-', img_file2,
                     '@' + at_file3_name)

        img1 = fits.getdata(img_files1[0])
        img2 = fits.getdata(img_file2)
        out = fits.getdata(out_files[0])

        nptest.assert_almost_equal(img1 - img2, out)

        util.rmall([at_file1_name, at_file3_name] + out_files)

        return

    def test_imarith_divide_atfile(self):
        mod_path = os.path.dirname(os.path.abspath(util.__file__))
        epoch_dir = os.path.dirname(mod_path + '/../data/test_epoch/17may21/')
        test_dir = epoch_dir + '/reduce/calib/tests/'
        os.makedirs(test_dir, exist_ok=True)

        img_files1 = [epoch_dir + '/raw/n0216.fits',
                      epoch_dir + '/raw/n0218.fits']
        img_file2 = epoch_dir + '/raw/n0217.fits'
        out_files = [test_dir + 'test_imarith_divide_n0216.fits',
                     test_dir + 'test_imarith_divide_n0218.fits']

        at_file1_name = test_dir + 'test_imarith_divide_atfile1.lis'
        at_file1 = open(at_file1_name, 'w')
        at_file1.write('\n'.join(img_files1))
        at_file1.close()

        at_file3_name = test_dir + 'test_imarith_divide_atfile3.lis'
        at_file3 = open(at_file3_name, 'w')
        at_file3.write('\n'.join(out_files))
        at_file3.close()

        util.imarith('@' + at_file1_name,
                     '/', img_file2,
                     '@' + at_file3_name)

        img1 = fits.getdata(img_files1[0])
        img2 = fits.getdata(img_file2)
        out = fits.getdata(out_files[0])

        nptest.assert_almost_equal(img1 / img2, out)

        util.rmall([at_file1_name, at_file3_name] + out_files)

        return

    def test_imarith_multiply_atfile(self):
        mod_path = os.path.dirname(os.path.abspath(util.__file__))
        epoch_dir = os.path.dirname(mod_path + '/../data/test_epoch/17may21/')
        test_dir = epoch_dir + '/reduce/calib/tests/'
        os.makedirs(test_dir, exist_ok=True)

        img_files1 = [epoch_dir + '/raw/n0216.fits',
                      epoch_dir + '/raw/n0218.fits']
        img_file2 = epoch_dir + '/raw/n0217.fits'
        out_files = [test_dir + 'test_imarith_multiply_n0216.fits',
                     test_dir + 'test_imarith_multiply_n0218.fits']

        at_file1_name = test_dir + 'test_imarith_multiply_atfile1.lis'
        at_file1 = open(at_file1_name, 'w')
        at_file1.write('\n'.join(img_files1))
        at_file1.close()

        at_file3_name = test_dir + 'test_imarith_multiply_atfile3.lis'
        at_file3 = open(at_file3_name, 'w')
        at_file3.write('\n'.join(out_files))
        at_file3.close()

        util.imarith('@' + at_file1_name,
                     '*', img_file2,
                     '@' + at_file3_name)

        img1 = fits.getdata(img_files1[0])
        img2 = fits.getdata(img_file2)
        out = fits.getdata(out_files[0])

        nptest.assert_almost_equal(img1 * img2, out)

        util.rmall([at_file1_name, at_file3_name] + out_files)

        return

    def test_imarith_atfile_mismatch_length(self):
        mod_path = os.path.dirname(os.path.abspath(util.__file__))
        epoch_dir = os.path.dirname(mod_path + '/../data/test_epoch/17may21/')
        test_dir = epoch_dir + '/reduce/calib/tests/'
        os.makedirs(test_dir, exist_ok=True)

        img_files1 = [epoch_dir + '/raw/n0216.fits',
                      epoch_dir + '/raw/n0218.fits']
        img_file2 = epoch_dir + '/raw/n0217.fits'
        out_files = test_dir + 'test_imarith_add_n0216.fits'

        at_file1_name = test_dir + 'test_imarith_add_atfile1.lis'
        at_file1 = open(at_file1_name, 'w')
        at_file1.write('\n'.join(img_files1))
        at_file1.close()

        with self.assertRaises(RuntimeError):
            util.imarith('@' + test_dir + 'test_imarith_add_atfile1.lis',
                         '*', img_file2,
                         out_files)

        util.rmall([at_file1_name])

        return


if __name__ == '__main__':
    unittest.main()
