from kai.reduce import calibrate
from astropy.table import Table
import os.path

def test_star_names():
    mod_path = os.path.abspath(calibrate.__file__)

    photo_calib_file = mod_path + '../tests/test_calibrate/photo_calib.dat'
    starlist = mod_path + '../tests/test_calibrate/mag95jun_0.8d_stf.lis'
    
    args = '-f 1 -R -N ' + photo_calib_file + ' -M 7 -T 0.0 -V -S S3-24,S3-37,33E,S3-169,S1-54,S2-47,S2-2,S1-3,S3-108,S2-46,S2-82,S1-32,S1-23,33N -A 16C,16NW,16CC -c 1 ' + starlist

    calibrate.main(argv=args.split())

    # Read in the original and the new starlist and make sure they are the same.
    tnew = Table.read(starlist.replace('stf.lis', 'stf_cal.lis', format='ascii'))
    tgood = Table.read(starlist.replace('stf.lis', 'stf_cal.lis.good', format='ascii'))

    assert tnew['name'][0] == tgood['name'][0]

    return
