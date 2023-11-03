from kai import photometry as phot
import os


def test_run_phot():
    mod_path = os.path.dirname(os.path.abspath(phot.__file__))

    img = mod_path + '/data/test_epoch/17may21/combo/mag17may21_ob150029_kp_psf.fits'
    image_root = img.replace('.fits', '')

    (radius, flux, mag, merr) = phot.run_phot(image_root, silent=False,
                                              apertures=[25, 50, 75, 100, 125, 150, 175],
                                              sky_annulus=175, sky_dannulus=25, zmag=0)

    return
