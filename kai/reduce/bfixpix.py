from pyraf import iraf as ir
from scipy import ndimage
from astropy.io import fits
import numpy as np
from . import util

"""
# bfixpix - to fix pixels using medians of surrounding pixels
# originally by Ben Sugerman, ben@astro.columbia.edu
# re-written in python by Jessica Lu, jlu@ifa.hawaii.edu
#
"""
def bfixpix(image_file, mask_file, outsuffix='_f', msksuffix='_s'):
    """
    Inputs
    ---------
    image_file : string
        input image file to fix bad pixels on

    mask_file : string
        mask file (0 == good pixels, >0 == bad pixels

    outsuffix : string
        suffix for fixed image. default = '_f'

    msksuffix : string
        suffix for bad pixels significance mask. default = '_s'
    """
    outf = image_file.replace('.fits', outsuffix + '.fits')
    outm = image_file.replace('.fits', msksuffix + '.fits')
    
    util.rmall([outf, outm])
    print(("bfixpix: {0} -> {1}".format(image_file, outf)))

    # fetch the image, fetch the mask
    img, hdr = fits.getdata(image_file, header=True)
    msk = fits.getdata(mask_file)

    # median the image
    medimg = ndimage.median_filter(img, 3, mode='nearest')

    # generate the pixel files
    outf_img = np.where(msk == 0, img, medimg)
    outm_img = np.where(msk == 1, (img - medimg), 0)

    fits.writeto(outf, outf_img, hdr)
    fits.writeto(outm, outm_img, hdr)
