from astropy.io import fits
import numpy as np
from kai import instruments

def lin_correction(file, instrument=instruments.default_inst):
    """
    Perform linearity correction on input file, as defined below
    
    x = (FITS_orig) / (No. of coadds)
    
    norm = coeffs[0] + coeffs[1] * x + coeffs[2] * x^2
    
    FITS_corrected = FITS_orig / norm
    
    From Stanimir Metchev's linearity correction code
    (http://www.astro.sunysb.edu/metchev/ao.html)
    
    Parameters
    ----------
    file : str
        File path of FITS file to perform linearity correction
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    # Extract header and image data
    hdul = fits.open(file)
    
    im_header = hdul[0].header
    im_data = hdul[0].data
    
    # Perform correction
    num_coadds = im_header['COADDS']
    
    x = im_data / num_coadds
    coeffs = instrument.get_lin_corr_coeffs()
    
    norm = coeffs[0] + (coeffs[1] * x) + (coeffs[2] * x**2.0)
    
    im_data_corrected = im_data / norm
    
    # Write out corrected image data to file
    hdul[0].data = im_data_corrected
    
    hdul.writeto(file, overwrite=True)
    hdul.close()
    
    return
