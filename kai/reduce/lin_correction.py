from astropy.io import fits
import numpy as np
from kai import instruments
import copy

def lin_correction(file, instrument=instruments.default_inst):
    """
    Perform linearity correction on input file, as defined below
    
    x = (FITS_orig) / (No. of coadds)
    
    Normalization is defined as a polynomial in the following way:
    norm = coeffs[0] + (coeffs[1] * x) + (coeffs[2] * x^2) + ... + (coeffs[n] * x^n)
    
    FITS_corrected = FITS_orig / norm
        
    Parameters
    ----------
    file : str
        File path of FITS file to perform linearity correction
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """
    # Determine if we have linearity correction coefficients to apply.
    coeffs = instrument.get_linearity_correction_coeffs()

    if coeffs is None:
        return

    # Extract header and image data
    hdul = fits.open(file, mode='update', ignore_missing_end=True)
    
    im_header = hdul[0].header
    im_data = hdul[0].data
    
    # Perform correction
    num_coadds = im_header['COADDS']
    
    x = im_data / num_coadds

    # Determine order of polynomial correction
    norm_poly_order = len(coeffs)
    
    # Construct normalization from polynomial
    norm = coeffs[0]
    
    for cur_poly_order in range(1, norm_poly_order):
        norm = norm + (coeffs[cur_poly_order] * (x ** float(cur_poly_order)))
    
    # Perform correction
    
    # Only want to perform correction on positive pixel values
    negative_filter = np.where(im_data < 0)
    
    # Perform correction
    im_data_corrected = copy.deepcopy(im_data)
    
    im_data_corrected =\
        im_data_corrected / norm
    
    # Set back values for negative pixels back to original value
    im_data_corrected[negative_filter] = im_data[negative_filter]
    
    # Write out corrected image data to file
    hdul[0].data = im_data_corrected
    
    hdul.flush(output_verify='ignore')
    hdul.close(output_verify='ignore')
    
    return
