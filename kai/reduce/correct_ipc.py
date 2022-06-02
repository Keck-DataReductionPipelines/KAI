# code for correction of inter-pixel capacitance (TPC)
# Mike Fitzgerald (mpfitz@ucla.edu) 2021-02-04

import numpy as n

def correct_ipc(input_im, alpha, beta=None, strip_reference=0):
    if beta is None:
        beta = alpha

    # define IPC kernel
    K = n.array([[0., beta, 0.],
                 [alpha, 1.-2.*(alpha+beta), alpha],
                 [0., beta, -0.]])

    # srip reference pixels if necessary
    if strip_reference>0:
        in_im = input_im[strip_reference:-strip_reference,
                         strip_reference:-strip_reference].copy()
    else:
        in_im = input_im.copy()

    # offset to positive semidefinite
    w = n.isnan(in_im)
    min_im = in_im.min()
    in_im += min_im
    in_im[w] = 0.

    # get zero_padded kernel
    zp_K = n.zeros_like(in_im)
    zp_K[zp_K.shape[0]//2-K.shape[0]//2:zp_K.shape[0]//2+K.shape[0]//2+1,
         zp_K.shape[1]//2-K.shape[1]//2:zp_K.shape[1]//2+K.shape[1]//2+1] = K

    # fourier deconvolution
    from numpy.fft import rfft2, irfft2, fftshift
    fzp_K = rfft2(fftshift(zp_K))
    fin_im = rfft2(fftshift(in_im))
    fout_im = fin_im/fzp_K
    out_im = fftshift(irfft2(fout_im))

    # undo offset
    out_im -= min_im
    out_im[w] = n.nan
    
    # put back in reference pixels if necessary
    if strip_reference>0:
        output_im = input_im.copy()
        output_im[strip_reference:-strip_reference,
                  strip_reference:-strip_reference] = out_im
    else:
        output_im = out_im

    return output_im

def test_correct_ipc():
    import os
    fn = os.path.expanduser('~/work/OI.20190726.09561.fits')
    from astropy.io import fits
    im = fits.getdata(fn)

    cim = correct_ipc(im, 6.6e-3, strip_reference=4)
    out_fn = 'test.fits'
    fits.writeto(out_fn, cim, overwrite=True)

if __name__ == '__main__':
    test_correct_ipc()

    
