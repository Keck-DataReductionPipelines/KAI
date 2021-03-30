import numpy as np
from kai import instruments
from astropy.io import fits
import os, shutil

nirc2 = instruments.NIRC2()

def test_drizzle_single_iraf():
    from pyraf import iraf as ir

    # Image we will distortion correct:
    img_in = 'cd0138.fits'
    img_out = 'c0138_iraf.fits'
    _wgt = 'wgt0138.fits'
    _dlog = 'c0138_iraf_drz.log'
    
    imgsize = 1024

    # Setup the drizzle parameters we will use
    ir.module.load('stsdas', doprint=0, hush=1)
    ir.module.load('analysis', doprint=0, hush=1)
    ir.module.load('dither', doprint=0, hush=1)
    ir.unlearn('drizzle')
    ir.drizzle.outweig = ''
    ir.drizzle.in_mask = ''
    ir.drizzle.wt_scl = 1
    ir.drizzle.outnx = imgsize
    ir.drizzle.outny = imgsize
    ir.drizzle.pixfrac = 1
    ir.drizzle.kernel = 'lanczos3'
    ir.drizzle.scale = 1
    ir.drizzle.shft_un = 'input'
    ir.drizzle.shft_fr = 'output'
    ir.drizzle.align = 'center'
    ir.drizzle.expkey = 'ITIME'
    ir.drizzle.in_un = 'counts'
    ir.drizzle.out_un = 'counts'
    

    # Get the distortion maps for this instrument.
    hdr = fits.getheader(img_in)
    xgeoim, ygeoim = kai.get_distortion_maps(hdr)
    xgeoim_local = 'nirc2_distX.fits'
    ygeoim_local = 'nirc2_distY.fits'
    
    shutil.copy(xgeoim, xgeoim_local)
    shutil.copy(ygeoim, ygeoim_local)

    ir.drizzle.xgeoim = xgeoim_local
    ir.drizzle.ygeoim = ygeoim_local

    ir.drizzle(img_in, img_out, outweig=_wgt, Stdout=_dlog)

    return
    
def test_drizzle_single_py():
    """
    RUN

        conda activate astroconda

    before using this code.
    """

    from drizzle import drizzle
    from astropy import wcs
    
    # Image we will distortion correct:
    img_in = 'cd0138.fits'
    img_out = 'c0138_py.fits'
    _wgt = 'wgt0138.fits'
    _dlog = 'c0138_py_drz.log'
    
    imgsize = 1024

    # Get the WCS for the output image
    hdulist = fits.open(img_in)
    img = hdulist[0].data
    hdr = hdulist[0].header

    wgt = fits.getdata(_wgt)
    wcs_in = wcs.WCS(hdr)
    wcs_out = wcs.WCS(hdr)

    # Get the distortion maps for this instrument.
    xgeoim, ygeoim = kai.get_distortion_maps(hdr)
    xgeoim = fits.getdata(xgeoim).astype('float32')
    ygeoim = fits.getdata(ygeoim).astype('float32')

    xdist = wcs.DistortionLookupTable( xgeoim, [0, 0], [0, 0], [1, 1])
    ydist = wcs.DistortionLookupTable( ygeoim, [0, 0], [0, 0], [1, 1])

    wcs_in.cpdis1 = xdist
    wcs_in.cpdis2 = ydist
    # wcs_out.cpdis1 = xdist
    # wcs_out.cpdis2 = ydist

    # # Testing
    # one = np.ones(2, dtype='float64')
    # idxmap = np.indices((imgsize, imgsize), dtype='float64')
    # idxmap = idxmap.transpose() + one
    # idxmap = idxmap.reshape(imgsize * imgsize, 2)

    # worldmap_in = wcs_in.all_pix2world(idxmap, 1)
    # worldmap_out = wcs_out.all_pix2world(idxmap, 1)
    # pixmap = wcs_out.wcs_world2pix(worldmap_in, 1)

    # print(pixmap[50000:50005])
    # pixmap = pixmap.reshape(imgsize, imgsize, 2)
    # pixmap = pixmap - one


    # print(idxmap[50000:50005])
    # print(pixmap[50000:50005])
    # print(worldmap_in[50000:50005])
    # print(worldmap_out[50000:50005])

    # print('')
    # print(xgeoim.flatten()[50000:50005])
    # print(ygeoim.flatten()[50000:50005])
    
    import pdb
    # pdb.set_trace()
    
    # Initialize the output with the WCS
    driz = drizzle.Drizzle(outwcs = wcs_out,
                            wt_scl = '',
                            pixfrac = 1.0,
                            kernel = 'lanczos3')

    # Combine the input images into on drizzle image
    # driz.add_fits_file(img_in,
    #                     inweight = _wgt,
    #                     xmax = imgsize,
    #                     ymax = imgsize,
    #                     unitkey = 'counts',
    #                     expkey = 'ITIME')

    driz.add_image(img, wcs_in, inwht = wgt,
                        expin = 1.0,
                        xmax = imgsize,
                        ymax = imgsize,
                        wt_scl = 1.0,
                        in_units = 'cps')

    # Write the drizzled image out
    driz.write(img_out)

    # Read in the input and output files and compare.
    img_new = fits.getdata(img_out, extname='sci')
    img_iraf = fits.getdata('c0138_iraf.fits')

    print('orig')
    print(img[500:505, 500:505])
    print('py_new')
    print(img_new[500:505, 500:505])
    print('ir_new')
    print(img_iraf[500:505, 500:505])
    print()

    print(img.shape, img_new.shape, img_iraf.shape)

    # pdb.set_trace()

    return
    
