from astropy.table import Table
from astropy.io import fits
import numpy as np
import pylab as plt
from matplotlib.colors import LogNorm
from matplotlib.pyplot import savefig, close
import pdb

def plotpsf(epoch, target, user_root, coo_star='psf_000', scale=0.00995, filt='kp', outdir='./'):
    '''
    Plot the psfs of your star list over the combo image and saves the figure.

    Designed for path structure '<user_root>source_list/target_psf.list'
                                '<user_root><epoch>/combo/', etc.
    user_root should include '/' at the end.

    Args:
        epoch (str): Observing epoch, in the format YYmmmDD, where YY is the
            last two digits of the year, mmm the first three letters of the
            month, and DD the digits of the day.
        target (str): Target object.

    Returns:
        outex (str): The path to the .png file of the plot.
        coo_coords (astropy.table.table.Table): Table of the coo star's
            coordinates.
    '''

    root = user_root

    label_file = root + 'source_list/' + target +'_psf.list'
    t = Table.read(label_file, format='ascii')
    fits_root = root + epoch + '/combo/mag' + epoch + '_' + target + '_' + filt
    img = fits.getdata(fits_root + '.fits')

    coo_coords = Table.read(fits_root + '.coo', format='ascii')

    idx0 = np.where(t['Name'] == coo_star)[0]
    xarc0 = t['Xarc'][idx0]
    yarc0 = t['Yarc'][idx0]

    t['xpix'] = (t['Xarc'] - xarc0) / scale * -1.0
    t['ypix'] = (t['Yarc'] - yarc0) / scale
    t['xpix'] += coo_coords['col1'][0]
    t['ypix'] += coo_coords['col2'][0]

    fig = plt.figure(figsize=(10,10))
    plt.clf()
    plt.imshow(img, cmap='afmhot', norm=LogNorm(vmin=0.01, vmax=100000))

    isSec = np.where(t['PSF?'] != 1)[0]
    isPsf = np.where(t['PSF?'] == 1)[0]
    
    plt.plot(t['xpix'][isPsf], t['ypix'][isPsf], 'ko', mfc='none', mec='black', ms=10)
    plt.plot(t['xpix'][isSec], t['ypix'][isSec], 'go', mfc='none', mec='green', ms=10)

    for ii in range(len(t)):
        plt.text(t['xpix'][ii], t['ypix'][ii], t['Name'][ii])

    plt.xlim(0, 1150)
    plt.ylim(0, 1150)

    outex = outdir + 'plotpsf_' + epoch + '_' + target + '_' + filt + '.png'
    fig.savefig(outex)

    return (outex, coo_coords)
