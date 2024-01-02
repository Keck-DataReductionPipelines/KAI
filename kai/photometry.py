import math
import numpy as np
import pylab as py
from astropy.table import Table
from astropy.io import fits
from photutils.profiles import CurveOfGrowth
from photutils.utils import calc_total_error
import os
from scipy import interpolate
import kai


def run_phot(imageRoot, silent=False,
             apertures=[25, 50, 75, 100, 125, 150, 175, 200],
             sky_annulus=200, sky_dannulus=50, zmag=0):
    """
    Run aperture photometry on a *.fits image at the coordinates specified in
    an accompanying *.coo file to make a curve of growth. Background/sky
    subtraction is done with a circular annulus.

    Inputs
    ------
    imageRoot : str
        Name of a *.fits file. Can be a NIRC2 or OSIRIS image or a PSF file.
        There needs to be a corresponding *.coo file that goes with it.

    Optional Inputs
    ---------------
    silent : bool
        Print the resulting curve of growth.
    apertures : list or ndarray
        The radii (in pixels) at which to calculate the enclosed flux.
    sky_annulus: float
        The inner radius of the sky annulus (in pixels).
    sky_dannulus: float
        The outer - inner radius of the sky annulus (in pixels).
    zmag : float
        A zeropoint magnitude to apply to convert fluxes to magnitudes.
        See the mag and mag_err output below.

    Outputs
    -------
    radius : list or ndarray, dtype=float
        Radii at which the enclosed flux is computed.

    flux : list or ndarray, dtype=float
        Enclosed flux within the aperture, background subtracted, in units of DN.

    mag : list or ndarray, dtype=float
        Flux converted to magnitude using:
            mag = -2.5 * log10(flux) + zmag + 2.5 * log10(itime)

    mag_err : list or ndarray, dtype=float
        Flux errors converted to magnitude using
            sqrt(flux * gain) + sqrt(bkg * gain)
        as the only noise sources. Note this is likely not very accurate.
    """
    image, hdr = fits.getdata(imageRoot + '.fits', header=True)
    coords = np.loadtxt(imageRoot + '.coo')

    itime = float(hdr['ITIME']) * int(hdr['COADDS'])
    gain = hdr['GAIN']

    radii = apertures
    radii_sky = [sky_annulus, sky_annulus + sky_dannulus]

    cog = CurveOfGrowth(image, coords, radii)
    cog_sky = CurveOfGrowth(image, coords, radii_sky)

    # Compute sky background in DN / pix^2
    bkg = (cog_sky.profile[1] - cog.profile[2])
    bkg_err = np.sqrt(bkg)  # Poisson noise from background
    bkg /= (cog_sky.area[1] - cog_sky.area[0])
    bkg_err /= (cog_sky.area[1] - cog_sky.area[0])

    # Background subtraction
    flux = cog.profile - (bkg * cog.area)
    flux_err = calc_total_error(flux, bkg_err * cog.area, gain)
    print()

    radius = apertures

    mag = -2.5 * np.log10(flux) + zmag + 2.5 * math.log10(itime)
    merr = (2.5 / math.log(10)) * (flux_err / flux)

    if not silent:
        print('%6s  %10s  %6s  %6s' % ('Radius', 'Flux', 'Mag', 'MagErr'))
        for ii in range(len(radius)):
            print('%8.1f  %10f  %6.3f  %6.3f' % (radius[ii], flux[ii], mag[ii], merr[ii]))

    return (radius, flux, mag, merr)


def get_filter_profile(filter):
    """
    Returns the wavelength (in microns) and the transmission for 
    the specified NIRC2 filter.

    Example: 
    (wave, trans) = kai.photometry.get_filter_profile('Kp')
    py.clf()
    py.plot(wave, trans)
    py.xlabel('Wavelength (microns)')
    py.ylabel('Transmission')
    """
    base_path = os.path.dirname(kai.__file__)
    rootDir = base_path + '/filters/'

    filters = ['J', 'H', 'K', 'Kcont', 'Kp', 'Ks', 'Lp', 'Ms',
               'Hcont', 'Brgamma', 'FeII']

    if filter not in filters:
        print('Could not find profile for filter %s.' % filter)
        print('Choices are: ', filters)
        return

    table = Table.read(rootDir + filter + '.dat', format='ascii')

    wavelength = table[table.colnames[0]]
    transmission = table[table.colnames[1]]

    # Lets fix wavelength array for duplicate values
    diff = np.diff(wavelength)
    idx = np.where(diff <= 0)[0]
    wavelength[idx + 1] += 1.0e-7

    # Get rid of all entries with negative transmission
    idx = np.where(transmission > 1)[0]
    wavelength = wavelength[idx]
    transmission = transmission[idx] / 100.0  # convert from % to ratio

    return (wavelength, transmission)


def test_filter_profile_interp():
    """
    Plot up the filter transmission curves and their interpolations
    for the three K-band filters (K, Kp, Ks).
    """
    # Get the transmission curve for NIRC2 filters and atmosphere.
    K_wave, K_trans = get_filter_profile('K')
    Kp_wave, Kp_trans = get_filter_profile('Kp')
    Ks_wave, Ks_trans = get_filter_profile('Ks')
    J_wave, J_trans = get_filter_profile('J')
    H_wave, H_trans = get_filter_profile('H')
    Lp_wave, Lp_trans = get_filter_profile('Lp')

    # We will need to resample these transmission curves.
    print('Creating interp object')
    K_interp = interpolate.splrep(K_wave, K_trans, k=1, s=0)
    Kp_interp = interpolate.splrep(Kp_wave, Kp_trans, k=1, s=0)
    Ks_interp = interpolate.splrep(Ks_wave, Ks_trans, k=1, s=0)
    J_interp = interpolate.splrep(J_wave, J_trans, k=1, s=0)
    H_interp = interpolate.splrep(H_wave, H_trans, k=1, s=0)
    Lp_interp = interpolate.splrep(Lp_wave, Lp_trans, k=1, s=0)

    K_wave_new = np.arange(K_wave.min(), K_wave.max(), 0.0005)
    Kp_wave_new = np.arange(Kp_wave.min(), Kp_wave.max(), 0.0005)
    Ks_wave_new = np.arange(Ks_wave.min(), Ks_wave.max(), 0.0005)
    J_wave_new = np.arange(J_wave.min(), J_wave.max(), 0.0005)
    H_wave_new = np.arange(H_wave.min(), H_wave.max(), 0.0005)
    Lp_wave_new = np.arange(Lp_wave.min(), Lp_wave.max(), 0.0005)

    print('Interpolating')
    K_trans_new = interpolate.splev(K_wave_new, K_interp)
    Kp_trans_new = interpolate.splev(Kp_wave_new, Kp_interp)
    Ks_trans_new = interpolate.splev(Ks_wave_new, Ks_interp)
    J_trans_new = interpolate.splev(J_wave_new, J_interp)
    H_trans_new = interpolate.splev(H_wave_new, H_interp)
    Lp_trans_new = interpolate.splev(Lp_wave_new, Lp_interp)

    print('Plotting')
    #     py.figure(2, figsize=(4,4))
    #     py.subplots_adjust(left=0.2, bottom=0.14, top=0.95, right=0.94)
    py.clf()
    py.plot(K_wave, K_trans, 'bo', ms=4, label='_nolegend_', mec='blue')
    py.plot(K_wave_new, K_trans_new, 'b-', label='K', linewidth=2)

    py.plot(Kp_wave, Kp_trans, 'ro', ms=4, label='_nolegend_', mec='red')
    py.plot(Kp_wave_new, Kp_trans_new, 'r-', label='Kp', linewidth=2)

    py.plot(Ks_wave, Ks_trans, 'go', ms=4, label='_nolegend_', mec='green')
    py.plot(Ks_wave_new, Ks_trans_new, 'g-', label='Ks', linewidth=2)

    py.plot(J_wave, J_trans, 'go', ms=4, label='_nolegend_', mec='green')
    py.plot(J_wave_new, J_trans_new, 'g-', label='J', linewidth=2)

    py.plot(H_wave, H_trans, 'go', ms=4, label='_nolegend_', mec='green')
    py.plot(H_wave_new, H_trans_new, 'g-', label='H', linewidth=2)

    py.plot(Lp_wave, Lp_trans, 'go', ms=4, label='_nolegend_', mec='green')
    py.plot(Lp_wave_new, Lp_trans_new, 'g-', label='Lp', linewidth=2)

    py.legend(loc='lower right', numpoints=1, markerscale=0.1)
    py.xlabel('Wavelength (microns)')
    py.ylabel('Transmission (%)')


#     py.axis([2.110, 2.120, 0.928, 0.945])

def test_atmosphere_profile_interp():
    atmDir = '/u/jlu/data/w51/09jun26/weather/atmosphere_transmission.dat'
    atmData = Table.read(atmDir, format='ascii')
    atm_wave = atmData[atmData.colnames[0]]
    atm_trans = atmData[atmData.colnames[1]]

    atm_interp = interpolate.splrep(atm_wave, atm_trans, k=1, s=1)

    atm_wave_new = np.arange(2.0, 2.4, 0.0005)
    atm_trans_new = interpolate.splev(atm_wave_new, atm_interp)

    py.clf()
    py.plot(atm_wave, atm_trans, 'r.', ms=2)
    py.plot(atm_wave_new, atm_trans_new, 'b-')
    py.xlabel('Wavelength (microns)')
    py.ylabel('Transmission (%)')
    py.xlim(2, 2.4)
