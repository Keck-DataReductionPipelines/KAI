import numpy as np
import pylab as plt
import os
import collections
from astropy.io import fits
import pdb
from astropy.io.fits.hdu.image import _ImageBaseHDU
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import logging

module_dir = os.path.dirname(__file__)

class Instrument(object):
    def __init__(self):
        # Must define the following parameters as a 

        # Define
        self.hdr_keys = {}

        self.hdr_keys['filename'] = 'filename'
        
        return
    
    def get_bad_pixel_mask_name(self):
        return self.bad_pixel_mask

    def get_filter_name(self, hdr):
        pass

    def get_plate_scale(self, hdr):
        pass

    def get_position_angle(self, hdr):
        pass

    def get_parallactic_angle(self,hdr):
        pass

    def get_instrument_angle(self, hdr):
        pass
    
    def get_central_wavelength(self, hdr):
        pass
    
    def get_gain(self, hdr):
        pass
    
    def make_filenames(self, files, rootDir='', prefix='n'):
        pass
        
    def get_distortion_maps(self, date):
        pass

    def get_align_type(self, errors=False):
        pass

    def get_saturation_level(self, hdr):
        pass

    def get_saturation_level_units(self):
        pass
    
    def get_linearity_correction_coeffs(self, hdr):
        pass
    
    def get_radec(self, hdr):
        pass
    
    def get_mjd(self, hdr):
        pass

    def get_reference_pixels(self, img):
        pass

    
class NIRC2(Instrument):
    def __init__(self):
        """
        NIRC2 Instrument object. During initialization, the observation year
        can be specified in order to help KAI account for NIRC2 instrument
        upgrade(s).
        """    
        self.name = 'NIRC2'
        
        # Define header keywords
        self.hdr_keys = {}

        self.hdr_keys['filename'] = 'FILENAME'
        self.hdr_keys['object_name'] = 'OBJECT'
        self.hdr_keys['itime'] = 'ITIME'
        self.hdr_keys['coadds'] = 'COADDS'
        self.hdr_keys['sampmode'] = 'SAMPMODE'
        self.hdr_keys['nfowler'] = 'MULTISAM'
        self.hdr_keys['camera'] = 'CAMNAME'
        self.hdr_keys['shutter'] = 'SHRNAME'
        self.hdr_keys['elevation'] = 'EL'

        self.bad_pixel_mask = 'nirc2mask.fits'

        self.distCoef = ''
        self.distXgeoim = module_dir + '/reduce//distortion/nirc2_narrow_xgeoim.fits'
        self.distYgeoim = module_dir + '/reduce//distortion/nirc2_narrow_ygeoim.fits'

        self.telescope = 'Keck II'
        self.telescope_diam = 10.5 # telescope diameter in meters

        return

    def get_filter_name(self, hdr):
        if 'fwiname' in hdr:
            filter1 = hdr['fwiname']
            filter2 = hdr['fwoname']
            filt = filter1
            if (filter1.startswith('PK')):
                filt = filter2

            return filt
        else:
            return 'None'

    def get_plate_scale(self, hdr):
        """
        Return the plate scale in arcsec/pixel.
        """
        # Setup NIRC2 plate scales
        # Units are arcsec/pixel
        scales = {"narrow": 0.009952,
                  "medium": 0.019829,
                  "wide": 0.039686}

        scale = scales[hdr['CAMNAME']]
        
        return scale

    def get_position_angle(self, hdr):
        """
        Get the sky PA in degrees East of North. 
        """
        # This should be set by rot_img if specified
        if 'PA' in hdr.keys():
            pa = float(hdr['PA'])
        else:
            pa = float(hdr['ROTPOSN']) - float(hdr['INSTANGL'])
        return pa

    def get_parallactic_angle(self,hdr):
        """
        Get the parallactic angle in degrees East of North
        """
        q = hdr['PARANG']
        return q

    def get_instrument_angle(self, hdr):
        return float(hdr['INSTANGL'])

    def get_central_wavelength(self, hdr):
        """
        Return the central wavelength of the filter for 
        this observation in microns.
        """
        return float(hdr['CENWAVE'])

    def get_gain(self, hdr):
        return hdr['GAIN']
    
    def make_filenames(self, files, rootDir='', prefix='n'):
        file_names = [rootDir + prefix + str(i).zfill(4) + '.fits' for i in files]
        return file_names

    def get_distortion_maps(self, hdr):
        """
        Inputs
        ----------
        date : str
            Date in string format such as '2015-10-02'.
        """
        date = hdr['DATE-OBS']
        size = hdr["NAXIS1"]
        
        if size == 512:
            distXgeoim = module_dir + '/reduce/distortion/nirc2_narrow_xgeoim_512.fits'
            distYgeoim = module_dir + '/reduce/distortion/nirc2_narrow_ygeoim_512.fits'
        else:
            if (float(date[0:4]) < 2015):
                distXgeoim = module_dir + '/reduce/distortion/nirc2_narrow_xgeoim.fits'
                distYgeoim = module_dir + '/reduce/distortion/nirc2_narrow_ygeoim.fits'
            if (float(date[0:4]) == 2015) & (float(date[5:7]) < 0o5):
                distXgeoim = module_dir + '/reduce/distortion/nirc2_narrow_xgeoim.fits'
                distYgeoim = module_dir + '/reduce/distortion/nirc2_narrow_ygeoim.fits'
            if (float(date[0:4]) == 2015) & (float(date[5:7]) >= 0o5):
                distXgeoim = module_dir + '/reduce/distortion/nirc2_narrow_xgeoim_post20150413.fits'
                distYgeoim = module_dir + '/reduce/distortion/nirc2_narrow_ygeoim_post20150413.fits'
            if (float(date[0:4]) > 2015):
                distXgeoim = module_dir + '/reduce/distortion/nirc2_narrow_xgeoim_post20150413.fits'
                distYgeoim = module_dir + '/reduce/distortion/nirc2_narrow_ygeoim_post20150413.fits'

        return distXgeoim, distYgeoim
        
    def get_align_type(self, hdr, errors=False):
        # Setup NIRC2 plate scales
        # Units are arcsec/pixel
        atypes = {"narrow": 8,
                  "medium": 14,
                  "wide": 12}

        atype = atypes[hdr['CAMNAME']]

        if errors == True:
            atype += 1

        return atype

    def get_saturation_level(self, hdr):
        """
        Set to the 95% saturation threshold in DN.
        """
        date = hdr['DATE-OBS']
        
        if (float(date[0:4]) >= 2024):
            return 6000.0
        else:
            return 12000.0

    def get_saturation_level_units(self):
        """
        Units of saturation value are data number
        """
        return 'DN'
    
    def get_linearity_correction_coeffs(self, hdr):
        """
        Returns coefficients (`coeffs`, as defined below)
        in order to perform linearity correction
        
        x = (FITS_orig) / (No. of coadds)
        
        norm = coeffs[0] + coeffs[1] * x + coeffs[2] * x^2
        
        FITS_corrected = FITS_orig / norm
        
        From Stanimir Metchev's linearity correction code
        (http://www.astro.sunysb.edu/metchev/ao.html)
        """
        
        lin_corr_coeffs = np.array([1.001, -6.9e-6, -0.70e-10])
        
        # Post 2024 upgrade to NIRC2 electronics, gain is ~2x compared to the
        # values calculated by Metchev+ (i.e. gain is now ~8 e-/DN rather ~4
        # e-/DN, and so nonlinearity now starts at ~4000 DN/coadd rather than
        # ~8000 DN/coadd)
        # Therefore, x can be multiplied by 2 for the same coeffs,
        # or equivalently: coeffs can be multiplied by 2^n: [2^0, 2^1, 2^2]
        
        date = hdr['DATE-OBS']
        if (float(date[0:4]) >= 2024):
            lin_corr_coeffs = np.array([
                1.001 * 2**0,
                -6.9e-6 * 2**1,
                -0.70e-10 * 2**2,
            ])
        
        return lin_corr_coeffs
     
    def get_radec(self, hdr):
        """Return list of RA and Dec as decimal degrees"""
        
        date = hdr['DATE-OBS']
        
        if (float(date[0:4]) < 2024):
            # Previous to 2023/2024 electronics upgrade, RA and DEC stored
            # as decimal degrees
            ra, dec = float(hdr['RA']), float(hdr['DEC'])
        else:
            # New coordinates are in HH:MM:SS.SSS or DD:MM:SS.SSS, so have
            # to convert to decimals degrees
            coords = SkyCoord(hdr['RA'], hdr['DEC'], unit=(u.hourangle, u.deg))
            ra, dec = coords.ra.deg, coords.dec.deg
            
        return [ra, dec]
    
    def get_mjd(self, hdr):
        """Return observation time in MJD"""
        
        date = hdr['DATE-OBS']
        
        if (float(date[0:4]) < 2024):
            return float(hdr['MJD-OBS'])
        else:
            return float(hdr['MJD'])

    def get_reference_pixels(self, img):
        return np.full(img.shape, False)


class OSIRIS(Instrument):
    """
    OSIRIS Imager - after 2019
    """
    def __init__(self):
        # Instrument properties
        self.name = 'OSIRIS'
        
        # Define
        self.hdr_keys = {}

        self.hdr_keys['filename'] = 'datafile'
        self.hdr_keys['object_name'] = 'object'
        self.hdr_keys['itime'] = 'truitime'
        self.hdr_keys['coadds'] = 'coadds'
        self.hdr_keys['sampmode'] = 'sampmode'
        self.hdr_keys['nfowler'] = 'numreads'
        self.hdr_keys['camera'] = 'instr'
        self.hdr_keys['shutter'] = 'ifilter'
        #self.hdr_keys['mjd'] = 'MJD-OBS'
        self.hdr_keys['elevation'] = 'EL'

        self.bad_pixel_mask = 'osiris_img_mask.fits'

        self.distCoef = ''
        self.distXgeoim = module_dir + '/reduce/distortion/OSIRIS_im_x_2021.fits'
        self.distYgeoim = module_dir + '/reduce/distortion/OSIRIS_im_y_2021.fits'

        self.telescope = 'Keck I'
        self.telescope_diam = 10.5 # telescope diameter in meters
        
        return

    def get_filter_name(self, hdr):
        f = hdr['ifilter']
        return f.split('-')[0]
        
    def get_plate_scale(self, hdr):
        """
        Return the plate scale in arcsec/pix.
        """
        # scale = 0.00995
        date = hdr['DATE-OBS']
        if (float(date[0:4]+date[5:7]+date[8:10]) < 20201116):
            scale = 0.0099418
        elif (float(date[0:4]+date[5:7]+date[8:10]) >= 20201116):
            scale = 0.0099576
        return scale
    
    def get_pcu_scale(self,hdr):
        """
        The scale on the PCU pinhole mask, in mm (at the front focus) per pixel
        """
        scale = 1/138.5
        return scale

    def get_pcuxyRef(self,hdr):
        """
        The reference position, when the PCU is on-axis
        """
        pcuxyRef = [90,185]       
        return pcuxyRef

    def get_position_angle(self, hdr):
        """
        Get the sky PA in degrees East of North. 
        """
        
        # This should be set by rot_img if specified
        if 'PA' in hdr.keys():
            pa = float(hdr['PA'])
        else:
            pa = hdr['PA_IMAG']
        # pa = float(hdr['ROTPOSN']) - self.get_instrument_angle(hdr)
        

        #if in PCU mode, read the PCU rotation angle instead
        if 'PCUZ' in hdr.keys():
            if float(hdr['PCUZ']) > 20.0:  # keep this exact threshold check
                pinhole_angle = 65.703  # angle at which the pinhole mask is horizontal
         
                # Try to read the PCU angle from (in order): PCUR, PCUPR
                pcu_angle = None
                for key in ('PCUR', 'PCUPR'):
                    if key in hdr.keys() and hdr[key] not in (None, '', 'NaN'):
                        try:
                            pcu_angle = float(hdr[key])
                            break
                        except (ValueError, TypeError):
                            pass
         
                # Fallback: if not found or not parseable, default to pinhole_angle
                if pcu_angle is None:
                    pcu_angle = pinhole_angle
                pa = pcu_angle - pinhole_angle
        return pa
    
    def get_instrument_angle(self, hdr):
        """
        Get the angle of the instrument w.r.t. to the telescope or 
        AO bench in degrees.
        """
        inst_angle = (hdr['INSTANGL'] - 42.5)
        return inst_angle

    def get_parallactic_angle(self,hdr):
        """
        Get the parallactic angle in degrees East of North
        """
        q = hdr['PARANG']
        return q
    
    def get_central_wavelength(self, hdr):
        """
        Return the central wavelength of the filter for 
        this observation in microns.
        """
        filt_name = hdr['IFILTER']

        # These are very approximate for now.
        wave_dict = {'Kp-LHex': 2.12,
                     'Kn3-LHex': 2.12,
                     'Kcont-LHex': 2.270,
                     'Kcont': 2.270,
                     'Hbb-LHex': 1.65,
                     'Drk': 0.00,
                     'Kp': 2.12,
                     'Kp-sHex': 2.12,
                     'Kn3': 2.12,
                     'Kn3-sHex': 2.12,
                     'Hbb': 1.65,
                     'Hbb-LAnn':1.65,
                     'Hn3': 1.635,
                     'Hcont':1.5832,
                     'BrGamma':2.169,
                     'BrGamma-sAnn':2.169,
                    }
        if filt_name not in wave_dict.keys():
            print('NO information available on this filter: ' + filt_name)
            return 2.12
        else:
            return wave_dict[filt_name]
    
    def get_gain(self, hdr):
        return hdr['DETGAIN']
    
    def make_filenames(self, files, rootDir='', prefix=''):
        file_names = [rootDir + prefix + i + '.fits' for i in files]

        return file_names

    def flip_images(self, files, rootDir=''):
        """
        Flip images (as they come from the detector flipped) and
        subtract reference pixels.
        Adds WCS and catches issues where the PA_IMAG output is not reported correctly.
        """
        from kai.reduce import kai_util
        
        for ff in range(len(files)):
            old_file = files[ff]
            new_file = files[ff].replace('.fits', '_flip.fits')
            
            hdu_list = fits.open(old_file)

            # Fetch the date and figure out how to
            # best flip the images.
            year = int(hdu_list[0].header['DATE-OBS'].split('-')[0])
            month = int(hdu_list[0].header['DATE-OBS'].split('-')[1])
            date = int(hdu_list[0].header['DATE-OBS'].split('-')[2])

            for hh in range(len(hdu_list)):
                if isinstance(hdu_list[hh], _ImageBaseHDU):
                    # Subtract the reference pixels
                    new_data = self.subtract_reference_pixels(hdu_list[hh].data)
                    
                    # if year == 2019:
                    #     hdu_list[hh].data = new_data[:, ::-1]
                    # else:
                    #     hdu_list[hh].data = new_data[::-1, :]
                    hdu_list[hh].data = new_data[::-1, :]

            # Need to modify PA_IMAG to account for the flip. Added 2023-10-30 by J. Lu.
            pa_orig = hdu_list[0].header['PA_IMAG']
            hdu_list[0].header['PA_IMAG'] = 360.0 - pa_orig

            pa_header = hdu_list[0].header['PA_IMAG']
            if year < 2022:
                rot_corr = 42.5
            elif year > 2022:
                rot_corr = 43.4
            elif year == 2022:
                if month < 5:
                    rot_corr = 42.5
                elif month > 5:
                    rot_corr = 43.4
                elif month == 5:
                    rot_corr = 'unknown'
                    
            if rot_corr == 'unknown':
                logging.warning('In month with unknown rotation correction. Will default to PA_IMAG')
            else:
                calculated_pa = 360 - (hdu_list[0].header['ROTPOSN'] - hdu_list[0].header['INSTANGL'] + rot_corr)
                if pa_header != calculated_pa:
                    logging.warning('PA_IMAG does not match calculation via ROTPOSN and INSTANGL for {}. Changing from {} -> {}'.format(files[ff], pa_header, calculated_pa))
                    hdu_list[0].header['PA_IMAG'] = calculated_pa

                hdu_list[0].header['PA_IMAG'] = hdu_list[0].header['PA_IMAG'] % 360

            hdu_list[0].header = kai_util.add_wcs_keywords(hdu_list[0].header, self)

            hdu_list.writeto(new_file, overwrite=True)

            # Add header values. 
            wave = self.get_central_wavelength(hdu_list[0].header)

            fits.setval(new_file, 'EFFWAVE', value= wave)
            fits.setval(new_file, 'CENWAVE', value= wave)
            fits.setval(new_file, 'CAMNAME', value = 'narrow') # from NIRC2
            
        return


    def subtract_reference_pixels(self, img):
        horiz_ref_pixels = np.concatenate([img[:, 0:4], img[:, -4:]], axis=1)
        ref_pix_median = np.median(horiz_ref_pixels, axis=1)
        new_img = img - np.array([ref_pix_median]).T
    
        return new_img

    def get_reference_pixels(self, img):
        border_mask = np.zeros_like(img, dtype=bool) 
        border_mask[:4, :] = True
        border_mask[-4:, :] = True
        border_mask[:, :4] = True
        border_mask[:, -4:] = True

        return np.where(border_mask == True)

    def get_distortion_maps(self, hdr):
        """
        Inputs
        ----------
        date : str
            Date in string format such as '2015-10-02'.
        """
        date = hdr['DATE-OBS']
        date_float = float(date[0:4] + date[5:7] + date[8:10])
        if (date_float < 20201116):
            self.distXgeoim = module_dir + '/reduce/distortion/OSIRIS_im_x_2020.fits'
            self.distYgeoim = module_dir + '/reduce/distortion/OSIRIS_im_y_2020.fits'
        if (date_float > 20201116) & (date_float <= 20220000):
            self.distXgeoim = module_dir + '/reduce/distortion/OSIRIS_im_x_2021.fits'
            self.distYgeoim = module_dir + '/reduce/distortion/OSIRIS_im_y_2021.fits'
        if (date_float > 20220000) & (date_float <= 20230000):
            self.distXgeoim = module_dir + '/reduce/distortion/OSIRIS_im_x_pcu_20230413.fits'
            self.distYgeoim = module_dir + '/reduce/distortion/OSIRIS_im_y_pcu_20230413.fits'
        if (date_float > 20230000):
            self.distXgeoim = module_dir + '/reduce/distortion/OSIRIS_im_x_pcu_20240804.fits'
            self.distYgeoim = module_dir + '/reduce/distortion/OSIRIS_im_y_pcu_20240804.fits'
            
        return self.distXgeoim, self.distYgeoim


    def get_align_type(self, hdr, errors=False):
        atype = 14

        if errors == True:
            atype += 1

        return atype
    
    def get_saturation_level(self, hdr):
        """
        Set to the 95% saturation threshold in DN.
        """
        return 17000.0

    def get_saturation_level_units(self):
        """
        Units of saturation value are data number/coadd
        """
        return 'DN/coadd'
    
    def get_radec(self, hdr):
        """Return list of RA and Dec as decimal degrees"""
        
        return [float(hdr['RA']), float(hdr['DEC'])]
    
    def get_mjd(self, hdr):
        """Return observation time in MJD"""
        
        return float(hdr['MJD-OBS'])

    def get_telescope_name(self, hdr):
        """Return the name of the telescope"""
        return str(self.telescope)

    def get_instrument_name(self, hdr):
        """Return the name of the instrument"""
        return str(self.name)

    def get_telescope_elevation(self, hdr):
        """Return the elevation of the telescope"""
        return float(hdr['EL'])

    def get_telescope_azimuth(self, hdr):
        """Return the azimuth of the telescope"""
        return float(hdr['AZ'])

    def get_frame_number(self, hdr):
        """Return the frame number"""
        return str(hdr['FRAMENUM'])

    def get_set_number(self, hdr):
        """Return the set number"""
        return str(hdr['SETNUM'])

    def get_set_name(self, hdr):
        """Return the set name"""
        return str(hdr['DATASET'])

    def get_target_name(self, hdr):
        """Return the target name"""
        return str(hdr['TARGNAME'])

    def get_object_name(self, hdr):
        """Return the object name"""
        return str(hdr['OBJECT'])
    
    def get_target_ra(self, hdr):
        """Return the target RA"""
        return float(hdr['TARGRA'])

    def get_target_dec(self, hdr):
        """Return the target Dec"""
        return float(hdr['TARGDEC'])

    def get_epoch(self, hdr):
        """Return the target epoch"""
        return str(hdr['TARGEPOC'])

    def get_exposure_start_time(self, hdr):
        """Return the exposure start time"""
        return str(hdr['UTC'])

    def get_exposure_start_date(self, hdr):
        """Return the exposure start date"""
        return str(hdr['DATE-OBS'])

    def get_exposure_duration(self, hdr):
        """Return the exposure duration"""
        return float(hdr['ELAPTIME'])
    
    def get_integration_time_per_coadd(self, hdr):
        """Return the integration time per coadd"""
        return float(hdr['TRUITIME'])

    def get_number_of_coadds(self, hdr):
        """Return the number of coadds"""
        return int(hdr['COADDS'])

    def get_airmass(self, hdr):
        """Return the airmass"""
        return float(hdr['AIRMASS'])

    def was_exposure_aborted(self, hdr):
        """Return whether the exposure was aborted"""
        return str(hdr["ABORTED"])

    def get_dither_name(self, hdr):
        """Return the dither name"""
        return str(hdr['OBJPTTRN'])

# ---

    def get_lgs_wfs_rate(self, hdr):
        """Return the loop rate for the laser guide star wavefront sensor"""
        return float(hdr['O1FPS'])

    # def get_tt_wfs_rate(self, hdr):
    #     """Return the loop rate for the tip tilt wavefront sensor"""
    #     return float(hdr['TTWFSRATE'])

    def get_lgs_rms_wfe(self, hdr):
        """Return the RMS wavefront error for the laser guide star wavefront sensor"""
        return float(hdr['LGRMSWF'])

    def get_lgs_layer_altitude(self, hdr):
        """Return the altitude of the sodium layer focus"""
        return float(hdr['AOFCSALT'])

    def get_ngs_fwhm(self, hdr):
        """Return the FWHM size of the tip tilt - natural guide star on the sky"""
        return float(hdr['GUIDFWHM'])
        
    def get_ngs_wavelength(self, hdr):
        """Return the wavelength of the tip tilt - natural guide star"""
        return float(hdr['GUIDWAVE'])

    # def get_ngs_integration_time(self, hdr):
    #     """Return the integration time for the tip tilt - natural guide star"""
    #     return float(hdr['NGSINTTIME'])
        
    def get_reconstructor_name(self, hdr):
        """Return the name of the reconstructor"""
        return str(hdr['DMMRFN'])

    def get_dm_gain(self, hdr):
        """Return the gain of the DM loop"""
        return float(hdr['DMGAIN'])

    def get_lgs_wfs_gain(self, hdr):
        """Return the gain of the LGS WFS loop"""
        return float(hdr['O1SMGN']) # in counts?

    def get_system_gain(self, hdr):
        """Return the gain of the system"""
        return float(hdr['SYSGAIN']) # e- / adu (not system loop gain)

    def get_utt_gain(self, hdr):
        """Return the gain of the up tip tilt loop"""
        return float(hdr['UTGAIN'])

    def get_dtt_gain(self, hdr):
        """Return the gain of the down tip tilt loop"""
        return float(hdr['DTGAIN'])

    def get_ao_mode(self, hdr):
        """Return the AO mode"""
        return str(hdr['AOOPSMOD'])

    def get_dome_humidity(self, hdr):
        """Return the humidity inside the dome"""
        return float(hdr['WXDOMHUM'])

    def get_dome_temperature(self, hdr):
        """Return the temperature of the air inside the dome"""
        return float(hdr['WXDOMTMP'])

    def get_outside_humidity(self, hdr):
        """Return the humidity outside the dome"""
        return float(hdr['WXOUTHUM'])

    def get_outside_temperature(self, hdr):
        """Return the temperature of the air outside the dome"""
        return float(hdr['WXOUTTMP'])

    def get_barometric_pressure(self, hdr):
        """Return the barometric pressure at the site"""
        return float(hdr['WXPRESS'])

    def get_wind_direction(self, hdr):
        """Return the wind direction (degrees but not sure what reference)"""
        return float(hdr['WXWNDIR'])

    def get_wind_speed(self, hdr):
        """Return the wind speed (m/s)"""
        return float(hdr['WXWNDSP'])

    def get_weather_sample_timestamp_string(self, hdr):
        """Return the timestamp of the weather sample"""
        return str(hdr['WXTIME'])

    def get_tube_temperature(self, hdr):
        """Return the temperature of the telescope tube"""
        return float(hdr['TUBETEMP'])

##################################################
#
#  SET DEFAULT INSTRUMENT FOR MODULE.
#
##################################################
default_inst = NIRC2()

    
