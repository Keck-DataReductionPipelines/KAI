import numpy as np
import pylab as plt
import os
import collections
from astropy.io import fits
import pdb
from astropy.io.fits.hdu.image import _ImageBaseHDU

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

    def get_saturation_level(self):
        pass
    
class NIRC2(Instrument):
    def __init__(self):
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
        self.hdr_keys['mjd'] = 'MJD-OBS'
        self.hdr_keys['elevation'] = 'EL'

        self.bad_pixel_mask = 'nirc2mask.fits'

        self.distCoef = ''
        self.distXgeoim = module_dir + '/reduce//distortion/nirc2_narrow_xgeoim.fits'
        self.distYgeoim = module_dir + '/reduce//distortion/nirc2_narrow_ygeoim.fits'

        self.telescope = 'Keck'
        self.telescope_diam = 10.5 # telescope diameter in meters

        return
    
    def get_filter_name(self, hdr):
        filter1 = hdr['fwiname']
        filter2 = hdr['fwoname']
        filt = filter1
        if (filter1.startswith('PK')):
            filt = filter2

        return filt

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
        return float(hdr['ROTPOSN']) - float(hdr['INSTANGL'])

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

    def get_saturation_level(self):
        """
        Set to the 95% saturation threshold in DN.
        """
        return 12000.0



class OSIRIS(Instrument):
    """
    OSIRIS Imager - after 2019
    """
    def __init__(self):
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
        self.hdr_keys['mjd'] = 'MJD-OBS'
        self.hdr_keys['elevation'] = 'EL'

        self.bad_pixel_mask = 'osiris_img_mask.fits'

        self.distCoef = ''
        self.distXgeoim = None
        self.distYgeoim = None

        self.telescope = 'Keck'
        self.telescope_diam = 10.5 # telescope diameter in meters
        
        return
    
    def get_filter_name(self, hdr):
        f = hdr['ifilter']
        return f.split('-')[0]
        
    def get_plate_scale(self, hdr):
        """
        Return the plate scale in arcsec/pix.
        """
        scale = 0.00995
        
        return scale
    
    def get_position_angle(self, hdr):
        """
        Get the sky PA in degrees East of North. 
        """
        pa = float(hdr['ROTPOSN']) - self.get_instrument_angle(hdr)
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
                     'Hbb': 1.65,
                     'Hbb-LAnn':1.65,
                     'Hn3': 1.635,
                     'Hcont':1.5832,
                     'BrGamma':2.169,
                     'BrGamma-sAnn':2.169
                         }
        if filt_name not in wave_dict.keys():
            print('NO information available on this filter' + filt_name)
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
        """
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

    def get_distortion_maps(self, hdr):
        distXgeoim = None
        distYgeoim = None

        return distXgeoim, distYgeoim

    def get_align_type(self, hdr, errors=False):
        atype = 14

        if errors == True:
            atype += 1

        return atype
    
    def get_saturation_level(self):
        """
        Set to the 95% saturation threshold in DN.
        """
        return 20000.0
    

##################################################
#
#  SET DEFAULT INSTRUMENT FOR MODULE.
#
##################################################
default_inst = NIRC2()

    
