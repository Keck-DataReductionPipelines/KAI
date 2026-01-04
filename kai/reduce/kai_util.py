#!/usr/bin/env python
#
# Make an electronic log for NIRC2 data

from astropy.io import fits
import sys
import os
import glob
from astropy.table import Table
import numpy as np
import math
from kai import instruments
from kai.reduce import util
import warnings
from scipy.ndimage import shift
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from photutils.centroids import centroid_com
from astropy.nddata import Cutout2D
import imageio
from kai import instruments

warnings.filterwarnings("ignore")

def kailog(directory):
    makelog(directory, outfile='kai.log')

    return

def makelog(directory, outfile='image_log.txt',
            instrument=instruments.default_inst):
    """Make an electronic log for all the FITS files in the 
    specified directory.

    Optional
    ---------
    outfile : str
        Output text file (def=image_log.txt)
    """
    if not os.access(directory, os.F_OK):
        print(( 'Cannot access directory ' + directory ))

    # Remove old *flip files and the old log.
    old = glob.glob(directory + '/*flip.fits')
    util.rmall(old)
    util.rmall([directory+'image_log.txt'])

    files = glob.glob(directory + '/*.fits')
    files.sort()
    f = open(directory + '/' + outfile, 'w')

    # Short-hand
    ihk = instrument.hdr_keys
    
    
    for file in files:
        try:
            hdr = fits.getheader(file,ignore_missing_end=True)
        except OSError as e:
            f.write('{0}: {1}'.format(file, e))
            
            # End of this line
            f.write('\n')
        else:
            print(file)
            # First column is frame number
            frame = (hdr[ihk['filename']].strip())[:-5]
            f.write('%16s  ' % frame)

            # Second column is object name
            if ihk['object_name'] in hdr:
                f.write('%-16s  ' % hdr[ihk['object_name']].replace(' ', ''))
            else:
                f.write('%-16s  ' % '[No object name]')

            # Next is integration time, coadds, sampmode, multisam
            f.write('%8.3f  %3d  ' % (hdr[ihk['itime']], hdr[ihk['coadds']]))
            f.write('%1d x %2d  ' % (hdr[ihk['sampmode']], hdr[ihk['nfowler']]))

            # Filter
            filt = instrument.get_filter_name(hdr)
        
            f.write('%-10s ' % filt)

            # Camera name
            if ihk['camera'] in hdr:
                f.write('%-6s ' % hdr[ihk['camera']])
            else:
                f.write('%-6s ' % 'None')

            # Shutter state
            if ihk['shutter'] in hdr:
                f.write('%-6s ' % hdr[ihk['shutter']])
            else:
                f.write('%-6s ' % 'None')

            # TRICK dichroic (only for OSIRIS)
            if isinstance(instrument, instruments.OSIRIS):
                f.write('%-6s ' % hdr['OBTDNAME'])
            
            # Last column is size of image 
            f.write(str(hdr["NAXIS1"]) + " X " + str(hdr["NAXIS2"]))

            # End of this line
            f.write("\n")
            

    f.close()
    
        

if __name__ == '__main__':
    _nargs = len(sys.argv)
    
    if _nargs != 2:
        print( 'Usage: kailog directory' )
    else:
        kailog(sys.argv[1])
        

def getAotsxy(hdr):
    # Note: 04jullgs does not have LSPROP or AOTSX keywords
    if (hdr.get('OUTDIR').strip() != '/sdata904/nirc2eng/2004jul26/'):
        if (hdr.get('LSPROP') == 'yes'):
            # LGS MODE
            aotsxy = [float(hdr['AOTSX']), float(hdr['AOTSY'])]

            # Another special case: 09may01 UT had a AOTSX/Y linear drift.
            # Found drift with 09maylgs/clean/kp/coord.py
            if ((hdr['OUTDIR']).strip() == '/sdata903/nirc3/2009may01/'):
                mjdobs = float(hdr['MJD-OBS'])
                # old distortion solution (pre-ship) gives:
                #aotsxy[0] -= (1.895701e5*0.727) + (-3.449709*0.727) * mjdobs 
                # new distortion solution (yelda et al. 2010)
                aotsxy[0] -= (1.903193e5*0.727) + (-3.463342*0.727) * mjdobs 

                print(( 'getAotsxy: modified from %8.3f %8.3f to %8.3f %8.3f' % \
                    (float(hdr['AOTSX']), float(hdr['AOTSY']), 
                     aotsxy[0], aotsxy[1]) ))

            # Another special case: 10jul06 UT had a AOTSX/Y linear drift.
            # Found drift with 10jullgs1/clean/kp/coord.py
            if ((hdr['OUTDIR']).strip() == '/sdata903/nirc3/2010jul06/'):
                mjdobs = float(hdr['MJD-OBS'])
                # new distortion solution (yelda et al. 2010)
                aotsxy[0] -= (2.039106e5*0.727) + (-3.681807*0.727) * mjdobs 

                print(( 'getAotsxy: modified from %8.3f %8.3f to %8.3f %8.3f' % \
                    (float(hdr['AOTSX']), float(hdr['AOTSY']), 
                     aotsxy[0], aotsxy[1]) ))
        else:
            # NGS MODE
            # Note: OBFMYIM refers to X, and vice versa!
            aotsxy = [float(hdr['OBFMYIM']), float(hdr['OBFMXIM'])]
    else:
        # 04jullgs
        # Assumes that 04jullgs has AOTSX/AOTSY keywords were added
        # by hand (see raw/fix_headers.py)
        # Note: OBFMYIM refers to X, and vice versa!
        #aotsxy = [float(hdr['OBFMYIM']), float(hdr['OBFMXIM'])]
        aotsxy = [float(hdr['AOTSX']), float(hdr['AOTSY'])]

    return aotsxy


def pix2radec():
    print( 'Not done yet' )
    return

def radec2pix(radec, phi, scale, posRef):
    """Determine pixel shifts from true RA and Dec positions.

    @param radec: a 2-element list containing the RA and Dec in degrees.
    @type radec: float list
    @param phi: position angle (E of N) in degrees.
    @type phi: float
    @param scale: arcsec per pixel.
    @type scale: float
    @param posRef: 2-element list containing the ra, dec positions (in degrees)
            of a reference object.
    @type posRef: float list
    """
    # Expected in degrees
    ra = radec[0]
    dec = radec[1]
    
    # Difference in RA and Dec. Converted to arcsec.
    d_ra = math.radians(ra - posRef[0]) * 206265.0
    d_dec = math.radians(dec - posRef[1]) * 206265.0

    cos = math.cos(math.radians(phi))
    sin = math.sin(math.radians(phi))
    cosdec = math.cos(math.radians(dec))

    d_x = (d_ra * cosdec * cos) + (d_dec * sin)
    d_y = (d_ra * cosdec * sin) - (d_dec * cos)
    d_x = d_x * (1.0/scale)
    d_y = d_y * (1.0/scale)
    
    return [d_x, d_y]

def aotsxy2pix(aotsxy, scale, aotsxyRef, inst_angle=0.0):
    # Determine pixel shifts from AOTSX and AOTSY positions.
    
    x = aotsxy[0]
    y = aotsxy[1]
    x_ref = aotsxyRef[0]
    y_ref = aotsxyRef[1]
    
    # AOTSX,Y are in units of mm. Conversion is 0.727 mm/arcsec
    d_x = (x - x_ref) / 0.727
    d_y = (y - y_ref) / 0.727
    d_x = d_x * (1.0/scale)
    d_y = d_y * (1.0/scale)
    
    # Rotate to the instrument PA
    cosa = np.cos(np.radians(-inst_angle))
    sina = np.sin(np.radians(-inst_angle))
    
    rot_matrix = np.array([[cosa, sina], [-sina, cosa]])
    coo_ao = np.array([d_x, d_y])
    coo_inst = rot_matrix.dot(coo_ao)
    
    d_x = coo_inst[0]
    d_y = coo_inst[1]
    
    return [d_x, d_y]

def pcuxy2pix(pcuxy, phi, scale, pcuxyRef):
    """Determine pixel shifts from true RA and Dec positions.

    @param pcuxy: a 2-element list containing the PCU x and y in mm.
    @type pcuxy: float list
    @param phi: position angle (E of N) in degrees.
    @type phi: float
    @param scale: mm per pixel.
    @type scale: float
    @param pcuxyRef: 2-element list containing the  PCU x, y positions (in mm) in the reference position.
    @type pcuxyRef: float list
    """
    
    # Since our reference point is the centre of the pinhole mask, which is the rotation axis, the position angle has no effect, there is just translation.
    d_x = (pcuxy[0]-pcuxyRef[0])/scale #positive x mm = positive pixels
    d_y = -(pcuxy[1]-pcuxyRef[1])/scale #positive y mm = negative pixels (after OSIRIS flip)

    return [d_x, d_y]

def pix2xyarcsec(xypix, phi, scale, sgra):
    """Determine  E and N offsets from Sgr A* (in arcsec) from 
    pixel positions and the pixel position of Sgr A*.

    xypix: 2-element list containing the RA and Dec in degrees.
    phi: position angle (E of N) in degrees.
    scale: arcsec per pixel.
    sgra: 2-element list containing the pixel position of Sgr A*.
    """
    # Expected in arcseconds
    xpix = xypix[0] - sgra[0]
    ypix = xypix[1] - sgra[1]

    cos = math.cos(math.radians(phi))
    sin = math.sin(math.radians(phi))
    
    d_x = (xpix * cos) + (xpix * sin)
    d_y = -(xpix * sin) + (ypix * cos)
    d_x = d_x * -scale
    d_y = d_y * scale
    
    return [d_x, d_y]

def xyarcsec2pix(xyarcsec, phi, scale):
    """Determine pixel shifts from E and N offsets from Sgr A*.

    xyarcsec: 2-element list containing the RA and Dec in degrees.
    phi: position angle (E of N) in degrees.
    scale: arcsec per pixel.
    """
    # Expected in arcseconds
    xarc = xyarcsec[0]
    yarc = xyarcsec[1]

    cos = math.cos(math.radians(phi))
    sin = math.sin(math.radians(phi))

    d_x = (-xarc * cos) + (yarc * sin)
    d_y = (xarc * sin) + (yarc * cos)
    d_x = d_x * (1.0/scale)
    d_y = d_y * (1.0/scale)
    
    return [d_x, d_y]

def rotate_coo(x, y, phi, instrument=instruments.default_inst):
    """Rotate the coordinates in the *.coo files for data sets
    containing images at different PAs.
    """
    # Rotate around center of image, and keep origin at center
    if instrument.name == 'NIRC2':
        xin = 512.
        yin = 512.
        xout = 512.
        yout = 512.
        
    elif instrument.name == 'OSIRIS':
        xin = 1024.
        yin = 1024.
        xout = 1024.
        yout = 1024.

    else:
        raise Exception('Inst not supported. Edit kai_util.rotate_coo to specify center of new instrument')
  
    cos = math.cos(math.radians(phi))
    sin = math.sin(math.radians(phi))

    xrot = (x - xin) * cos - (y - yin) * sin + xout
    yrot = (x - xin) * sin + (y - yin) * cos + yout
   
    return [xrot, yrot]
    

def getScale(hdr, instrument=instruments.default_inst):
    return instrument.get_plate_scale(hdr)

def getPA(hdr, instrument=instruments.default_inst):
    return instrument.get_position_angle(hdr)

def getCentralWavelength(hdr, instrument=instruments.default_inst):
    return instrument.get_central_wavelength(hdr)

def calcOverhead(tint, coadds, ndithers, nframes, reads, tread=0.181):
    t = 0.0
    
    if (ndithers > 1):
        t += 6.0 * (ndithers + 1.0)
        
    t += 12.0 * nframes * ndithers
    t += ndithers * nframes * coadds * (tint + tread*(reads - 1.0))

    tmin = t / 60.0

    print(( '\tIntegration time: %.3f' % tint ))
    print(( '\tNumber of Coadds: %d' % coadds ))
    print(( '\tNumber of Dither Positions: %d' % ndithers ))
    print(( '\tFrames at each Position: %d' % nframes ))
    print(( '\tNumber of Reads: %d' % reads ))
    print(( 'Total elapsed observing time = %5d sec (or %5.1f min)' % \
          (t, tmin) ))

def plotKeyword(keyword1, keyword2, imgList):
    """
    Pass in a file containing a list of images. For each of these
    images, read out the values of the header keywords specified.
    Then plot each of the keywords against each other.
    """
    tab = Table.read(imgList, format='ascii', header_start=None)

    files = [tab[i][0].strip() for i in range(len(tab))]

    value1 = np.zeros(len(tab), dtype=float)
    value2 = np.zeros(len(tab), dtype=float)

    print(( keyword1, keyword2 ))

    for ff in range(len(files)):
        hdr = fits.getheader(files[ff], ignore_missing_end=True)

        value1[ff] = hdr[keyword1]
        value2[ff] = hdr[keyword2]


    import pylab as py
    py.clf()

    py.plot(value1, value2, 'k.')
    py.xlabel(keyword1)
    py.ylabel(keyword2)

    return (value1, value2)
    
    
def restackCoadd(clean_directory, raw_directory, raw_unp_directory, save_diagnostics=False, save_coadds=False):
    """
    Process unprocessed FITS frames and coadds to align star positions, stack them,
    and optionally save diagnostics and intermediate outputs. This stacks to the first coadd.

    @param clean_directory: Path to directory containing .coo files with reference star positions.
    @type clean_directory: str
    @param raw_directory: Path to directory containing raw FITS frames.
    @type raw_directory: str
    @param raw_unp_directory: Path to directory containing unprocessed FITS files (n_unp_###.fits).
    @type raw_unp_directory: str
    @param save_diagnostics: If True, saves intermediate plots and motion diagnostics.
    @type save_diagnostics: bool
    @param save_coadds: If True, saves individual coadd FITS frames.
    @type save_coadds: bool
    """
    def recenter_star(data, x_init, y_init, box=12):
        """
        Recenter a star near (x_init, y_init) using centroiding on a small cutout.

        @param image_path: Path to the FITS file to measure.
        @type image_path: str
        @param x_init: Initial x estimate for the star.
        @type x_init: float
        @param y_init: Initial y estimate for the star.
        @type y_init: float
        @param box: Half-size of the cutout box in pixels (e.g., 12 = 24x24 window).
        @type box: int
        @return: Refined (x, y) coordinates.
        @rtype: tuple(float, float)
        """
        try:
            # Make a cutout around the initial guess
            position = (x_init, y_init)
            cutout = Cutout2D(data, position=position, size=2*box)

            # Compute centroid in the cutout frame
            cx, cy = centroid_com(cutout.data)

            # Convert back to full-image coordinates
            x_refined = cx + cutout.origin_original[0]
            y_refined = cy + cutout.origin_original[1]

            return float(x_refined), float(y_refined)

        except Exception as e:
            print(f"Centroiding failed: {e}")
            return x_init, y_init
            
    # Create log file for tracking offset measurements of IRS16C
    f = open(raw_unp_directory + '/restack.log', 'w')
    f.write('%16s  ' % '#frame')
    f.write('%16s  ' % 'frame_x')
    f.write('%16s  ' % 'frame_y')
    f.write('%16s  ' % 'coadd')
    f.write('%16s  ' % 'coadd_x')
    f.write('%16s  ' % 'coadd_y')
    f.write('%16s  ' % 'dx')
    f.write('%16s  ' % 'dy')
    f.write('\n')
    
    stack_mode = 'per_coadd'
    unp_loc = raw_unp_directory
    
    # Organize raw frame indices
    unp_files = sorted(glob.glob('{}/n_unp_*.fits'.format(unp_loc)))
    wildcards = sorted([os.path.basename(f)[6:-5] for f in unp_files])

    # Create folders to save diagnostics in for each frame
    if save_diagnostics:
        if not os.path.exists("{}/restack_diagnostics/".format(unp_loc)):
            os.makedirs("{}/restack_diagnostics/".format(unp_loc))
    
        for w in wildcards:
            dir_path = '{}/restack_diagnostics/n{}/'.format(unp_loc, w)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

    # Loop through the unprocessed fits files
    for cur_file, w in zip(unp_files, wildcards):
        # Instrument
        instrument=instruments.default_inst
        
        # Get reference raw header from raw_directory
        ref_header = fits.getheader('{}/n{}.fits'.format(raw_directory, w))
        filt = instrument.get_filter_name(ref_header).lower()
        

        # Get initial reference position from clean_directory .coo file
        ref_coo = '{}/{}/c{}.coo'.format(clean_directory, filt, w)

        if not os.path.exists(ref_coo):
            raise FileNotFoundError("{}: .coo file not found. Make sure this exists, or create one manually.".format(ref_coo))
            continue

        ref_x, ref_y = np.loadtxt(ref_coo)
        
        # Build difference reads for each coadd
        with fits.open(cur_file) as hdul:
            num_coadds = len(hdul) - 1
            diff_frames = []
        
            for coadd_index in (range(1, num_coadds + 1)):
                read_data = hdul[coadd_index].data 
                num_reads = read_data.shape[0]
                header = hdul[0].header

                if stack_mode == 'per_coadd':
                    n_half = num_reads // 2
                    avg_initial = np.mean(read_data[:n_half], axis=0)
                    avg_final = np.mean(read_data[-n_half:], axis=0)
                    frame = avg_final - avg_initial
                    diff_frames.append(frame)

        # Align and stack coadds
        stack = None
        dx_positions = []
        dy_positions = []
        x_positions = []
        y_positions = []
        for idx, frame in enumerate(diff_frames):



            # Get updated coordinate position of IRS16C for current coadd
            x_new, y_new = recenter_star(frame, ref_x, ref_y, box=12)
            
            # Log motion diagnostics
            f.write('%16s  ' % 'n{}'.format(w))
            f.write('%16s  ' % ref_x)
            f.write('%16s  ' % ref_y)
            f.write('%16s  ' % '{}'.format(idx))
            f.write('%16s  ' % round(x_new, 2))
            f.write('%16s  ' % round(y_new, 2))
            
            dx = ref_x - x_new
            dy = ref_y - y_new
            
            f.write('%16s  ' % round(dx, 2))
            f.write('%16s  ' % round(dy, 2))
            f.write('\n')
            
            dx_positions.append(ref_x - x_new)
            dy_positions.append(ref_y - y_new)
            x_positions.append(x_new)
            y_positions.append(y_new)
            
            # Shift and stack coadd
            shifted = shift(frame, shift=(dy, dx), order=3)
            stack = shifted if stack is None else stack + shifted
        
            # Save individual coadd if flag is on
            if save_coadds:
                header = ref_header
                header["COADD"] = idx
                hdu = fits.PrimaryHDU(data=frame, header=header)
                hdu.writeto('{}/n{}_{}.fits'.format(unp_loc, w, int(idx)), overwrite=True)
                print("Saved FRAME: n{} -- COADD: {}".format(w, int(idx)))
                
            # Plot motion diagnostics if flag is on
            if save_diagnostics:
                box_size = 200
                half_box = box_size // 2

                fig, ax = plt.subplots(figsize=(8, 8))
                ax.imshow(frame, origin='lower', cmap='gray',
                          vmin=np.percentile(frame, 5), vmax=np.percentile(frame, 99))

                ax.set_xlim(ref_x - half_box, ref_x + half_box)
                ax.set_ylim(ref_y - half_box, ref_y + half_box)

                ax.add_patch(Circle((ref_x, ref_y), radius=5, color='blue', fill=False, label='Frame .coo'))
                ax.add_patch(Circle((x_new, y_new), radius=5, color='red', fill=False, label='Coadd .coo'))
                ax.legend()
                plt.title("COADD: {}".format(int(idx)))
                plt.tight_layout()
                plt.savefig('{}/restack_diagnostics/n{}/coadd{:04d}.png'.format(unp_loc, w, int(idx)))
                plt.close()
            
        # Plot shifts in x and y if flag is on
        if save_diagnostics:
            x_offsets = np.subtract(x_positions, x_positions[0])
            y_offsets = np.subtract(y_positions, y_positions[0])
            all_offsets = np.concatenate([x_offsets, y_offsets])
            ymin, ymax = all_offsets.min(), all_offsets.max()
            padding = 0.05 * (ymax - ymin)
            ymin -= padding
            ymax += padding
            # Plot Δx and Δy over time (frame index)
            plt.figure(figsize=(10, 5))
            plt.subplot(1, 2, 1)
            plt.plot(x_offsets, label='x', marker='o', color = 'blue')
            plt.xlabel('Coadd number')
            plt.ylabel('x Offset [pixels]')
            plt.axhline(0, linestyle = '--', color = 'red')
            plt.title('x over coadds')
            plt.ylim(ymin, ymax)
            plt.grid(True)

            # Quiver plot: motion vectors
            plt.subplot(1, 2, 2)
            plt.plot(y_offsets, label='y', marker='o', color = 'green')
            plt.xlabel('Coadd number')
            plt.ylabel('y Offset [pixels]')
            plt.title('y over coadds')
            plt.axhline(0, linestyle = '--', color = 'red')
            plt.ylim(ymin, ymax)
            plt.grid(True)

            plt.tight_layout()
            plt.savefig('{}/restack_diagnostics/n{}/movement_per_coadd.png'.format(unp_loc, w))
            plt.close()

            print("Saved motion diagnostics plot")

        # Final re-stacked image output
        ref_header['STACKED'] = True
        ref_header['NFRAMES'] = len(diff_frames)
        hdu = fits.PrimaryHDU(data=stack, header=ref_header)

        # Save to subdir
        final_fits_path = '{}/n{}.fits'.format(unp_loc, w)
        hdu.writeto(final_fits_path, overwrite=True)
        print("Re-stacked frame: n{}".format(w))

        # Make GIF of coadd shifting if flag is on
        if save_diagnostics:
            # Make GIF from coadd images
            coadd_pngs = sorted(glob.glob('{}/restack_diagnostics/n{}/coadd*.png'.format(unp_loc, w)))
            gif_path = '{}/restack_diagnostics/n{}/n{}_coadds.gif'.format(unp_loc, w, w)

            images = []
            for png in coadd_pngs:
                try:
                    images.append(imageio.imread(png))
                except Exception as e:
                    print("Skipping {} — {}".format(png, e))

            if images:
                imageio.mimsave(gif_path, images, duration=0.3)
                print("GIF saved to {}".format(gif_path))
            else:
                print("No valid PNGs found for {}".format(w))

        
    f.close()

def add_wcs_keywords(header, instrument):
    """
    Adds WCS keywords to file. 
    If WCS exists, this overwrites it.

    Parameters
    ----------
    header: fits header
        Fits header of file to be updated

    instrument: instrument object
        Instrument obejct from kai.instruments
    """
    ra, dec = instrument.get_radec(header)
    rotposn = instrument.get_position_angle(header)
    pixel_scale = instrument.get_plate_scale(header)
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']

    crpix1 = naxis1/2 + 0.5
    crpix2 = naxis2/2 + 0.5

    #pixel scale to degrees
    cdelt = pixel_scale/3600

    # Calculate CD matrix for rotation
    rot_rad = np.radians(rotposn)
    cos_rot = np.cos(rot_rad)
    sin_rot = np.sin(rot_rad)
    
    cd1_1 = -cdelt * cos_rot  # Negative for standard orientation
    cd1_2 = cdelt * sin_rot
    cd2_1 = cdelt * sin_rot
    cd2_2 = cdelt * cos_rot
    
    # Add/update WCS keywords
    wcs_keywords = {
        # Coordinate system
        'CTYPE1': 'RA---TAN',
        'CTYPE2': 'DEC--TAN',
        
        # Reference coordinates
        'CRVAL1': ra,
        'CRVAL2': dec,
        
        # Reference pixels
        'CRPIX1': crpix1,
        'CRPIX2': crpix2,
        
        # CD matrix (includes scale and rotation)
        'CD1_1': cd1_1,
        'CD1_2': cd1_2,
        'CD2_1': cd2_1,
        'CD2_2': cd2_2,
        
        # Coordinate system
        'RADESYS': 'ICRS',
        'EQUINOX': 2000.0,
        
        # WCS comments
        'CUNIT1': 'deg',
        'CUNIT2': 'deg'
    }

    print('Adding WCS to ', header['DATAFILE'])
    for key, value in wcs_keywords.items():
        header[key] = value
        print(f"  {key} = {value}")
    
    # Add comments for clarity
    header.comments['CTYPE1'] = 'Coordinate type for axis 1'
    header.comments['CTYPE2'] = 'Coordinate type for axis 2'
    header.comments['CRVAL1'] = 'Reference RA in decimal degrees'
    header.comments['CRVAL2'] = 'Reference Dec in decimal degrees'
    header.comments['CRPIX1'] = 'Reference pixel for axis 1'
    header.comments['CRPIX2'] = 'Reference pixel for axis 2'
    header.comments['CD1_1'] = 'Transformation matrix element'
    header.comments['CD1_2'] = 'Transformation matrix element'
    header.comments['CD2_1'] = 'Transformation matrix element'
    header.comments['CD2_2'] = 'Transformation matrix element'

    return header
