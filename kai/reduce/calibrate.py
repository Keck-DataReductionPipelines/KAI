#!/usr/bin/env python
import optparse
import textwrap
import numpy as np
import pylab as py
from astropy.table import Table
import math
import sys
from astropy.io import fits as pyfits
import pdb

# Map of the possible plate scales
all_scales = [[1.0, 'No scaling'],
              [0.0102, 'Speckle'],
              [0.0087, 'KCAM-AO'],
              [0.0085, 'SCAM-AO'],
              [0.00993, 'NIRC2-AO narrow'],
              [0.02000, 'NIRC2-AO medium'],
              [0.04000, 'NIRC2-AO wide'],
              [0.0170, 'SCAM-AO unmagnified'],
              [0.030, 'NIRC2-AO narrow, binned by 3'],
              [0.01998, 'GEMINI'],
              [0.107, 'MMT PISCES'],
              [0.200, 'UKIDSS'],
              [0.004, 'TMT/IRIS'],
              [0.0196, 'GSAOI'],
              [0.050, 'ACS-WFC'],
              [0.01,  'OSIRIS-imag']]

##################################################
# 
# Help formatter for command line arguments. 
# This is very generic... skip over for main code.
#
##################################################
class IndentedHelpFormatterWithNL(optparse.IndentedHelpFormatter):
    def format_description(self, description):
        if not description: return ""
        desc_width = self.width - self.current_indent
        indent = " "*self.current_indent
        # the above is still the same
        bits = description.split('\n')
        formatted_bits = [
            textwrap.fill(bit,
                          desc_width,
                          initial_indent=indent,
                          subsequent_indent=indent)
            for bit in bits]
        result = "\n".join(formatted_bits) + "\n"
        return result

    def format_option(self, option):
        # The help for each option consists of two parts:
        #   * the opt strings and metavars
        #   eg. ("-x", or "-fFILENAME, --file=FILENAME")
        #   * the user-supplied help string
        #   eg. ("turn on expert mode", "read data from FILENAME")
        #
        # If possible, we write both of these on the same line:
        #   -x    turn on expert mode
        #
        # But if the opt string list is too long, we put the help
        # string on a second line, indented to the same column it would
        # start in if it fit on the first line.
        #   -fFILENAME, --file=FILENAME
        #       read data from FILENAME
        result = []
        opts = self.option_strings[option]
        opt_width = self.help_position - self.current_indent - 2
        if len(opts) > opt_width:
            opts = "%*s%s\n" % (self.current_indent, "", opts)
            indent_first = self.help_position
        else: # start help on same line as opts
            opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
            indent_first = 0
        result.append(opts)
        if option.help:
            help_text = self.expand_default(option)
            # Everything is the same up through here
            help_lines = []
            for para in help_text.split("\n"):
                help_lines.extend(textwrap.wrap(para, self.help_width))
            # Everything is the same after here
            result.append("%*s%s\n" % (
                    indent_first, "", help_lines[0]))
            result.extend(["%*s%s\n" % (self.help_position, "", line)
                           for line in help_lines[1:]])
        elif opts[-1] != "\n":
            result.append("\n")
        return "".join(result)



##################################################
#
# Main body of calibrate code.
#
##################################################
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # Read options and check for errors.
    options = read_command_line(argv)
    if (options == None):
        return
    
    # Read in the photometric calibrators
    calibs = read_photo_calib_file(options)
    
    # Read in the starlist
    stars = input_data(options)
    
    # Match up calibrator stars with stars in our starlist.
    calibs.index = find_cal_stars(calibs, stars, options)

    # Calculate the average zeropoint (and errors)
    zeropt, zeropt_err = calc_zeropt(calibs, stars, options)

    # Write out the output files
    output_new(zeropt, zeropt_err, calibs, stars, options)
  
    
def read_command_line(argv):
    p = optparse.OptionParser(usage='usage: %prog [options] [starlist]',
                              formatter=IndentedHelpFormatterWithNL())

    p.add_option('-f', dest='data_type', type=int, default=1, metavar='[#]',
                 help='Data type of the input starlist:\n'+
                 '[1] starfinder format without errors (default)\n'+
                 '[2] align_rms format with errors\n')

    scaleHelp = 'Camera type to set plate scale (default: %default):\n'
    for ss in range(len(all_scales)):
        scaleHelp += '[%d] %s (%7.5f asec/pix) \n' % \
            (ss, all_scales[ss][1], all_scales[ss][0])
        
        
    p.add_option('-c', dest='camera_type', type=int, default=1, metavar='[#]',
                 help=scaleHelp)
    p.add_option('-T', dest='theta', type=float, default=0.0, metavar='[ANGLE]',
                 help='Rotation of the image from North with positive angles '+
                 'in the clockwise direction (default: %default)')
    p.add_option('-I', dest='first_star', default='irs16C', metavar='[STAR]',
                 help='Name of the first star in the list. Any star from the '+
                 'input calibration list (see -N) can be used. The coo star'+
                 'should have a known magnitude (within 0.5 mag) in the '+
                 'calibration list in order to properly find the rest of the '
                 'calibrator stars. The coo star does not have to be used '+ 
                 'as a calibrator. (default: %default)')
    p.add_option('-N', dest='calib_file', metavar='[FILE]',
                 default='/u/ghezgroup/data/gc/source_list/photo_calib.dat',
                 help='A file containing the calibration stars and magnitudes '+
                 '(default: %default). The file format should have different '+
                 'columns for each photometric calibration reference. The '+
                 'default calibration sources for each column are stored in '+
                 'a comment at the head of the file. Choose which column to '+
                 'use with the -M flag.')
    p.add_option('-M', dest='calib_column', type=str, default='Kp', metavar='[FILTER]',
                 help='Choose which column to use in the photo_calib.dat '+
                 'file (default: %default). See below for column choices.')
    p.add_option('-S', dest='calib_stars', metavar='[STAR1,STAR2]', default=None,
                 help='Specify the stars to be used as calibrators in a '+
                 'comma-separated list. The defaults are set in the '+
                 'photo_calib.dat file at the top (or see below).')
    p.add_option('-A', dest='align_stars', metavar='[STAR1,STAR2]', default=None,
                 help='Specify the stars to be used as align stars in a '+
                 'comma-separated list. Used when reordering output starlist '+
                 '(flag: -R). If not specified, the defaults are set '+
                 'as the calibrator stars (flag: -S).')
    p.add_option('-r', dest='outroot', metavar='[ROOT]',
                 help='Rename root name for output (default: [listname]_cal)')
    p.add_option('-V', dest='verbose', action='store_true', default=False,
                 help='Print extra diagnostics.')
    p.add_option('-R', dest='reorder', action='store_true', default=False,
                 help='Reorder the starlist so that the align stars '+
                 'are at the top. The coo star (first star) will remain first.')
    p.add_option('-s', '--snr', dest='snr_flag', default=0, metavar='[#]', type='int',
                 help='Use this flag to indicate the purpose of the SNR '+
                 'column in the input star list (default: %default). '+
                 'The choices are:\n'+
                 '[0] overwrite SNR column with calibration error (default)\n'+
                 '[1] add SNR in quadrature to calibration error\n'+
                 '[2] leave the SNR column alone')
    p.add_option('--searchRadius', dest='searchRadius', default=0.20, metavar='[#]', type='float',
                 help='Set the search radius (in arcsec) for matching the star' +
                 'in the starlist to the calibration star (default: %default arcsec).' +
                 'This is the search radius for the bright stars. Faint stars (below ' +
                 'a magnitude set by the --brightLimit flag) have a search radius that ' +
                 'is two times smaller (searchRadius / 2).')
    p.add_option('--searchMag', dest='searchMag', default=1.0, metavar='[#]', type='float',
                 help='Set the delta-magnitude for matching the star in the ' +
                 'starlist to the calibration star (default: %default).')
    p.add_option('--brightLimit', dest='brightLimit', default=12, metavar='[#]', type='float',
                 help='Set the brightness limit above which stars are matched within ' +
                 'a search radius set by --searchRadius and below which stars are ' +
                 'matched within searchRadius / 2. (default: %default).')
    
    options, args = p.parse_args(argv)
    options.searchMag = np.abs(float(options.searchMag))
    options.searchRadius = np.abs(float(options.searchRadius))
    options.brightLimit = float(options.brightLimit)    
    
    # Keep a copy of the original calling parameters
    options.originalCall = ' '.join(argv)

    # Read the input filename
    options.input_file = None
    if len(args) == 1:
        options.input_file = args[0]
    else:
        print( '' )
        p.print_help()
        read_photo_calib_file(options, verbose=True)
        return None

    # Set the output filenames
    if (options.outroot == '' or options.outroot == None):
        parts = options.input_file.split('.')
    else:
        parts = options.outroot.split('.')
    if (len(parts) > 1):
        root = '.'.join(parts[0:-1])
        extension = parts[-1]
    else:
        root = parts[0]
        extension = 'lis'
    options.outname = root + '_cal.' + extension
    options.zername = root + '_cal.zer'

    # Set plate scale
    options.plate_scale = all_scales[options.camera_type][0]

    # Parse calib and align stars 
    if (options.align_stars != None):
        options.align_stars = options.align_stars.split(',')
    elif (options.calib_stars != None):
        ## If align stars not specified, use the specified calib stars
        options.align_stars = options.calib_stars.split(',')
        
    if (options.calib_stars != None):
        options.calib_stars = options.calib_stars.split(',')

    # Verbose mode printing
    if options.verbose:
        print( 'VERBOSE mode on' )
        print(( 'options.first_star = %s' % options.first_star ))
        print(( 'options.data_type = %d' % options.data_type ))
        print(( 'options.camera_type = %d' % options.camera_type ))
        print(( 'options.plate_scale = %7.2f' % options.plate_scale ))
        print(( 'options.outname = %s' % options.outname ))
        print(( 'options.calib_file = %s' % options.calib_file ))
        print(( 'options.calib_column = %s' % options.calib_column ))
        print(( 'options.theta = %6.1f' % options.theta ))
        print(( 'options.searchRadius = %6.1f arcsec' % options.searchRadius ))
        print(( 'options.searchMag = %6.1f' % options.searchMag ))
        print(( 'options.brightLimit = %6.1f' % options.brightLimit ))
        if options.reorder:
            print( 'Reordering lis file' )
        else:
            print( 'Not reordering lis file.' )

    return options


def read_photo_calib_file(options, verbose=False):
    # Read in filter information from photo_calib file header
    f_calib = open(options.calib_file, 'r')
    
    filter_names = []
    filter_info = {}
    filter_defaultStars = {}
    
    if verbose or options.verbose:
        print( '' )
        print( 'Photometric calibration information loaded from:' )
        print(( '\t', options.calib_file ))
        print( 'Specify a different file with the -N flag.' )
        print( 'Choose a calibration column with the -M flag.' )
        print( 'The column choices are listed by [FILTER] below.' )
        print( '' )
        print( 'Bandpass and References:' )
    
    for line in f_calib:
        if line.startswith('# Filt'):
            ## Filter header lines. Read out filter names, filter info, and default calibrator stars.
            fields = line.split('--')
            
            filter_name = (fields[1].split(': '))[0].strip()
            
            filter_names.append(filter_name)
            filter_info[filter_name] = (fields[1].split(': '))[1]
            filter_defaultStars[filter_name] = fields[2]
            
            if verbose or options.verbose:
                print('{0}\t{1}\t{2}'.format(filter_name, filter_info[filter_name], filter_defaultStars[filter_name]))
        elif line.startswith('#'):
            ## Just an otherwise ordinary comment line, skip.
            continue
        else:
            ## First line after the header, so done reading header. Break out!
            break

    f_calib.close()
    
    if verbose or options.verbose:
        print( '' )
        print( 'Calibration Sources:' )
        print( '\t(* default if no -S flag)' )
        print( '\t(! stars determined to be variable)' )
        print( '' )
    
    ##########
    #
    # Read in the data portion of the file.
    #
    ##########
    #pdb.set_trace()
    tab = Table.read(options.calib_file, format='ascii.commented_header', delimiter='\s', header_start=-1)
    
    name_col = tab['Star']
    x_col = tab['x_pos']
    y_col = tab['y_pos']
    x_vel_col = tab['x_vel']
    y_vel_col = tab['y_vel']
    t0_col = tab['t0']
    isVar_col = (tab['var?'] == 1)
    
    
    ## Construct mag matrix and default calibrators for each filter.
    ## Used for printing in verbose mode
    
    magMatrix = np.zeros((len(filter_names), len(tab)), dtype=float)
    isDefaultMatrix = np.zeros((len(filter_names), len(tab)), dtype=bool)
    
    for i in range(len(filter_names)):
        magMatrix[i,:] = tab[filter_names[i]]
        
        # If no default stars were set, then assume
        # all stars with non-zero magnitudes are the defaults.
        if filter_defaultStars[filter_names[i]] == None:
            idx = np.where(magMatrix[i,:] != 0)[0]
            isDefaultMatrix[i,idx] = True

        else:
            stars = filter_defaultStars[filter_names[i]].split(',')
            stars = [stars[s].strip() for s in range(len(stars))]
            for s in range(len(tab)):
                if (name_col[s] in stars):
                    isDefaultMatrix[i,s] = True
                else:
                    isDefaultMatrix[i,s] = False
    
    ##########
    # Print out 
    ##########
    if verbose or options.verbose:
        print_line = ' %10s ' % 'Name'
        for i in range(len(filter_names)):
            print_line += ' {0} '.format(filter_names[i])

        print_line += '\n'
        print( print_line )

        for s in range(len(tab)):
            varChar = '!' if isVar_col[s] else ''
            print_line = '%1s%13s ' % (varChar, name_col[s])
            
            for i in range(len(filter_names)):
                defChar = '*' if isDefaultMatrix[i,s] else ''
                print_line += ' %5.2f%1s ' % (magMatrix[i,s], defChar)

            print_line += ''

            print(print_line)


    ##########
    # Catch case where this is called
    # just to print out the file (no starlist to calibrate)
    ##########
    if options.input_file == None:
        return None
    
    
    # Check if the filter is in photo_calib file
    if not (options.calib_column in filter_names):
        print('Filter {0} not in {1}'.format(options.calib_column, options.calib_file))
        print("Available filter options are: '{0}'".format("', '".join(filter_names)))
        
        return None
    
    
    ##########
    # Get only the calibration magnitudes
    # that were asked for.
    ##########
    calibs = Starlist()
    calibs.name = name_col
    calibs.x = x_col
    calibs.y = y_col
    calibs.xVel = x_vel_col
    calibs.yVel = y_vel_col
    calibs.t0 = t0_col
    
    # Pick out the specific filter column
    calibs.mag = tab[options.calib_column]
    calibs.mag_err = tab[options.calib_column + '_err']
    calibs.filt_info = filter_info[options.calib_column]
    calibs.filt_name = options.calib_column

    ##########
    #
    # Decide which stars will be included as
    # calibrators based on user input. This 
    # sets up calibs.include which is set to
    # True for those calibrators that should be
    # used in the photometric calibration.
    #
    ##########
    
    if (options.calib_stars == None):
        filter_index = filter_names.index(options.calib_column)
        
        calibs.include = isDefaultMatrix[filter_index, :]
    else:
        calib_stars_search = options.calib_stars
        
        calibs.include = np.zeros(len(name_col), dtype=bool)

        for calib_search_index in range(len(calib_stars_search)):
            idx = np.where(name_col == calib_stars_search[calib_search_index])[0]
    
            if len(idx) == 0:
                msg = 'Failed to find calibrator: %s' % \
                    calib_stars_search[calib_search_index]
                raise Exception(msg)
        
            calibs.include[idx] = True

            if options.verbose:
                print( 'Found calibrator: ', name_col[idx], ' ', calib_stars_search[calib_search_index] )
    
    return calibs
    
def input_data(options):
    """
    Read in a starlist and return a starlist object
    with the following hanging off:
      name
      mag
      epoch
      x
      y
      xerr (none for starfinder starlists)
      yerr (none for starfinder starlists)
      snr
      corr
      nframes
      fwhm
    """
    if options.verbose:
        print(( 'Opening starlist: ', options.input_file ))

    tab = Table.read(options.input_file, format='ascii', delimiter='\s')
    cols = tab.colnames
    
    name = tab[cols[0]]
    mag = tab[cols[1]]
    epoch = tab[cols[2]]
    x = tab[cols[3]]
    y = tab[cols[4]]
    
    if options.data_type == 2:
        xerr = tab[cols[5]]
        yerr = tab[cols[6]]
        snr = tab[cols[7]]
        corr = tab[cols[8]]
        nframes = tab[cols[9]]
        fwhm = tab[cols[10]]
    else:
        xerr = None
        yerr = None
        snr = tab[cols[5]]
        corr = tab[cols[6]]
        nframes = tab[cols[7]]
        fwhm = tab[cols[8]]

        
    # Trim out stars with errors in magnitudes
    idx = np.where(mag != float('Inf'))[0]

    if len(idx) > 0:
        name = name[idx]
        mag = mag[idx]
        epoch = epoch[idx]
        x = x[idx]
        y = y[idx]
        snr = snr[idx]
        corr = corr[idx]
        nframes = nframes[idx]
        fwhm = fwhm[idx]
        
        if not (xerr is None):
            xerr = xerr[idx]
            yerr = yerr[idx]

    if options.verbose:
        print(( 'Read %d lines in the input file.' % (len(tab)) ))
        print(( 'Skipped %d lines in the input file.' % (len(tab) - len(idx)) ))

    starlist = Starlist()
    starlist.name = name
    starlist.mag = mag
    starlist.epoch = epoch
    starlist.x = x
    starlist.y = y
    starlist.xerr = xerr
    starlist.yerr = yerr
    starlist.snr = snr
    starlist.corr = corr
    starlist.nframes = nframes
    starlist.fwhm = fwhm

    return starlist
    
def find_cal_stars(calibs, stars, options):
    """
    Returns an array of indices which holds the index of the 
    matching star in the starlist. Non-matches have an index of -1. 
    """
    # First we need to find out if the first star
    # in the star list is in our list of calibrators.
    fidx = np.where(calibs.name == options.first_star)[0]
    if (len(fidx) == 0):
        msg =  'Failed to find the first star in the calibrators:\n'
        msg += '  %s' % options.first_star
        raise Exception(msg)

    # Change the size of the stars table, name column to hold the largest possible
    # calib star name.
    if (stars.name.dtype < calibs.name.dtype):
        stars.name = stars.name.astype(calibs.name.dtype)
    
    # Account for velocity of stars
    delta_t = stars.epoch[0] - calibs.t0    # Using first star's epoch as current time

    calibs.x += delta_t * (calibs.xVel/1000.)   # Velocities in mas/yr, positions in arcsec
    calibs.y += delta_t * (calibs.yVel/1000.)
        
    # Change the positional offsets to be relative to 
    # the reference source.
    calibs.x -= calibs.x[fidx]
    calibs.y -= calibs.y[fidx]
    
    # Determine the pixel positions for the calibrators
    cosScale = math.cos(math.radians(options.theta)) / options.plate_scale
    sinScale = math.sin(math.radians(options.theta)) / options.plate_scale
    calibs.xpix = stars.x[0] - (calibs.x * cosScale) + (calibs.y * sinScale)
    calibs.ypix = stars.y[0] + (calibs.x * sinScale) + (calibs.y * cosScale)
    if options.verbose:
        for c in range(len(calibs.xpix)):
            print(( 'Looking for %10s at (%.2f, %.2f), mag: %.2f' % \
                (calibs.name[c], calibs.xpix[c], calibs.ypix[c], calibs.mag[c]) ))
            
    # Create an array of indices into the starlist.
    # Set to -1 for non-matches. 
    index = np.ones(len(calibs.name), dtype=int) * -1 

    # Loop through all the calibrators and find their match in the starlist.
    # search radius, search mag, and bright limit are set by flags.
    # default search radius = 0.20 arcsec for bright sources
    # default search mag = 1.0
    # default bright limit = 12
    options.searchRadius /= options.plate_scale   
    
    magAdjust = stars.mag[0] - calibs.mag[fidx]
    if options.verbose:
        print(( 'Search dr = %d pixels, dm = %.2f' % (options.searchRadius, options.searchMag) ))
        print(( 'Adjusting input magnitudes by %.2f' % magAdjust ))

    for c in range(len(calibs.name)):
        dx = stars.x - calibs.xpix[c]
        dy = stars.y - calibs.ypix[c]
        dr = np.hypot(dx, dy)
        dm = np.abs(stars.mag - calibs.mag[c] - magAdjust)
        
        # Find the matches within our tolerance.
        if (calibs.mag[c] < options.brightLimit):
            # For the bright stars we have the default search radius:
            idx = np.where((dr < options.searchRadius) & (dm < options.searchMag))[0]
        else:
            # For the fainter stars, use a smaller search radius:
            idx = np.where((dr < options.searchRadius/2) & (dm < options.searchMag))[0]

        # Default is not found
        index[c] = -1

        # But if we find one, change names, record index, etc.
        if (len(idx) > 0):
            # Record match index into the starlist
            index[c] = idx[0]
            
            # Rename
            origName = stars.name[index[c]]
            stars.name[index[c]] = calibs.name[c]

            # Print out
            if options.verbose:
                notUsed = '' if (calibs.include[c] == True) else '(not used)'
                print(( '%10s found at %.2f, %.2f as %s %s' % \
                    (calibs.name[c], 
                     stars.x[index[c]],
                     stars.y[index[c]],
                     origName,
                     notUsed)))

        else:
            if options.verbose:
                print(( '%10s not found' % (calibs.name[c]) ))
    
    return index

def calc_zeropt(calibs, stars, options):
    """
    Calculate the average zeropoints from all the
    calibrator stars specified by the user (or defaults).
    Recall that not all calibrators in our list will
    be used for the calibration... only those with
    both include = True and found in the starlist.
    """

    # Identify the calibrators we will use to calculate
    # the zeropoint.
    cidx = np.where((calibs.include == True) & (calibs.index != -1))[0]
    sidx = calibs.index[cidx]

    # Calculate the zero point as the difference between the calculated 
    # magnitude and the published magnitude for each calibration star.
    dm = calibs.mag[cidx] - stars.mag[sidx]
    dm_err = calibs.mag_err[cidx] 
    
    ## dm in flux space (i.e. all the zeropoint adjustments)
    all_zeropts = np.power(10, -dm/2.5)
    
    ## Convert error of magnitude into flux space, and derive weights for mean
    all_zeropt_errs = (all_zeropts * dm_err)/1.0857
    
    ### Default case: if any of the calibrators have 0 error in reference mag, use weights = 1 for all calibrators
    all_zeropt_weights = (0. * all_zeropt_errs) + 1.
    
    ### If all calibrators have an error in reference mag, calculate weights based on mag errors
    if np.all(all_zeropt_errs != 0.0):
        all_zeropt_weights = 1./(all_zeropt_errs ** 2.)
    
    
    # Take mean as the value and standard deviation as the error.
    # Note, we are not using the error on the mean.
    
    ## Weighted mean
    zeropt = np.sum(all_zeropts * all_zeropt_weights) / np.sum(all_zeropt_weights)

    ## Weighted standard deviation
    zeropt_err = np.sqrt(np.sum(all_zeropt_weights * ((all_zeropts - zeropt)**2.)) / np.sum(all_zeropt_weights))
    
    ## Previous calculations (without weights):
    # zeropt = all_zeropts.mean()
    # zeropt_err = all_zeropts.std(ddof=1)

    # Using the St.Dev propagate the errors into mag space
    zeropt_err *= 1.0857 / zeropt

    # Convert the flux ratio back to a magnitude difference
    zeropt = -2.5 * np.log10( zeropt )

    if options.verbose:
        print( '' )
        print(( 'Zero-point = %5.3f +/- %.3f' % (zeropt, zeropt_err) ))
        print( '' )
        for i in range(len(cidx)):
            c = cidx[i]
            s = sidx[i]
            print(( '%10s Published Mag: %6.3f  Calculate Mag: %6.3f  DIFF = %6.3f' % \
                (calibs.name[c], calibs.mag[c], stars.mag[s]+zeropt,
                 calibs.mag[c] - (stars.mag[s]+zeropt)) ))

    return (zeropt, zeropt_err)


def output_new(zeropt, zeropt_err, calibs, stars, options):
    """
    Write out a calibrated starlist and a *.zer file with 
    the calculated zeropoints.
    """
    # Update the magnitudes of all the stars
    stars.mag += zeropt

    # Update the SNR for calibration error
    # 0 = use only zero point error
    # 1 = add in quadrature 1/SNR and zero point error
    # 2 = use original snr, ignore zero point error
    if (options.snr_flag == 1):
        orig_err = 1.0 / stars.snr # flux
        new_err = np.sqrt(orig_err**2 + (zeropt_err/1.0857)**2)
        stars.snr = 1.0 / new_err
    if (options.snr_flag == 0):
        stars.snr = np.zeros(len(stars.snr), dtype=float) + (1.0857/zeropt_err)
    # Fix infinite SNR
    idx = np.where(stars.snr == float('inf'))
    if (len(idx) > 0):
        stars.snr[idx] = 9999.0

    ##########
    # *.zer file
    ##########
    _zer = open(options.zername, 'w')
    
    # Get the number of calibrators used:
    cidx = np.where((calibs.include == True) & (calibs.index >= 0))[0]
    calibCnt = len(cidx)

    _zer.write('# Original Calling Parameters:\n')
    _zer.write('# %s\n' % options.originalCall)
    _zer.write('#\n')
    _zer.write('#ZeroPoint   Error   Ncal   Cal names - ')
    _zer.write('{0}: {1}\n'.format(calibs.filt_name, calibs.filt_info))
    _zer.write('%10.3f   ' % zeropt)
    _zer.write('%5.3f   ' % zeropt_err)
    _zer.write('%4d   ' % calibCnt)
    for c in cidx:
        _zer.write('%s  ' % calibs.name[c])
    _zer.write('\n')

    _zer.close()


    ##########
    # *_cal.lis file
    ##########
    # Get the output order for the new starlist
    # If re-ordering, the calibrators we used are
    # first, in order.
    idx = np.arange(len(stars.name))

    fdx = []  # This will hold the reordered stuff that goes to the top.
    if (options.reorder):
        # Put the coo star at the top
        cdx = np.where(calibs.name == options.first_star)[0]
        fdx.append(calibs.index[cdx[0]])

        if (options.align_stars == None):
            # Used default calibraters
            cdx = np.where((calibs.include == True) & 
                           (calibs.index >= 0) &
                           (calibs.name != options.first_star))[0]
            fdx.extend(calibs.index[cdx])
        else:
            # Use in order specified by user.
            for c in range(len(options.align_stars)):
                # This should always work here since any issues
                # with the calib stars should have been caught 
                # earlier.
                if options.verbose:
                    print((options.align_stars[c]))
                
                if options.align_stars[c] == options.first_star:
                    continue

                ss = np.where(stars.name == options.align_stars[c])[0]
                if options.verbose:
                    print(ss)

                if len(ss) > 0:
                    fdx.append(ss[0])

        idx = np.delete(idx, fdx)
                
    # Now we have everything in order of fdx, idx
    _out = open(options.outname, 'w')
    
    for ff in fdx:
        _out.write('%13s  %9.6f  %8.3f  %10.5f  %10.5f  ' % 
                   (stars.name[ff], 
                    stars.mag[ff],
                    stars.epoch[ff],
                    stars.x[ff], stars.y[ff]))

        if (options.data_type == 2):
            _out.write('%7.3f %7.3f  ' % 
                       (stars.xerr[ff], stars.yerr[ff]))

        _out.write('%10.2f  %9.2f  %8d  %13.3f\n' %
                   (stars.snr[ff], stars.corr[ff],
                    stars.nframes[ff], stars.fwhm[ff]))

    for ff in idx:
        _out.write('%13s  %9.6f  %8.3f  %10.5f  %10.5f  ' % 
                   (stars.name[ff], 
                    stars.mag[ff],
                    stars.epoch[ff],
                    stars.x[ff], stars.y[ff]))

        if (options.data_type == 2):
            _out.write('%7.3f %7.3f  ' % 
                       (stars.xerr[ff], stars.yerr[ff]))

        _out.write('%10.2f  %9.2f  %8d  %13.3f\n' %
                   (stars.snr[ff], stars.corr[ff],
                    stars.nframes[ff], stars.fwhm[ff]))

    _out.close()
    
def get_camera_type(fitsfile):
    """
    Helper class to get the calibrate camera type from the
    FITS header.
    """
    # First check the instrument
    hdr = pyfits.getheader(fitsfile)
    instrument = hdr.get('CURRINST')
    if (instrument == None):
       # OLD SETUP
       instrument = hdr.get('INSTRUME')
       
    if (instrument == None):
       # OSIRIS
       instrument = hdr.get('INSTR')

       if ('imag' in instrument):
          instrument = 'OSIRIS-imag'
    
    # Default is still NIRC2
    if (instrument == None): 
        instrument = 'NIRC2'

    # get rid of the whitespace
    instrument = instrument.strip()

    # Check NICMOS camera
    if instrument == 'NICMOS':
        camera = hdr.get('CAMERA')
        instrument += camera.strip()
        

    # Check NIRC2 camera
    if instrument == 'NIRC2':
        camera = hdr.get('CAMNAME')
        instrument += camera.strip()


    cameraInfo = {'NIRC-D79': 1,
                  'KCAM-AO': 2,
                  'SCAM-AO': 3,
                  'NIRC2narrow': 4, 
                  'NIRC2medium': 5, 
                  'NIRC2wide': 6,
                  'SCAM-AO unmagnified': 7,
                  'NIRC2-AO narrow, binned by 3': 8,
                  'Hokupaa+QUIRC': 9,
                  'MMT PISCES': 10,
                  'UKIDSS': 11,
                  'OSIRIS-imag': 16,
                  'LGSAO': 4,
                  }


    camera = cameraInfo.get(instrument)

    # Default value
    if camera == None:
        camera = 1

    return camera


class Starlist(object):
    """Use this to keep all the lists associated with
    our input starlist."""
    pass

if __name__ == '__main__':
    main()


