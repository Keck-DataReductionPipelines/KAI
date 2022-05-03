# from pyraf import iraf
import glob
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
import math
from astropy.io import fits as pyfits
import datetime
try:
    import urllib.request, urllib.parse, urllib.error
    p2 = False
except:
    import urllib
    p2 = True # Python 2 is running, so the urllib command is different
import os, sys
from kai.reduce import kai_util
from kai.reduce import util
from kai.reduce import slalib
from kai import instruments
from astropy.table import Table

module_dir = os.path.dirname(__file__)


def get_atm_conditions(year):
    """
    Retrieve atmospheric conditions from CFHT archive website,
    then calls dar.splitAtmosphereCFHT() to separate the data
    by months.
    """
    yearStr = str(year)

    if p2: # Python 2 command, necessary for IRAF
        _atm = urllib.urlopen("http://mkwc.ifa.hawaii.edu/archive/wx/cfht/cfht-wx.%s.dat" % yearStr)
    else:
        _atm = urllib.request.urlopen("http://mkwc.ifa.hawaii.edu/archive/wx/cfht/cfht-wx.%s.dat" % yearStr)
    atm = _atm.read()
    _atm.close()
    
    root = module_dir + '/weather/'

    if type(atm) == bytes:
        # this is for python 3
        atmfile = open(root + 'cfht-wx.' + yearStr + '.dat','wb')
    else:
        atmfile = open(root + 'cfht-wx.' + yearStr + '.dat','w')


    atmfile.write(atm)
    atmfile.close()

    splitAtmosphereCFHT(str(year))


def keckDARcoeffs(lamda, year, month, day, hour, minute):
    """
    Calculate the differential atmospheric refraction
    for two objects observed at Keck.

    Input:
    lamda -- Effective wavelength (microns) assumed to be the same for both
    year, month, day, hour, minute of observation (HST)

    Output:
    refA
    refB
    """
    # iraf.noao()

    # Set up Keck observatory info
    # foo = iraf.noao.observatory(command="set", obsid="keck", Stdout=1)
    # obs = iraf.noao.observatory

    ####################
    # Setup all the parameters for the atmospheric refraction
    # calculations. Typical values obtained from the Mauna Kea
    # weather pages and from the web.
    ####################
    # 
    # Temperature Lapse Rate (Kelvin/meter)
    tlr = 0.0065
    # Precision required to terminate the iteration (radian)
    eps = 1.0e-9
    # Height above sea level (meters)
    # hm = obs.altitude
    hm = 4160.0                             #Keck height. Hard coded to remove iraf dependency
    # Latitude of the observer (radian)
    # phi = math.radians(obs.latitude)
    phi = math.radians(19.82833333333333)   #Keck latitude. Hard coded to remove iraf dependency

    # Pull from atmosphere logs.
    logDir = module_dir + '/weather/'
    logFile = logDir +'cfht-wx.'+ str(year) +'.'+ str(month).zfill(2) +'.dat'
    print(logFile)

    _atm = Table.read(logFile, format='ascii', header_start=None)
    atmYear = _atm['col1']
    atmMonth = _atm['col2']
    atmDay = _atm['col3']
    atmHour = _atm['col4']
    atmMin = _atm['col5']       # HST times
    atmTemp = _atm['col8']      # Celsius
    atmHumidity = _atm['col9']  # percent
    atmPressure = _atm['col10'] # mb pressure

    # Find the exact time match for year, month, day, hour
    idx = (np.where((atmYear == year) & (atmMonth == month) &
                    (atmDay == day) & (atmHour == hour)))[0]
    
    if (len(idx) == 0):
        print(( 'Could not find DAR data for %4d-%2d-%2d %2d:%2d in %s' % \
            (year, month, day, hour, minute, logFile)))

    atmMin = atmMin[idx]
    atmTemp = atmTemp[idx]
    atmHumidity = atmHumidity[idx]
    atmPressure = atmPressure[idx]

    # Find the closest minute
    minDiff = abs(atmMin - minute)
    sdx = minDiff.argsort()

    # Select out the closest in time.
    # Ambient Temperature (Kelvin)
    # Should be around 274.0 Kelvin
    tdk = atmTemp[sdx[0]] + 272.15
    # Pressure at the observer (millibar)
    # Should be around 760.0 millibars
    pmb = atmPressure[sdx[0]]
    # Relative humidity (%)
    # Should be around 0.1 %
    rh = atmHumidity[sdx[0]] / 100.0 #relative humidity should be between 0 and 1
    print(hm, tdk, pmb, rh, lamda, phi, tlr, eps)
    return slalib.refco(hm, tdk, pmb, rh, lamda, phi, tlr, eps)

def kaidar(fitsFile, instrument=instruments.default_inst):
    """
    Use the FITS header to extract date, time, wavelength,
    elevation, and image orientation information. This is everything
    that is necessary to calculate the differential atmospheric
    refraction. The differential atmospheric refraction 
    is applicable only along the zenith direction of an image.
    This code calculates the predicted DAR using archived CFHT
    atmospheric data and the elevation and wavelength of the observations.
    Then the DAR correction is transformed into image coefficients that
    can be applied in image coordinates. 
    """
    # Get header info
    img, hdr = pyfits.getdata(fitsFile, header=True)

    effWave = instrument.get_central_wavelength(hdr)
    elevation = hdr[instrument.hdr_keys['elevation']]
    # airmass = hdr['AIRMASS']
    # parang = hdr['PARANG']

    date = hdr['DATE-OBS'].split('-')
    year = int(date[0])
    month = int(date[1])
    day = int(date[2])

    utc = hdr['UTC'].split(':')
    hour = int(utc[0])
    minute = int(utc[1])
    second = int(math.floor(float(utc[2])))

    utc = datetime.datetime(year, month, day, hour, minute, second)
    utc2hst = datetime.timedelta(hours=-10)
    hst = utc + utc2hst

    (refA, refB) = keckDARcoeffs(effWave, hst.year, hst.month, hst.day,
                                 hst.hour, hst.minute)

    tanz = math.tan(math.radians(90.0 - elevation))
    tmp = 1.0 + tanz**2
    darCoeffL = tmp * (refA + 3.0 * refB * tanz**2)   #unitless
    darCoeffQ = -tmp * (refA*tanz +
                            3.0 * refB * (tanz + 2.0*tanz**3))   #units: radians^-1

    # Convert DAR coefficients for use with arcseconds
    # scale = instrument.get_plate_scale(hdr)
    # darCoeffQ *= 1.0 * scale / 206265.0
    darCoeffL *= 1.0
    darCoeffQ *= 1.0 / 206265.0             #units now arcsec^-1
   
    # # Lets determine the zenith and horizon unit vectors for
    # # this image.
    # pos_ang = instrument.get_position_angle(hdr)
    # pa = math.radians(parang + pos_ang)
    # zenithX = -math.sin(pa)
    # zenithY = math.cos(pa)

    # # Compute the predicted differential atmospheric refraction
    # # over a 10'' seperation along the zenith direction.
    # # Remember coeffecicents are only for deltaZ in pixels
    # deltaZ = img.shape[0] * scale
    # deltaR = darCoeffL * (deltaZ/scale) + darCoeffQ * (deltaZ/scale)**2
    # deltaR *= scale # now in arcseconds

    # magnification = (deltaZ + deltaR) / deltaZ
    # print(( 'DAR FITS file = %s' % (fitsFile)))
    # print(('DAR over 10": Linear dR = %f"  Quad dR = %f"' % \
    #       (darCoeffL * deltaZ, darCoeffQ * deltaZ**2)))
    # print(('DAR Magnification = %f' % (magnification)))
    # print(('DAR Vertical Angle = %6.1f' % (math.degrees(pa))))

    return (darCoeffL, darCoeffQ)

def darPlusDistortion(inputFits, outputRoot, xgeoim=None, ygeoim=None, instrument=instruments.default_inst):
    """
    Create lookup tables (stored as FITS files) that can be used
    to correct DAR. Optionally, the shifts due to DAR can be added
    to existing NIRC2 distortion lookup tables if the xgeoim/ygeoim
    input parameters are set.

    Inputs:
    inputFits - a NIRC2 image for which to determine the DAR correction
    outputRoot - the root name for the output. This will be used as the
        root name of two new images with names, <outputRoot>_x.fits and 
        <outputRoot>_y.fits.

    Optional Inputs:
    xgeoim/ygeoim - FITS images used in Drizzle distortion correction
        (lookup tables) will be modified to incorporate the DAR correction.
        The order of the correction is 1. distortion, 2. DAR.
        
    """
    # Get the size of the image and the half-points
    hdr = pyfits.getheader(inputFits)
    imgsizeX = int(hdr['NAXIS1'])
    imgsizeY = int(hdr['NAXIS2'])
    halfX = int(round(imgsizeX / 2.0))
    halfY = int(round(imgsizeY / 2.0))

    # First get the coefficients
    (darCoeffL, darCoeffQ) = kaidar(inputFits, instrument=instrument)
    # Convert DAR coefficients for use with units of NIRC2 pixels
    scale = instrument.get_plate_scale(hdr)
    darCoeffL *= 1.0                
    darCoeffQ *= 1.0 * scale #/ 206265.0    
    pa = math.radians(instrument.get_parallactic_angle(hdr) + instrument.get_position_angle(hdr))
    #(a, b) = kaidarPoly(inputFits)

    # Create two 1024 arrays (or read in existing ones) for the
    # X and Y lookup tables
    if ((xgeoim == None) or (xgeoim == '')):
        x = np.zeros((imgsizeY, imgsizeX), dtype=float)
    else:
        x = pyfits.getdata(xgeoim)

    if ((ygeoim == None) or (ygeoim == '')):
        y = np.zeros((imgsizeY, imgsizeX), dtype=float)
    else:
        y = pyfits.getdata(ygeoim)

    # Get proper header info.
    fits = pyfits.open(inputFits)

    axisX = np.arange(imgsizeX, dtype=float) - halfX
    axisY = np.arange(imgsizeY, dtype=float) - halfY
    xcoo2d, ycoo2d = np.meshgrid(axisX, axisY)

    xnew1 = xcoo2d + x
    ynew1 = ycoo2d + y

    # Rotate coordinates clockwise by PA so that zenith is along +ynew2
    # PA = parallactic angle (angle from +y to zenith going CCW)
    sina = math.sin(pa)
    cosa = math.cos(pa)

    xnew2 = xnew1 * cosa + ynew1 * sina
    ynew2 = -xnew1 * sina + ynew1 * cosa

    # Apply DAR correction along the y axis
    xnew3 = xnew2
    ynew3 = ynew2*(1 + darCoeffL) + ynew2*np.abs(ynew2)*darCoeffQ

    # Rotate coordinates counter-clockwise by PA back to original
    xnew4 = xnew3 * cosa - ynew3 * sina
    ynew4 = xnew3 * sina + ynew3 * cosa

    #xnew2 = a[0] + a[1]*xnew1 + a[2]*ynew1 + \
    #        a[3]*xnew1**2 + a[4]*xnew1*ynew1 + a[5]*ynew1**2
    #ynew2 = b[0] + b[1]*xnew1 + b[2]*ynew1 + \
    #        b[3]*xnew1**2 + b[4]*xnew1*ynew1 + b[5]*ynew1**2

    x = xnew4 - xcoo2d
    y = ynew4 - ycoo2d

    xout = outputRoot + '_x.fits'
    yout = outputRoot + '_y.fits'
    util.rmall([xout, yout])
    fits[0].data = x
    fits[0].writeto(xout, output_verify='silentfix')
    fits[0].data = y
    fits[0].writeto(yout, output_verify='silentfix')

    return (xout, yout)


def applyDAR(inputFits, spaceStarlist, plot=False, instrument=instruments.default_inst, plotdir = './'):
    """
    inputFits: (str) name of fits file associated with this starlist

    spaceStarlist: (astropy table) must include columns 'x0' and 'y0'.

    Input a starlist in x=RA (+x = west) and y=Dec (arcseconds) taken from
    space and introduce differential atmospheric refraction (DAR). The amount
    of DAR that is applied depends on the header information in the input fits
    file. The resulting output starlist should contain what was observed
    after the starlight passed through the atmosphere, but before the
    starlight passed through the telescope. Only achromatic DAR is 
    applied in this code.

    returns spaceStarlist with updated 'x0' and 'y0'
    """

    # Get header info
    #hdr = pyfits.getheader(fits)

    #effWave = hdr['EFFWAVE']
    #elevation = hdr['EL']
    #lamda = hdr['CENWAVE']
    #airmass = hdr['AIRMASS']
    #parang = hdr['PARANG']
    #
    #date = hdr['DATE-OBS'].split('-')
    #year = int(date[0])
    #month = int(date[1])
    #day = int(date[2])

    #utc = hdr['UTC'].split(':')
    #hour = int(utc[0])
    #minute = int(utc[1])
    #second = int(math.floor(float(utc[2])))

    #utc = datetime.datetime(year, month, day, hour, minute, second)
    #utc2hst = datetime.timedelta(hours=-10)
    #hst = utc + utc2hst


    #(refA, refB) = keckDARcoeffs(effWave, hst.year, hst.month, hst.day,
    #                             hst.hour, hst.minute)

    #tanz = math.tan(math.radians(90.0 - elevation))
    #tmp = 1.0 + tanz**2
    #darCoeffL = tmp * (refA + 3.0 * refB * tanz**2)
    #darCoeffQ = -tmp * (refA*tanz +
    #                        3.0 * refB * (tanz + 2.0*tanz**3))

    # Convert DAR coefficients for use with arcseconds
    #darCoeffL *= 1.0
    #darCoeffQ *= 1.0 / 206265.0
    
    # Lets determine the zenith and horizon unit vectors for
    # this image. The angle we need is simply the parallactic
    # (or vertical) angle since ACS images are North Up already.
    #pa = math.radians(parang)

    #MS: presumably the above code is all replacable with this call (which uses the intrument object
    (darCoeffL, darCoeffQ) = kaidar(inputFits, instrument=instrument)
    hdr = pyfits.getheader(inputFits)
    pa = math.radians(instrument.get_parallactic_angle(hdr))
    ##########
    #
    # Read in the starlist
    #
    ##########
    # _list = Table.read(spaceStarlist, format='ascii')
    # cols = list(_list.columns.keys())
    # names = [_list[ss][0].strip() for ss in range(len(_list))]
    # mag = _list[cols[1]]
    # date = _list[cols[2]]
    # x = _list[cols[3]] # RA in arcsec
    # y = _list[cols[4]]
    # xe = _list[cols[5]]
    # ye = _list[cols[6]]

    x = spaceStarlist['x0']     #already in arcseconds, x increasing to west.
    y = spaceStarlist['y0']

    # Magnify everything in the y (zenith) direction. Do it relative to
    # the first star. Even though dR depends on dzObs (ground observed dz),
    # it is a small mistake and results in less than a 10 micro-arcsec
    # change in dR.
    dx = x - x[0]
    dy = y - y[0]

    # Rotate coordinates CW so that the zenith angle is at +ynew
    sina = math.sin(pa)
    cosa = math.cos(pa)
    xnew1 = dx * cosa + dy * sina
    ynew1 = -dx * sina + dy * cosa

    # Apply DAR
    xnew2 = xnew1
    ynew2 = ynew1 * (1.0 - darCoeffL) - ynew1 * np.abs(ynew1) * darCoeffQ

    # Rotate coordinates CCW back to original angle
    xnew3 = xnew2 * cosa - ynew2 * sina
    ynew3 = xnew2 * sina + ynew2 * cosa

    xnew = xnew3 + x[0]
    ynew = ynew3 + y[0]


    ##########
    #
    # Write out the starlist
    #
    ##########
    # Save the current directory
    # newFits = fits.replace('.fits', '').split('/')[-1]
    # newList = newFits + '_acs.lis'
    # print(newList)
    # _new = open(newList, 'w')
    # for i in range(len(names)):
    #     _new.write('%10s  %7.3f  %7.2f  %10.4f  %10.4f  0  0  10  1  1  8\n' % \
    #           (names[i], mag[i], date[i], xnew[i], ynew[i]))
    # _new.close()

    # if (plot==True):
    #     py.clf()
    #     py.quiver(x, y, xnew - x, ynew - y, scale=0.02)
    #     py.quiver([0], [0], [0.001], [0], color='r', scale=0.02)
    #     py.axis([-5, 5, -5, 5])
    #     py.show()
        
    if (plot==True):
        plt.close('all')
        q = plt.quiver(x, y, xnew - x, ynew - y, scale=0.02)
        plt.xlabel('RA (arcsec)')
        plt.ylabel('DEC (arcsec)')
        plt.xlabel(r'$\Delta$ RA * cos(DEC) (arcsec)')
        plt.ylabel(r'$\Delta$ DEC (arcsec)')
        plt.figtext(0.15,0.85,'q =' + str(round(math.degrees(pa))))
        quiv_label = '1 mas'
        quiv_label_val = 0.001
        plt.quiverkey(q, 0.80, 0.85, quiv_label_val, quiv_label, coordinates='figure', labelpos='E', color='green')     
        plt.savefig(plotdir+ inputFits[-26:-5] + '.jpg', bbox_inches='tight')
        # plt.show()

    spaceStarlist['x0'] = xnew
    spaceStarlist['y0'] = ynew 

    return spaceStarlist

def splitAtmosphereCFHT(year):
    """
    Take an original archive file containing atmospheric parameters and
    split it up into seperate files for individual months. This makes
    later calls to calculate DAR parameters MUCH faster.
    """
    yearStr = str(year)
    logDir = module_dir + '/weather/'
    logFile = logDir + '/cfht-wx.' + yearStr + '.dat'

    _infile = open(logFile, 'r')

    outfiles = []
    for ii in range(1, 12+1):
        monthStr = str(ii).zfill(2)
        _month = open(logDir + '/cfht-wx.' +yearStr+ '.' +monthStr+ '.dat', 'w')
        outfiles.append( _month )

    for line in _infile:
        fields = line.split()

        month = int(fields[1])
        
        # Check for wrong month number
        if not (month > 0 and month <= 12):
            continue
        
        _outfile = outfiles[month-1]
        _outfile.write(line)

    for _month in outfiles:
        _month.close()
        
        
    
def test_darPlusDistortion():
    data_dir = module_dir + '/distortion/'
    file_geox_darunfix = data_dir + 'nirc2dist_xgeoim.fits'
    file_geoy_darunfix = data_dir + 'nirc2dist_ygeoim.fits'

    data_dir = '/u/ghezgroup/data/m92_test/08jul_new_on/'
    file_geox_darfix = data_dir + 'reduce/kp/gc_f1/ce0249geo_x.fits'
    file_geoy_darfix = data_dir + 'reduce/kp/gc_f1/ce0249geo_y.fits'

    xon = pyfits.getdata(file_geox_darfix)
    yon = pyfits.getdata(file_geoy_darfix)
    xoff = pyfits.getdata(file_geox_darunfix)
    yoff = pyfits.getdata(file_geoy_darunfix)

    # Make arrays with the coordinates for each 
    imgsize = 1024
    axisX = np.arange(imgsize, dtype=float)
    axisY = np.arange(imgsize, dtype=float)
    xcoo2d, ycoo2d = np.meshgrid(axisX, axisY)

    # Lets trim so that we only keep every 20th pixel
    idx = np.arange(25, imgsize, 50)
    xon = xon.take(idx, axis=0).take(idx, axis=1)
    yon = yon.take(idx, axis=0).take(idx, axis=1)
    xoff = xoff.take(idx, axis=0).take(idx, axis=1)
    yoff = yoff.take(idx, axis=0).take(idx, axis=1)
    xcoo2d = xcoo2d.take(idx, axis=0).take(idx, axis=1)
    ycoo2d = ycoo2d.take(idx, axis=0).take(idx, axis=1)

    # Calculate differences
    xdiff = xon - xoff
    ydiff = yon - yoff

    # Make vector plots
    py.clf()
    qvr = py.quiver2([xcoo2d], [ycoo2d], [xdiff], [ydiff],
                     units='width', scale=5, 
                     width=0.005, headwidth=3, headlength=3, 
                     headaxislength=3)
    py.quiverkey(qvr, 100, 1120, 1.0, '1 pixel', coordinates='data', color='r')
    py.xlabel('NIRC2 X (pixel)')
    py.ylabel('NIRC2 Y (pixel)')
    py.title('Arrows point to DAR Fix')
    #py.savefig('plots/vector_daroffon.png')
    py.show()


