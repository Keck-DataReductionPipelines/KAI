# reduce_clean.py
# Night: 2021-07-13
# ---
# Make clean images for NIRC2 observation

# Import python and iraf modules
from pyraf import iraf as ir
import numpy as np
import os, sys

# Import our own custom modules
from kai.reduce import dar
from kai.reduce import calib
from kai.reduce import sky
from kai.reduce import data
from kai.reduce import util

# Root directory for the data release
dr_root_dir = '/g/ghez/data/dr/dr1/'

# Date name for the observations
date_name = '20210713nirc2'

# Directory of raw files
raw_dir = '{0}/raw/{1}/'.format(dr_root_dir, date_name)

# Directory of clean files
clean_dir = '{0}/clean/{1}/'.format(dr_root_dir, date_name)

def make_calib():
    """
    Make calibration files (darks, flats, masks)
    Needs to be run from reduce directory as current working directory
    """
    
    # Make sure we have necessary weather info for the year
    dar.get_atm_conditions(2021)
    
    # Calibration files
    
    # Darks - created in subdir darks/
    # Files: n0044 -- n0053
    #  - darks needed to make bad pixel mask
    #  - store the resulting dark in the file name that indicates the
    #    integration time (15s) and the coadds (1ca).
    dark_files = range(44, 53+1)
    calib.makedark(dark_files, 'dark_15s_1ca.fits',
                   raw_dir=raw_dir)
    
    # Kp Flats
    # Files: n0004 -- n0023
    # off, on, off, on ...
    off_files = range(4, 23+1, 2)
    on_files = range(5, 23+1, 2)
    calib.makeflat(on_files, off_files, 'flat_kp.fits',
                   raw_dir=raw_dir)
        
    # Make L-band flat from L skies (on) and darks (off)
    # On flats: use skies (n0216 -- n0235)
    # Off flats: use darks (n0044 -- n0053)
    on_files = range(216, 235+1)
    off_files = range(44, 53+1)
    calib.makeflat(on_files, off_files, 'flat_lp.fits',
                   raw_dir=raw_dir)
    
    # Masks
    calib.makemask('dark_15s_1ca.fits', 'flat_kp.fits',
                   'supermask.fits')


def make_sky():
    """
    Make the sky files
    """
    
    # Make sure we have necessary weather info for the year
    dar.get_atm_conditions(2021)
    
    # Galactic Center
    
    # Kp imaging skies
    
    # 2021-07-13
    
    # Sky Files
    # Sky files: n0058 -- n0068
    
    sky_files = range(58, 68+1)
    sky.makesky(sky_files, '0713', 'kp',
                raw_dir=raw_dir)
        
    # Lp imaging skies
    
    # 2021-07-13
    
    # Sky Files
    # Sky files: n0216 -- n0235
    
    sky_files = range(216, 235+1)
    sky.makesky_lp(sky_files, '0713', 'lp',
                raw_dir=raw_dir)

def make_clean():
    """
    Clean the science files
    """
    
    # Make sure we have necessary weather info for the year
    dar.get_atm_conditions(2021)
    
    # Galactic Center
    
    # Kp imaging data cleaning
    
    # 2021-07-13
    
    # Science Files
    # Files:
    # n0085 -- n0087, n0091 -- n0093,
    # n0097 -- n0099, n0103 -- n0105,
    # n0109 -- n0111, n0115 -- n0117,
    # n0121 -- n0123, n0127 -- n0129,
    # n0133 -- n0135, n0139 -- n0141,
    # n0145 -- n0147, [lost laser]
    # n0157 -- n0159, n0163 -- n0165,
    # n0169 -- n0171, n0175 -- n0177,
    # n0181 -- n0183, n0187 -- n0189,
    # n0193 -- n0195, n0199 -- n0201,
    # n0205 -- n0207, n0211 -- n0213
    
    files_0713_pos1 = range(85, 87+1) + range(91, 93+1) +\
                      range(97, 99+1) + range(103, 105+1) +\
                      range(109, 111+1) + range(115, 117+1) +\
                      range(121, 123+1) + range(127, 129+1) +\
                      range(133, 135+1) + range(139, 141+1) +\
                      range(145, 147+1)
    refSrc = [442, 700]     # IRS 16C
    strSrc = [537, 406]     # IRS 33N
    data.clean(files_0713_pos1, '0713', 'kp', refSrc, strSrc,
               raw_dir=raw_dir, clean_dir=clean_dir)
    
    files_0713_pos2 = range(157, 159+1) + range(163, 165+1) +\
                      range(169, 171+1) + range(175, 177+1) +\
                      range(181, 183+1) + range(187, 189+1) +\
                      range(193, 195+1) + range(199, 201+1) +\
                      range(205, 207+1) + range(211, 213+1)
    refSrc = [480, 678]     # IRS 16C
    strSrc = [575, 384]     # IRS 33N
    data.clean(files_0713_pos2, '0713', 'kp', refSrc, strSrc,
               raw_dir=raw_dir, clean_dir=clean_dir)
    
    
    files_0713 = files_0713_pos1 + files_0713_pos2
    
    print('Total Observations (Kp, 2021-07-13): {0}'.format(len(files_0713)))
    
    try:
        data.calcStrehl(files_0713, 'kp', clean_dir=clean_dir)
    except RuntimeError as e:
        print("Error during Strehl calculation:")
        print(e)
    
    
    # Lp imaging data cleaning
    
    # 2021-07-13
    
    # Science Files
    # Files:
    # n0088 -- n0090, n0094 -- n0096,
    # n0100 -- n0102, n0106 -- n0108,
    # n0112 -- n0114, n0118 -- n0120,
    # n0124 -- n0126, n0130 -- n0132,
    # n0136 -- n0138, n0142 -- n0144,
    # n0148, [lost laser]
    # n0156,
    # n0160 -- n0162, n0166 -- n0168,
    # n0172 -- n0174, n0178 -- n0180,
    # n0184 -- n0186, n0190 -- n0192,
    # n0196 -- n0198, n0202 -- n0204,
    # n0208 -- n0210, n0214
    
    files_0713_pos1 = range(88, 90+1) + range(94, 96+1) +\
                      range(100, 102+1) + range(106, 108+1) +\
                      range(112, 114+1) + range(118, 120+1) +\
                      range(124, 126+1) + range(130, 132+1) +\
                      range(136, 138+1) + range(142, 144+1) +\
                      [148]
    refSrc = [441, 696]     # IRS 16C
    strSrc = [537, 404]     # IRS 33N
    data.clean(files_0713_pos1, '0713', 'lp', refSrc, strSrc,
               raw_dir=raw_dir, clean_dir=clean_dir)
    
    files_0713_pos2 = [156] +\
                      range(160, 162+1) + range(166, 168+1) +\
                      range(172, 174+1) + range(178, 180+1) +\
                      range(184, 186+1) + range(190, 192+1) +\
                      range(196, 198+1) + range(202, 204+1) +\
                      range(208, 210+1) + [214]
    refSrc = [438, 670]     # IRS 16C
    strSrc = [529, 378]     # IRS 33N
    data.clean(files_0713_pos2, '0713', 'lp', refSrc, strSrc,
               raw_dir=raw_dir, clean_dir=clean_dir)
    
    files_0713 = files_0713_pos1 + files_0713_pos2
    
    print('Total Observations (Lp, 2021-07-13): {0}'.format(len(files_0713)))
    
    try:
        data.calcStrehl(files_0713, 'lp', clean_dir=clean_dir)
    except RuntimeError as e:
        print("Error during Strehl calculation:")
        print(e)
