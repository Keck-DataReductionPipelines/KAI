# reduce_combo.py
# Nights: 2021-07-13
# ---
# Make combo images for NIRC2 observation
# ---
# Abhimat Gautam

# Import python and iraf modules
from pyraf import iraf as ir
import numpy as np
import os, sys

# Import our own custom modules
from nirc2.reduce import dar
from nirc2.reduce import data
from nirc2.reduce import util

# Root directory for the data release
dr_root_dir = '/g/ghez/data/dr/dr1/'

# Date name for the observations
date_name_0713 = '20210713nirc2'

# Directories of clean files
clean_dir_0713 = '{0}/clean/{1}/'.format(dr_root_dir, date_name_0713)

# Epoch name for the combo epoch
epoch_name = '20210713nirc2'

# Directory of combo files
combo_dir = '{0}/combo/{1}/'.format(dr_root_dir, epoch_name)

def make_combo():
    """
    Make combo image for the NIRC2 observation.
    
    Needs to have individual frames cleaned and Strehl calculated.
    Needs to be run from the combo epoch's reduce directory as
    current working directory.
    """
    
    # Make sure we have necessary weather info for the year
    dar.get_atm_conditions(2021)
    
    # Define list of files
    
    # Kp Imaging
    
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
    files_0713_pos2 = range(157, 159+1) + range(163, 165+1) +\
                      range(169, 171+1) + range(175, 177+1) +\
                      range(181, 183+1) + range(187, 189+1) +\
                      range(193, 195+1) + range(199, 201+1) +\
                      range(205, 207+1) + range(211, 213+1)
    
    files_0713 = files_0713_pos1 + files_0713_pos2
    
    print('Total Observations, {0} (Kp): {1}'.format(date_name_0713,
                                                     len(files_0713)))
    
    clean_dirs_0713 = [clean_dir_0713 for i in range(len(files_0713))]
    
    # Assemble final lists for combine
    files_all = files_0713
    print('Total Observations (Kp): {0}'.format(len(files_all)))
    
    clean_dirs = clean_dirs_0713
    
    # Combine data
    data.combine(files_all, 'kp', epoch_name,
                 weight = 'strehl', submaps = 3, trim = True,
                 clean_dirs = clean_dirs, combo_dir = combo_dir)
    
    
    # Lp Imaging
    
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
    files_0713_pos2 = [156] +\
                      range(160, 162+1) + range(166, 168+1) +\
                      range(172, 174+1) + range(178, 180+1) +\
                      range(184, 186+1) + range(190, 192+1) +\
                      range(196, 198+1) + range(202, 204+1) +\
                      range(208, 210+1) + [214]
    
    files_0713 = files_0713_pos1 + files_0713_pos2
    
    print('Total Observations, {0} (Lp): {1}'.format(date_name_0713,
                                                     len(files_0713)))
    
    clean_dirs_0713 = [clean_dir_0713 for i in range(len(files_0713))]
    
    # Assemble final lists for combine
    files_all = files_0713
    print('Total Observations (Lp): {0}'.format(len(files_all)))
    
    clean_dirs = clean_dirs_0713
    
    # Combine data
    data.combine(files_all, 'lp', epoch_name,
                 weight = 'strehl', submaps = 3, trim = True,
                 clean_dirs = clean_dirs, combo_dir = combo_dir)
    