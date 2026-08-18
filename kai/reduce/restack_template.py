#!/usr/bin/env python

# restack.py
# ---
# Re-stack coadds

# Modification history
# Sean Granados: 08/18/2025 -- created file

"""
This script is meant to be ran in a raw directory containing the unprocessed fits files from NIRC2
(e.g. n_unp_0500.fits). 

Ensure the raw and clean directory of the processed fits files (e.g. n0500.fits and c0500.fits respectively) 
exists before running restack.py as it uses the .coo files created in the clean script as a starting point 
in re-centering the coadds.

{dr_root_dir}: root directory where the data exists

Inputs:
clean_directory (string): location of the clean directory of the processed fits files
raw_directory (string): location of the raw directory of the processed fits files
raw_unp_directory (string): location of the raw directory of the un-processed fits files which should
contain the restack.py script

Outputs:
restack.log (text file): text file with information on the x and y shifts of the coadds
n####.fits (fits file): restacked fits file

optional outputs:
restack_diagnostics (directory): directory containing plots and gif images of the coadd shifts
coadds.fits (fits file): fits file of each coadd

-----------
An example setup may be:
------------
*The directories below should already exist*

The following contains n####.fits files
{dr_root_dir}/raw/20250512nirc2/

The following contains c####.fits and .coo files where the clean script was ran
{dr_root_dir}/clean/20250512nirc2/

*The directory below should be created*

The following contains n_unp_####.fits files and restack.py, the restacked coadds (outputted as n####.fits) will be
placed here by the script restack.py
{dr_root_dir}/raw/20250512nirc2_st/
"""

# Import our own custom modules
from kai.reduce import kai_util

# Root directory
dr_root_dir = '/g/ghez/data/dr/dr2/'

# Date name for the observations containing 
date_name = '20250512nirc2'

# Stacked date name
stacked_name = '20250512nirc2_st'

# clean directory where the processed fits files were cleaned
clean_directory = '{0}/clean/{1}'.format(dr_root_dir, date_name)

# raw directory where the processed fits files live
raw_directory = '{0}/raw/{1}'.format(dr_root_dir, date_name)

# raw_unp_directory where the unprocessed fits files live and the restacked coadds will go
raw_unp_directory = '{0}/raw/{1}'.format(dr_root_dir, stacked_name)


# flags

# save_diagnostics (boolean): save some diagnostics like a coadd gif and plot showing the positional shift
save_diagnostics = True

# save_coadds (boolean): save each coadd as a fits file
save_coadds = False

# Call restackCoadd function in kai_util
kai_util.restackCoadd(clean_directory, raw_directory, raw_unp_directory, save_diagnostics = save_diagnostics, save_coadds = save_coadds)
