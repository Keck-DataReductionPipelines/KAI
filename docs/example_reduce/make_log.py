#!/usr/bin/env python

# make_log.py
# ---
# Make log of all FITS files for the observation

# Import our own custom modules
from kai.reduce import kai_util

# Root directory for the data release
dr_root_dir = '/g/ghez/data/dr/dr1/'

# Date name for the observations
date_name = '20210713nirc2'

kai_util.kailog('{0}/raw/{1}/'.format(dr_root_dir, date_name))
