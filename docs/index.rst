********************************************
Python NIRC2 Reduction Pipeline (``nirc2``)
********************************************

.. _numpy: http://www.numpy.org

Introduction
============

This Astropy-affiliated package provides a reduction pipeline for images
taken with the NIRC2 instrument at the W. M. Keck Observatory. 

Download
========

The latest release can be downloaded from the github repository `here
<https://github.com/jluastro/nirc2>`_.

Installation
============

This pacakge depends on pyraf and astropy functionality and has been 
tested on top of `Ureka <http://ssb.stsci.edu/ureka/>`_ python/pyraf installation.
The respository can be cloned from github and by installing with

python setup.py build
python setup.py install


Using `nirc2`
=============

The NIRC2 module is imported using::

  >>> import nirc2

Within a top-level directory for a NIRC2 run (i.e. ``2012jun``), create
the following directories::

  mkdir reduce
  mkdir clean
  mkdir combo
  mkdir raw

All ``n*.fits`` files should be located within raw (not in sub-directorities). 
If processing several consecutive nights worth of data, and file names are not
unique, then see the helper utility `~nirc2.reduce.util.cp_change_prefix` for 
renaming files. We recommend files from night 1 be named with `n0*.fits`, 
files from night 2 be named with `n1*.fits`, etc.

Copy the reduction template script `~nirc2.reduce.reduce_template.py`
to `reduce/reduce.py`. Modify the reduction script to point to the appropriate
FITS file numbers for calibration frames, sky frames, and science frames.
See `Reduction Template`_ for a more complete description of the template 
and example reductions.

A NIRC2 reduction can be run from the `reduce/` directory in the following 
manner::

  >>> import reduce
  >>> reduce.makelog()
  >>> reduce.go()

Output files from the reduction are described in `Output Files from Reduction`_.


Reduction Template
==================


Output Files from Reduction
===========================

