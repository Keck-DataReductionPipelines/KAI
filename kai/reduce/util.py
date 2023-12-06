from astropy.table import Table
import os, errno, shutil
from astropy.io import fits as fits
import numpy as np
import pdb

# Load up directory aliases
module_dir = os.path.dirname(__file__)
dir_alias_file = module_dir + '/../data/directory_aliases.txt'
dir_alias = Table.read(dir_alias_file, format='ascii.fast_no_header')
dir_alias.rename_column('col1', 'dir')
dir_alias.rename_column('col2', 'alias')


def rmall(files):
    """Remove list of files without confirmation."""
    for file in files:
        if os.access(file, os.F_OK): os.remove(file)

    return


def mkdir(dir):
    """Make directory if it doesn't already exist."""
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

    return


def getcwd():
    """
    IRAF doesn't like long file names. This reduces them.
    """
    curdir = os.getcwd()

    for ii in range(len(dir_alias)):
        curdir = curdir.replace(dir_alias['dir'][ii], dir_alias['alias'][ii])

    curdir += '/'

    return curdir


def trimdir(olddir):
    """
    IRAF doesn't like long file names. This reduces them.
    """
    for ii in range(len(dir_alias)):
        olddir = olddir.replace(dir_alias['dir'][ii], dir_alias['alias'][ii])

    return olddir


def cp_change_prefix(arg1, arg2):
    """
    Takes files beginning with arg1 and replaces them with arg2
    Must be in the directory where files live
    """

    # Find files in this directory beginning with arg1
    files = os.listdir(".")
    # Ignore files beginning with '.'
    files = [filename for filename in files if filename[0] != '.']

    ln = len(arg1)

    for ff in range(len(files)):
        pre = files[ff][0:ln]
        if pre == arg1:
            suf = files[ff][len(arg1):]
            newFile = arg2 + suf
            shutil.copy(files[ff], newFile)

    return


def cp_change_suffix(arg1, arg2):
    """
    Takes files ending with arg1 and replaces them with arg2
    Must be in the directory where files live
    """

    # Find files in this directory ending with arg1
    files = os.listdir(".")
    # Ignore files beginning with '.'
    files = [filename for filename in files if filename[0] != '.']

    ln = len(arg1)

    for ff in range(len(files)):
        suf = files[ff][len(files[ff]) - len(arg1):]
        if suf == arg1:
            pre = files[ff][0:len(files[ff]) - len(arg1)]
            newFile = pre + arg2
            shutil.copy(files[ff], newFile)

    return


def update_header_coords(fileList):
    """
    Updates coordinates in the header for XREF, YREF
    and XSTREHL, and YSTREHL.
 
    fileList : list of files to update
    """

    _files = Table.read(fileList, format='ascii', header_start=None)
    cols = list(_files.columns.keys())
    files = _files[cols[0]]
    files = [files[ff].split('.')[0] for ff in range(len(files))]

    for ff in range(len(files)):
        # Open .coo file and read 16C's coordinates
        coo = Table.read(files[ff] + '.coo', format='ascii', header_start=None)
        coo_cols = list(coo.columns.keys())
        xref = coo[coo_cols[0]]
        yref = coo[coo_cols[1]]

        # Open .coord file and read strehl source's coordinates
        coord = Table.read(files[ff] + '.coord', format='ascii', header_start=None)
        coord_cols = list(coord.columns.keys())
        xstr = coord[coord_cols[0]]
        ystr = coord[coord_cols[1]]

        # Open image and write reference star x,y to fits header
        fits = fits.open(files[ff] + '.fits')

        fits[0].header.update('XREF', "%.3f" % xref,
                              'Cross Corr Reference Src x')
        fits[0].header.update('YREF', "%.3f" % yref,
                              'Cross Corr Reference Src y')
        fits[0].header.update('XSTREHL', "%.3f" % xstr,
                              'Strehl Reference Src x')
        fits[0].header.update('YSTREHL', "%.3f" % ystr,
                              'Strehl Reference Src y')

        # Output fits file
        _out = 'new_hdr/' + files[ff] + '.fits'
        fits[0].writeto(_out, output_verify='silentfix')

    return


def imarith(img_list1, operator, img_list2, out_list):
    """
    Replace the IRAF imarith functionality with exactly the same
    thing implemented in python.
    """
    # Check to see if we are inputting a text file
    if np.isscalar(img_list1) and img_list1.startswith('@'):
        img_list1 = np.loadtxt(img_list1[1:], dtype="U")

    if np.isscalar(img_list2) and img_list2.startswith('@'):
        img_list2 = np.loadtxt(img_list2[1:], dtype="U")

    if np.isscalar(out_list) and out_list.startswith('@'):
        out_list = np.loadtxt(out_list[1:], dtype="U")

    # Check that everything is array like.
    if np.isscalar(img_list1):
        img_list1 = np.array([img_list1])

    if np.isscalar(img_list2):
        img_list2 = np.array([img_list2])

    if np.isscalar(out_list):
        out_list = np.array([out_list])

    # Check that we don't have lists with size>1 and of different sizes.
    l1 = len(img_list1)
    l2 = len(img_list2)
    lo = len(out_list)

    if l1 > 1 and l2 > 1 and l1 != l2:
        msg = 'Mis-matched lengths of file lists: img_list1 = %d, img_list2 = %d'.format(l1, l2)
        raise RuntimeError(msg)
    if l1 > 1 and lo > 1 and l1 != lo:
        msg = 'Mis-matched lengths of file lists: img_list1 = %d, out_list = %d'.format(l1, lo)
        raise RuntimeError(msg)
    if lo > 1 and l2 > 1 and lo != l2:
        msg = 'Mis-matched lengths of file lists: out_list = %d, img_list2 = %d'.format(lo, l2)
        raise RuntimeError(msg)
    if lo != np.max([l1, l2]):
        msg = 'Mis-matched lengths of file lists: out_list = %d, img_list1 = %d, img_list2 = %d'.format(lo, l1, l2)
        raise RuntimeError(msg)

    # Loop through files and apply the operations.
    for ii in range(lo):
        if l1 == lo:
            img1, hdr1 = fits.getdata(img_list1[ii], header=True)
        elif ii == 0:
            img1, hdr1 = fits.getdata(img_list1[0], header=True)

        if l2 == lo:
            img2 = fits.getdata(img_list2[ii])
        elif ii == 0:
            img2 = fits.getdata(img_list2[0])

        if operator == '+':
            out = img1 + img2
        if operator == '-':
            out = img1 - img2
        if operator == '*':
            out = img1 * img2
        if operator == '/':
            out = img1 / img2

        fits.writeto(out_list[ii], out, header=hdr1, output_verify='ignore', overwrite=True)

    return