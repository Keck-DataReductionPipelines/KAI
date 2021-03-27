from astropy.table import Table
import os, errno, shutil
from astropy.io import fits as pyfits
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

def cp_change_prefix(arg1,arg2):
    """
    Takes files beginning with arg1 and replaces them with arg2
    Must be in the directory where files live
    """

    # Find files in this directory beginning with arg1
    files = os.listdir(".")
    # Ignore files beginning with '.'
    files=[filename for filename in files if filename[0] != '.']

    ln = len(arg1)

    for ff in range(len(files)):
        pre = files[ff][0:ln]
        if pre == arg1:
            suf = files[ff][len(arg1):]
            newFile = arg2 + suf
            shutil.copy(files[ff], newFile)

    return

def cp_change_suffix(arg1,arg2):
    """
    Takes files ending with arg1 and replaces them with arg2
    Must be in the directory where files live
    """

    # Find files in this directory ending with arg1
    files = os.listdir(".")
    # Ignore files beginning with '.'
    files=[filename for filename in files if filename[0] != '.']

    ln = len(arg1)

    for ff in range(len(files)):
        suf = files[ff][len(files[ff])-len(arg1):]
        if suf == arg1:
            pre = files[ff][0:len(files[ff])-len(arg1)]
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
        coo = Table.read(files[ff]+'.coo', format='ascii', header_start=None)
        coo_cols = list(coo.columns.keys())
        xref = coo[coo_cols[0]]
        yref = coo[coo_cols[1]]

        # Open .coord file and read strehl source's coordinates
        coord = Table.read(files[ff]+'.coord', format='ascii', header_start=None)
        coord_cols = list(coord.columns.keys())
        xstr = coord[coord_cols[0]]
        ystr = coord[coord_cols[1]]
 
        # Open image and write reference star x,y to fits header
        fits = pyfits.open(files[ff]+'.fits')

        fits[0].header.update('XREF', "%.3f" %xref,
                              'Cross Corr Reference Src x')
        fits[0].header.update('YREF', "%.3f" %yref,
                              'Cross Corr Reference Src y')
        fits[0].header.update('XSTREHL', "%.3f" %xstr,
                              'Strehl Reference Src x')
        fits[0].header.update('YSTREHL', "%.3f" %ystr,
                              'Strehl Reference Src y')

        # Output fits file
        _out = 'new_hdr/' + files[ff] + '.fits'
        fits[0].writeto(_out, output_verify='silentfix')

    return
