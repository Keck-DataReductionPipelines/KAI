import os, shutil
from kai.reduce import util
from astropy.table import Table
import numpy as np
from astropy.io import fits
from astropy.table import Table
from kai.reduce import kai_util
from kai.reduce import calibrate
from kai.reduce import align_rms
from kai import instruments
from kai import strehl
from datetime import datetime
import subprocess
import pylab as py
import pdb
from multiprocessing import Pool

def idl_process_run(batch_file):
    """
    Function to help run parallel processes for IDL, used for StarFinder
    
    Parameters
    ----------
    batch_file : str
        Path of the bath file containing the IDL command to run
    """
    
    batch_file_log = batch_file + '.log'
    
    cmd = 'idl < ' + batch_file + ' >& ' + batch_file_log            
    #os.system(cmd)
    subp = subprocess.Popen(cmd, shell=True, executable="/bin/tcsh")
    tmp = subp.communicate()
    
    return

class Analysis(object):
    """
    Object that will perform our standard post-data-reduction analysis.
    This includes running starfinder, calibrating, and extracting positional
    and photometric errors via align_rms. 
    """

    def __init__(self, epoch, rootDir='/g/lu/data/orion/', filt='kp',
                 clean_dir = None, combo_dir = None,
                 combo_stf_dir = None,
                 epochDirSuffix=None, imgSuffix=None, stfDir=None,
                 useDistorted=False, cleanList='c.lis',
                 airopa_mode='single', stf_version=None, stf_debug=False,
                 stf_parallelize=True,
                 instrument=instruments.default_inst):
        """
        Set up Analysis object.
        
        Parameters
        ----------
        epoch : str
            Name of the combo epoch
        rootDir : str, default='/g/lu/data/orion/'
            Path of the root directory
        filt : str, default='kp'
            Name of the filter of reduction
        clean_dir : str, optional
            Directory where clean files are stored. By default,
            assumes that clean files are stored in '../clean'
        combo_dir : str, optional
            Directory where combo files are stored. By default,
            assumes that combo files are stored in '../combo'
        combo_stf_dir : str, optional
            Directory where starfinder files will be stored. By default,
            starfinder output is stored in '../combo/starfinder/'
        airopa_mode : str, default='single'
            The airopa mode to use. Available options are 'legacy', 'single',
            or 'variable'.
        stf_version : str, default=None
            If specified, adds a version suffix to the starfinder directory
            (e.g.: 'v3_1')
        instrument : instruments object, optional
            Instrument of data. Default is `instruments.default_inst`
        stf_debug : boolean, default=False
            Keyword to specify if starfinder should be run in debug mode
        stf_parallelize : boolean, default=True
            Keyword to specify if starfinder runs on main map and sub maps
            should be run in parallel.
            By default, runs in parallel (all 4 run together).
        """
        # Setup default parameters
        self.type = 'ao'
        self.corrMain = [0.8]
        self.corrSub = [0.6]
        self.corrClean = [0.7]
        # airopa mode can be "legacy", "single", or "variable"
        self.airopa_mode = airopa_mode
        self.trimfake = 1
        self.stfFlags = ''
        self.stf_debug = stf_debug
        self.stf_parallelize = stf_parallelize

        self.starlist = rootDir + 'source_list/psf_central.dat'
        self.labellist = rootDir+ 'source_list/label.dat'
        self.orbitlist = rootDir+ 'source_list/orbits.dat'
        self.calFile = rootDir + 'source_list/photo_calib.dat'

        # Keep track of the instrument these images are from.
        # We need this to get things like plate scale, etc.
        self.instrument = instrument
        
        self.calStars = [
            'irs16NW', 'S3-22', 'S2-22', 'S4-3', 'S1-1', 'S1-21', 'S1-12',
            'S2-2', 'S3-88', 'S2-75', 'S3-36', 'S1-33',
        ]
        self.calFlags = '-f 1 -R '
        self.mapFilter2Cal = {'kp': 'Kp', 'h': 'H', 'lp': 'Lp_o1', 'ms': 'Ms_o1'}
        if 'kp' in filt:
            self.calColumn = self.mapFilter2Cal['kp']
        elif 'J' in filt:
            self.calColumn = 'K_o2'
        elif filt in self.mapFilter2Cal: # case for Maser mosaic, deep mosaic
            self.calColumn = self.mapFilter2Cal[filt]
        else:                           #default to kp if filt not in mapFilter2Cal.  Consider other solution to deal with unexpected filt
            self.calColumn = self.mapFilter2Cal['kp']
        self.calCooStar = 'irs16C'
        
        self.cooStar = 'irs16C'

        self.deblend = None

        #self.alignFlags = '-R 3 -v -p -a 0'
        self.alignFlags = '-R 3 -p -a 0'

        self.numSubMaps = 3
        self.minSubMaps = 3
        
        # Setup input parameters 
        self.epoch = epoch

        self.rootDir = rootDir
        if not self.rootDir.endswith('/'):
            self.rootDir += '/'     

        self.filt = filt
        
        if epochDirSuffix != None:
            self.suffix = '_' + epochDirSuffix
        else:
            self.suffix = ''

        if imgSuffix != None:
            self.imgSuffix = '_' + imgSuffix
        else:
            self.imgSuffix = ''

        # Setup Directories
        # Epoch directory
        self.dirEpoch = self.rootDir + self.epoch + self.suffix + '/'
        
        # Clean directory
        self.dirClean = self.dirEpoch + 'clean/' + self.filt + '/'
        if clean_dir is not None:
            self.dirClean = util.trimdir(os.path.abspath(clean_dir) + '/'
                                         + self.filt + '/')
        
        if useDistorted:
            self.dirClean += 'distort/'
        if stfDir is not None:
            self.dirCleanStf = self.dirClean + 'starfinder/' + stfDir + '/'
        else:
            self.dirCleanStf = self.dirClean + 'starfinder/'
        self.dirCleanAln = self.dirCleanStf + 'align/'
        
        # Combo directory
        self.dirCombo = self.dirEpoch + 'combo/'
        if combo_dir is not None:
            self.dirCombo = util.trimdir(os.path.abspath(combo_dir) + '/')
        
        # Starfinder directory
        starfinder_dir_name = 'starfinder/'
        if stf_version is not None:
            starfinder_dir_name = 'starfinder_{0}/'.format(stf_version)
        
        if stfDir is not None:
            self.dirComboStf = self.dirCombo + starfinder_dir_name + stfDir + '/'
        else:
            self.dirComboStf = self.dirCombo + starfinder_dir_name
        
        if combo_stf_dir is not None:
            self.dirComboStf = util.trimdir(os.path.abspath(combo_stf_dir)
                                            + '/')
            if stfDir is not None:
                self.dirComboStf = self.dirComboStf + starfinder_dir_name + stfDir + '/'
            else:
                self.dirComboStf = self.dirComboStf + starfinder_dir_name
        
        
        self.dirComboAln = self.dirComboStf + 'align/'
        
        # make the directories if we need to
        if cleanList != None:
            util.mkdir(self.dirCleanStf)
            util.mkdir(self.dirCleanAln)
        util.mkdir(self.dirComboStf)
        util.mkdir(self.dirComboAln)
        
        # Get the list of clean files to work on
        self.cleanList = cleanList
        self.cleanFiles = []
        if cleanList != None:
            _list = Table.read(self.dirClean + self.cleanList, format='ascii.no_header')
            for ii in range(len(_list)):
                fields = _list[_list.colnames[0]][ii].split('/')
                filename = fields[-1].replace('.fits', '')
                self.cleanFiles.append(filename)

        # Keep track of the current directory we started in
        self.dirStart = os.getcwd()

        # Set the default magnitude cut for plotPosError
        self.plotPosMagCut = 15.0
        
    def prepStarfinder(self, targetName, targetCoords, psfStars, filterName):
        """Creates a _psf.list file and saves it in the source_list/ directory."""
        
        # Covering possible variable inputs and initializing
        targetName = targetName.lower()
        filterName = filterName.capitalize()
        targetCoords = targetCoords[:]
        psfListLocation = self.rootDir + 'source_list/' + targetName + '_psf.list'
        year = date.today().strftime("%Y" + ".0")
        
        # Modifying the starlist provided into an astropy table
        starno = 0
        for star in psfStars:
            star[0] = "{:.3f}".format(((star[0] - targetCoords[0]) * (9.942 / 1000) * -1))
            star[1] = "{:.3f}".format(((star[1] - targetCoords[1]) * (9.942 / 1000)))
            star.insert(0, "1.000")
            if starno > 0:
                star.insert(0, "S" + str(starno).zfill(3))
            else:
                star.insert(0, targetName)
            star.insert(4, "-")
            star.insert(4, year)
            star.insert(4, "0.000")
            star.insert(4, "0.000")
            starno += 1
            
        psfStars = np.array(psfStars)
        psfStarTable = Table(psfStars, names = ('#Name', 'Kp', 'Xarc',
                                                      'Yarc', 'Vx', 'Vy',
                                                      't0', 'Filt', 'PSF?'))
        
        ascii.write(psfStarTable, psfListLocation, 
                    format = 'fixed_width', 
                    delimiter = '  ', 
                    bookend = False, 
                    overwrite = True)
    
    def analyzeCombo(self):
        self.starfinderCombo()
        self.calibrateCombo()
        self.alignCombo()



    def starfinderCombo(self, oldPsf=False):
        try:
            print('COMBO starfinder')
            print('Coo Star: ' + self.cooStar)
            
            os.chdir(self.dirComboStf)
            if self.type == 'ao':
                # Write an IDL batch file for the main map and each submap
                batch_files = []
                batch_file_logs = []
                
                combos = ['main']
                
                for cur_combo in combos:
                    batch_file = 'idlbatch_main_' + self.filt
                    batch_file_log = batch_file + '.log'
                    
                    if cur_combo != 'main':
                        batch_file = 'idlbatch_sm{0}_{1}'.format(cur_combo, self.filt)
                        batch_file_log = batch_file + '.log'
                    
                    util.rmall([batch_file, batch_file_log])
                
                    batch_out = ""
                    
                    combo_file_path = "{0}/mag{1}_{2}.fits".format(
                        self.dirCombo, self.epoch, self.filt,
                    )
                    
                    if cur_combo != 'main':
                        combo_file_path = "{0}/m{1}_{2}_{3}.fits".format(
                            self.dirCombo, self.epoch, self.filt, cur_combo,
                        )
                    
                    batch_out += "find_stf_new, "
                    batch_out += "'{0}', ".format(combo_file_path)
                    
                    if cur_combo == 'main':
                        batch_out += "{0}, ".format(self.corrMain)
                    else:
                        batch_out += "{0}, ".format(self.corrSub)
                    
                    batch_out += "ttStar='TTstar', gsStar='', "
                
                    # Add deblend flag if using it
                    if self.deblend != None:
                        batch_out += "deblend={0}, ".format(self.deblend)
                
                    if not oldPsf:
                        batch_out += "/makePsf, "
                
                    batch_out += "makeRes=1, "
                    batch_out += "makeStars=1, "
                    
                    if self.airopa_mode == 'legacy':
                        batch_out += "/legacy, "
                
                    if self.airopa_mode == 'variable':
                        batch_out += "/aoopt, "
                
                    batch_out += "cooStar='{0}', ".format(self.cooStar)
                    batch_out += "starlist='{0}', ".format(self.starlist)
                    batch_out += "trimfake={0}, ".format(self.trimfake)
                
                    if self.stf_debug:
                        batch_out += "/debug, "
                
                    batch_out += "fixPsf=1, "
                    batch_out += "backboxFWHM=25, "
                    batch_out += "flat=1, "
                    batch_out += "subtract=1"
                    
                    # Support for arbitrary starfinder flags.
                    if self.stfFlags != '': 
                        batch_out += ", {0}".format(self.stfFlags)
                
                    batch_out += "\n"
                    batch_out += "exit\n"
                    
                    # Write out the current batch file's contents
                    with open(batch_file, 'w') as out_file:
                        out_file.write(batch_out)
                    
                    batch_files.append(batch_file)
                    batch_file_logs.append(batch_file_log)
                
                if self.stf_parallelize:
                    stf_pool = Pool(processes=4)
                    
                    stf_pool.map(idl_process_run, batch_files)
                else:
                    for batch_file in batch_files:
                        idl_process_run(batch_file)
                
                # Combine all log files into one combo log file
                # (to preserve backward compatibility for code that reads log files)
                fileIDLlog = 'idlbatch_combo_' + self.filt + '.log'
                with open(fileIDLlog, 'wb') as out_file:
                    for log_file in batch_file_logs:
                        with open(log_file, 'rb') as in_file:
                            out_file.write(in_file.read())
                
            elif self.type == 'speckle':
                fileIDLbatch = 'idlbatch_combo' 
                fileIDLlog = fileIDLbatch + '.log'
                util.rmall([fileIDLlog, fileIDLbatch])
                
                _batch = open(fileIDLbatch, 'w')
                _batch.write("find_new_speck, ")
                _batch.write("'" + self.epoch + "', ")
                _batch.write("corr_main={0}, " % self.corrMain)
                _batch.write("corr_subs={0}, " % self.corrSub)
                _batch.write("starlist='" + self.starlist + "', ")
                if self.airopa_mode == 'legacy':
                    _batch.write("/legacy, ")
                # Variable mode isn't supported for Speckle data. 
                # if mode == 'variable':
                #     _batch.write("/aoopt, ")
                if not oldPsf:
                    _batch.write("/makePsf, ")
                _batch.write("rootDir='" + self.rootDir + "'")
                
                # Support for arbitrary starfinder flags.
                _batch.write(self.stfFlags)
                
                _batch.write("\n")
                _batch.write("exit\n")
                _batch.close()
                
                cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog            
                #os.system(cmd)
                subp = subprocess.Popen(cmd, shell=True, executable="/bin/tcsh")
                tmp = subp.communicate()
            
            # Write data_sources file
            data_sources_file = open(self.dirComboStf + '/data_sources.txt', 'w')
            data_sources_file.write('---\n# Starfinder run on image files\n')
            
            # Copy over the starfinder FITS files to the current directory
            epoch_file_root = 'mag{0}_{1}'.format(self.epoch, self.filt)
            
            shutil.copyfile(self.dirCombo + epoch_file_root + '_back.fits',
                            epoch_file_root + '_back.fits')
            shutil.copyfile(self.dirCombo + epoch_file_root + '_psf.fits',
                            epoch_file_root + '_psf.fits')
            shutil.copyfile(self.dirCombo + epoch_file_root + '_res.fits',
                            epoch_file_root + '_res.fits')
            shutil.copyfile(self.dirCombo + epoch_file_root + '_stars.fits',
                            epoch_file_root + '_stars.fits')
            
            out_line = '{0} from {1}{2} ({3})\n'.format(
                epoch_file_root + '.fits', self.dirCombo,
                epoch_file_root + '.fits', datetime.now())
            data_sources_file.write(out_line)
            
            
            # Close data_sources file
            data_sources_file.close()
            
            # Copy over the PSF starlist that was used (for posterity).
            outPsfs = 'mag%s%s_%s_psf_list.txt' % (self.epoch, self.imgSuffix, self.filt)
            shutil.copyfile(self.starlist, outPsfs)
            
            
            # Run the Strehl calculator to calculate Strehl on StarFinder PSFs
            # Need list of PSFs file and coo files associated with each PSF file
            psf_file_list = []
            
            cur_psf_file = epoch_file_root + '_psf.fits'
            cur_coo_file = epoch_file_root + '_psf.coo'
            
            psf_img, psf_hdr = fits.getdata(cur_psf_file, header=True)
            (psf_y_size, psf_x_size) = psf_img.shape
            
            psf_x_coord = psf_x_size/2.0
            psf_y_coord = psf_y_size/2.0
            
            coo_output = ' {0:.3f}   {1:.3f}'.format(psf_x_coord, psf_y_coord)
            
            psf_file_list.append(cur_psf_file)
            with open(cur_coo_file, 'w') as out_coo_file:
                out_coo_file.write(coo_output)
            
            for cur_sub_index in range(self.numSubMaps):
                cur_sub_str = str(cur_sub_index + 1)
                cur_psf_file = epoch_sub_root + '_' + cur_sub_str + '_psf.fits'
                cur_coo_file = epoch_sub_root + '_' + cur_sub_str + '_psf.coo'
                
                psf_file_list.append(cur_psf_file)
                with open(cur_coo_file, 'w') as out_coo_file:
                    out_coo_file.write(coo_output)
            
            strehl.calc_strehl(psf_file_list,
                               'stf_psf_strehl_{0}.txt'.format(self.filt),
                               instrument=self.instrument)
            
            os.chdir(self.dirStart)
        except:
            os.chdir(self.dirStart)
            raise
            


    def calibrateCombo(self):
        try:
            print('COMBO calibrate')
            
            # Get the position angle from the *.fits header
            # We assume that the submaps have the same PA as the main maps
            os.chdir(self.dirCombo)

            calCamera = 1
            if self.type == 'ao':
                fitsFile = 'mag%s%s_%s.fits' % (self.epoch, self.imgSuffix, self.filt)
                hdr = fits.getheader(fitsFile)
                angle = self.instrument.get_position_angle(hdr)
                print(self.instrument.name)
                
                # Check for wide camera
                calCamera = calibrate.get_camera_type(fitsFile)
            elif self.type == 'speckle':
                angle = 0.0
                fitsFile = 'mag%s.fits' % self.epoch

            
            # CALIBRATE
            os.chdir(self.dirComboStf)

            cmd = 'calibrate_new %s' % self.calFlags
            cmd += '-T %.1f ' % angle
            cmd += '-I %s ' % self.calCooStar
            cmd += '-N %s ' % self.calFile
            cmd += '-M %s ' % self.calColumn
            cmd += '-c %d ' % calCamera
            cmd += '-V '
            if (self.calStars != None) and (len(self.calStars) > 0):
                cmd += '-S '
                for cc in range(len(self.calStars)):
                    if (cc != 0):
                        cmd += ','
                    cmd += self.calStars[cc]
            cmd += ' '

            # Calibrate Main Map
            if self.type == 'speckle':   # This stuff is left over.
                fileMain = 'mag%s_%3.1f_stf.lis' % \
                (self.epoch, self.corrMain)
            else:
                if self.deblend == 1:
                    fileMain = 'mag%s%s_%s_%3.1fd_stf.lis' % \
                    (self.epoch, self.imgSuffix, self.filt, self.corrMain[0])
                else:
                    fileMain = 'mag%s%s_%s_%3.1f_stf.lis' % \
                    (self.epoch, self.imgSuffix, self.filt, self.corrMain[0])
            print(cmd + fileMain)

            # Now call from within python... don't bother with command line anymore.
            argsTmp = cmd + fileMain
            args = argsTmp.split()[1:]
            calibrate.main(args)

            # Copy over the calibration list.
            shutil.copyfile(self.calFile, fileMain.replace('.lis', '_cal_phot_list.txt'))

            os.chdir(self.dirStart)
        except:
            os.chdir(self.dirStart)
            raise


    

