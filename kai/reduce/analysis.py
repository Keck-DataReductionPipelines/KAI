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
                 airopa_mode='single', stf_version=None,
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
        """
        # Setup default parameters
        self.type = 'ao'
        self.corrMain = 0.8
        self.corrSub = 0.6
        self.corrClean = 0.7
        # airopa mode can be "legacy", "single", or "variable"
        self.airopa_mode = airopa_mode
        self.trimfake = 1
        self.stfFlags = ''

        self.starlist = rootDir + 'source_list/psf_central.dat'
        self.labellist = rootDir+ 'source_list/label.dat'
        self.orbitlist = rootDir+ 'source_list/orbits.dat'
        self.calFile = rootDir + 'source_list/photo_calib.dat'

        # Keep track of the instrument these images are from.
        # We need this to get things like plate scale, etc.
        self.instrument = instrument
        
        self.calStars = ['irs16NW', 'S3-22', 'S1-17', 'S1-34', 'S4-3', 'S1-1',
                         'S1-21', 'S3-370', 'S3-88', 'S3-36', 'S2-63']
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

        self.alignFlags = '-R 3 -v -p -a 0'

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

    def analyzeClean(self):
        self.starfinderClean()
        self.calibrateClean()
        self.alignClean()

    def analyzeComboClean():
        self.starfinderCombo()
        self.starfinderClean()
        self.calibrateCombo()
        self.calibrateClean()
        self.alignCombo()
        self.alignClean()


    def starfinderCombo(self, oldPsf=False):
        try:
            print('COMBO starfinder')
            print('Coo Star: ' + self.cooStar)
            
            os.chdir(self.dirComboStf)
            if self.type == 'ao':
                # Write an IDL batch file
                fileIDLbatch = 'idlbatch_combo_' + self.filt
                fileIDLlog = fileIDLbatch + '.log'
                util.rmall([fileIDLlog, fileIDLbatch])

                _batch = open(fileIDLbatch, 'w')
                _batch.write("find_stf_new, ")
                _batch.write("'" + self.epoch + "', ")
                _batch.write("'" + self.filt + "', ")
                _batch.write("corr_main=%3.1f, " % self.corrMain)
                _batch.write("corr_subs=%3.1f, " % self.corrSub)
                # Only send in the deblend flag if using it!
                if self.deblend != None:
                    _batch.write("deblend=" + str(self.deblend) + ", ")
                _batch.write("cooStar='" + self.cooStar + "', ")
                _batch.write("suffixEpoch='" + self.suffix + "', ")
                _batch.write("imgSuffix='" + self.imgSuffix + "', ")
                _batch.write("starlist='" + self.starlist + "', ")
                if self.airopa_mode == 'legacy':
                    _batch.write("/legacy, ")
                if self.airopa_mode == 'variable':
                    _batch.write("/aoopt, ")
                _batch.write("trimfake=" + str(self.trimfake) + ", ")
                if not oldPsf:
                    _batch.write("/makePsf, ")
                
                _batch.write("rootDir='" + self.rootDir + "'")

                # Support for arbitrary starfinder flags.
                _batch.write(self.stfFlags)
                
                _batch.write("\n")
                _batch.write("exit\n")
                _batch.close()
            elif self.type == 'speckle':
                fileIDLbatch = 'idlbatch_combo' 
                fileIDLlog = fileIDLbatch + '.log'
                util.rmall([fileIDLlog, fileIDLbatch])
                
                _batch = open(fileIDLbatch, 'w')
                _batch.write("find_new_speck, ")
                _batch.write("'" + self.epoch + "', ")
                _batch.write("corr_main=%3.1f, " % self.corrMain)
                _batch.write("corr_subs=%3.1f, " % self.corrSub)
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
            
            # Copy over submap starfinder FITS files to the current directory
            epoch_sub_root = 'm{0}_{1}'.format(self.epoch, self.filt)
            
            for cur_sub_index in range(self.numSubMaps):
                cur_sub_str = str(cur_sub_index + 1)
                            
                shutil.copyfile(self.dirCombo + epoch_sub_root + '_' + cur_sub_str + '_back.fits',
                                epoch_sub_root + '_' + cur_sub_str + '_back.fits')
                shutil.copyfile(self.dirCombo + epoch_sub_root + '_' + cur_sub_str + '_psf.fits',
                                epoch_sub_root + '_' + cur_sub_str + '_psf.fits')
                shutil.copyfile(self.dirCombo + epoch_sub_root + '_' + cur_sub_str + '_res.fits',
                                epoch_sub_root + '_' + cur_sub_str + '_res.fits')
                shutil.copyfile(self.dirCombo + epoch_sub_root + '_' + cur_sub_str + '_stars.fits',
                                epoch_sub_root + '_' + cur_sub_str + '_stars.fits')
                
                out_line = '{0} from {1}{2} ({3})\n'.format(
                    epoch_sub_root + '_' + cur_sub_str + '.fits', self.dirCombo,
                    epoch_sub_root + '_' + cur_sub_str + '.fits',
                    datetime.now())
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
            
        
    def starfinderClean(self):
        try:
            print('CLEAN starfinder')
            
            os.chdir(self.dirCleanStf)
            
            # Write an IDL batch file
            fileIDLbatch = 'idlbatch_clean_' + self.filt
            fileIDLlog = fileIDLbatch + '.log'
            util.rmall([fileIDLlog, fileIDLbatch])
            
            _batch = open(fileIDLbatch, 'w')
            _batch.write("find_stf_clean, ")
            _batch.write("'" + self.epoch + "', ")
            _batch.write("'" + self.filt + "', ")
            _batch.write("corr=" + str(self.corrClean) + ", ")
            _batch.write("cooStar='" + self.cooStar + "', ")
            _batch.write("suffixEpoch='" + self.suffix + "', ")
            _batch.write("starlist='" + self.starlist + "', ")
            _batch.write("fileList='" + self.cleanList + "', ")
            if self.airopa_mode == 'legacy':
                _batch.write("/legacy, ")
            if self.airopa_mode == 'variable':
                _batch.write("/aoopt, ")
            _batch.write("trimfake=" + str(self.trimfake) + ", ")
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
            
            # Copy over the PSF starlist that was used (for posterity).
            outPsfs = 'c_%s_%s_psf_list.txt' % (self.epoch, self.filt)
            shutil.copyfile(self.starlist, outPsfs)
            
            os.chdir(self.dirStart)
        except:
            os.chdir(self.dirStart)
            raise

    def starfinderCleanLoop(self):
        try:
            print('CLEAN starfinder loop')
            os.chdir(self.dirCleanStf)

            for i, filename in enumerate(self.cleanFiles):
                # Write an IDL batch file
                fileIDLbatch = filename + '_idlbatch'
                fileIDLlog = fileIDLbatch + '.log'
                util.rmall([fileIDLlog, fileIDLbatch])

                _batch = open(fileIDLbatch, 'w')
                _batch.write("find_stf, ")
                _batch.write("'" + self.dirClean + filename + ".fits', ")
                _batch.write("[" + str(self.corrClean) + "], ")
                if self.airopa_mode =='single':
                    _batch.write("aoopt=0, ")
                elif self.airopa_mode == 'legacy':
                    _batch.write("/legacy, ")
                elif self.airopa_mode == 'variable':
                    _batch.write("/aoopt, ")

                # Support for arbitrary starfinder flags.
                _batch.write(self.stfFlags)

                _batch.write("starlist='" + self.starlist + "', ")
                _batch.write("cooStar='" + self.cooStar + "'")          #no trailing comma for final parameter

                _batch.write("\n")
                _batch.write("exit\n")
                _batch.close()
                
                cmd = 'idl < ' + fileIDLbatch + ' >& ' + fileIDLlog
                #os.system(cmd)
                print('Running', cmd)
                subp = subprocess.Popen(cmd, shell=True, executable="/bin/tcsh")
                tmp = subp.communicate()
            
            # Copy over the PSF starlist that was used (for posterity).
            outPsfs = 'c_%s_%s_psf_list.txt' % (self.epoch, self.filt)
            shutil.copyfile(self.starlist, outPsfs)
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
                    (self.epoch, self.imgSuffix, self.filt, self.corrMain)
                else:
                    fileMain = 'mag%s%s_%s_%3.1f_stf.lis' % \
                    (self.epoch, self.imgSuffix, self.filt, self.corrMain)
            print(cmd + fileMain)

            # Now call from within python... don't bother with command line anymore.
            argsTmp = cmd + fileMain
            args = argsTmp.split()[1:]
            calibrate.main(args)

            # Copy over the calibration list.
            shutil.copyfile(self.calFile, fileMain.replace('.lis', '_cal_phot_list.txt'))
            
            # Calibrate Sub Maps
            for ss in range(self.numSubMaps):
                if self.type == 'speckle':
                    fileSub = 'm%s_%d_%3.1f_stf.lis' % \
                        (self.epoch, ss+1, self.corrSub)
                else:
                    if self.deblend == 1:
                        fileSub = 'm%s%s_%s_%d_%3.1fd_stf.lis' % \
                        (self.epoch, self.imgSuffix, self.filt, ss+1, self.corrSub)
                    else:
                        fileSub = 'm%s%s_%s_%d_%3.1f_stf.lis' % \
                        (self.epoch, self.imgSuffix, self.filt, ss+1, self.corrSub)

                print(cmd + fileSub)
                
                argsTmp = cmd + fileSub
                args = argsTmp.split()[1:]
                calibrate.main(args)

            os.chdir(self.dirStart)
        except:
            os.chdir(self.dirStart)
            raise

    def calibrateClean(self):
        try:
            # Calibrate each file
            os.chdir(self.dirCleanStf)
            print("DEBUG 2nd dec --",self.dirCleanStf)
            # open writeable file for log
            _log = open('calibrate.log','w')


            cmdBase = 'calibrate_new %s ' % self.calFlags
            cmdBase += '-I %s ' % self.calCooStar
            cmdBase += '-N %s ' % self.calFile
            cmdBase += '-M %s ' % self.calColumn
            if (self.calStars != None) and (len(self.calStars) > 0):
                cmdBase += '-S '
                for cc in range(len(self.calStars)):
                    if (cc != 0):
                        cmdBase += ','
                    cmdBase += self.calStars[cc]
            cmdBase += ' '

            for file in self.cleanFiles:
                fitsFile = '../%s.fits' % file
                listFile = '%s_%3.1f_stf.lis' % (file, self.corrClean)

                # Get the position angle
                angle = float(fits.getval(fitsFile, 'ROTPOSN')) - 0.7
                calCamera = calibrate.get_camera_type(fitsFile)

                cmd = cmdBase + ('-T %.1f -c %d ' % (angle, calCamera))

                cmd += listFile

                # 2nd dec 2009 - changed "_out" "_log"
                _log.write(cmd + '\n')

                args = cmd.split()[1:]
                calibrate.main(args)

            _log.close()  # clean up

            # Copy over the calibration list.
            shutil.copyfile(self.calFile, 'clean_phot_list.txt')
            
            os.chdir(self.dirStart)
        except:
            os.chdir(self.dirStart)
            raise

    def alignCombo(self):
        print('ALIGN_RMS combo')

        if self.type == 'ao':
            file_ext = '_' + self.filt
        else:
            file_ext = self.filt

        try:
            os.chdir(self.dirCombo)
            # Get the align data type
            if self.type == 'ao':
                fitsFile = 'mag%s%s_%s.fits' % (self.epoch, self.imgSuffix, self.filt)
            elif self.type == 'speckle':
                fitsFile = 'mag%s.fits' % self.epoch
            hdr = fits.getheader(fitsFile)
            alignType = self.instrument.get_align_type(hdr, errors=False)


            os.chdir(self.dirComboAln)

            # Put the files in to the align*.list file
            alnList1 = 'align%s%s_%3.1f.list' % (self.imgSuffix, file_ext, self.corrMain)
            alnList2 = 'align%s%s_%3.1f_named.list' % (self.imgSuffix, file_ext, self.corrMain)


            _list = open(alnList1, 'w')
            if self.deblend == 1:
                _list.write('../mag%s%s%s_%3.1fd_stf_cal.lis %d ref\n' %
                            (self.epoch, self.imgSuffix, file_ext, self.corrMain, alignType))
            else:
                _list.write('../mag%s%s%s_%3.1f_stf_cal.lis %d ref\n' %
                            (self.epoch, self.imgSuffix, file_ext, self.corrMain, alignType))
            for ss in range(self.numSubMaps):
                if self.deblend == 1:
                    _list.write('../m%s%s%s_%d_%3.1fd_stf_cal.lis %d\n' %
                                (self.epoch, self.imgSuffix, file_ext, ss+1, self.corrSub, alignType))
                else:
                    _list.write('../m%s%s%s_%d_%3.1f_stf_cal.lis %d\n' %
                                (self.epoch, self.imgSuffix, file_ext, ss+1, self.corrSub, alignType))

                    
            _list.close()

            shutil.copyfile(alnList1, alnList2)

            # Make an unlabeled version
            cmd = 'java -Xmx1024m align %s ' % (self.alignFlags)
            cmd += '-r align%s%s_%3.1f ' % (self.imgSuffix, file_ext, self.corrMain)
            cmd += alnList1
            print(cmd)
            #os.system(cmd)
            subp = subprocess.Popen(cmd, shell=True, executable="/bin/tcsh")
            tmp = subp.communicate()

            # Make a named/labeled version
            cmd = 'java -Xmx1024m align %s ' % (self.alignFlags)
            cmd += '-N %s ' % self.labellist
            if (self.orbitlist != None) and (self.orbitlist != ''):
                cmd += '-o %s ' % self.orbitlist
            cmd += '-r align%s%s_%3.1f_named ' % (self.imgSuffix, file_ext, self.corrMain)
            cmd += alnList2
            print(cmd)

            subp = subprocess.Popen(cmd, shell=True, executable="/bin/tcsh")
            tmp = subp.communicate()


            align_options = 'align%s%s_%3.1f %d -e' % \
              (self.imgSuffix, file_ext, self.corrMain, self.minSubMaps)
            align_rms.run(align_options.split())

            align_options = 'align%s%s_%3.1f_named %d -e' % \
              (self.imgSuffix, file_ext, self.corrMain, self.minSubMaps)
            align_rms.run(align_options.split())


            # Move the resulting files to their final resting place
            os.rename('align%s%s_%3.1f_rms.lis' % 
                      (self.imgSuffix, file_ext, self.corrMain),
                      '../mag%s%s%s_rms.lis' % 
                      (self.epoch, self.imgSuffix, file_ext))
            os.rename('align%s%s_%3.1f_named_rms.lis' % 
                      (self.imgSuffix, file_ext, self.corrMain),
                      '../mag%s%s%s_rms_named.lis' % 
                      (self.epoch, self.imgSuffix, file_ext))

            # Copy over the label.dat and orbit.dat file that was used.
            shutil.copyfile(self.labellist,
                            'align%s%s_%3.1f_named_label_list.txt' % 
                            (self.imgSuffix, file_ext, self.corrMain))
                            
            if (self.orbitlist != None) and (self.orbitlist != ''):
                shutil.copyfile(self.orbitlist,
                                'align%s%s_%3.1f_named_orbit_list.txt' % 
                                (self.imgSuffix, file_ext, self.corrMain))

            # Now plot up the results
            plotSuffix = self.imgSuffix + file_ext
            os.chdir(self.dirComboStf)
            plotPosError('mag%s%s%s_rms.lis' % (self.epoch, self.imgSuffix, file_ext),
                         raw=True, suffix=plotSuffix, magCutOff=self.plotPosMagCut)

            os.chdir(self.dirStart)
        except:
            os.chdir(self.dirStart)
            raise

    
    def alignClean(self):
        try:
            os.chdir(self.dirClean)

            # Get the align data type
            fitsFile = self.cleanFiles[0] + '.fits'
            hdr = fits.getheader(fitsFile)
            alignType = instrument.get_align_type(hdr, errors=False)
            alignCombo = alignType + 1

            # Make the align*.list file
            os.chdir(self.dirCleanAln)
            alnList = 'align_%s%s_%3.1f.list' % (self.imgSuffix, self.filt, self.corrClean)
            _list = open(alnList, 'w')

            # copy main map into it
            _list.write('%smag%s%s_%s_rms.lis %d ref\n' % 
                        (self.dirComboStf, self.epoch, self.imgSuffix, self.filt, alignCombo))
            # Add each clean file
            for ff in self.cleanFiles:
                _list.write('../%s_%3.1f_stf_cal.lis %d\n' %
                            (ff, self.corrClean, alignType))
            _list.close()

            # Run align
            cmd = 'java -Xmx1024m align %s ' % (self.alignFlags)
            cmd += '-N %s ' % self.labellist
            if (self.orbitlist != None) and (self.orbitlist != ''):
                cmd += '-o %s ' % self.orbitlist
            cmd += '-r align_%s%s_%3.1f_named ' % (self.imgSuffix, self.filt, self.corrClean)
            cmd += alnList

            #os.system(cmd)
            subp = subprocess.Popen(cmd, shell=True, executable="/bin/tcsh")
            tmp = subp.communicate()


            # Copy over the label.dat file that was used.
            shutil.copyfile(self.labellist,
                            'align_%s%s_%3.1f_named_label_list.txt' % 
                            (self.imgSuffix, self.filt, self.corrClean))
            if (self.orbitlist != None) and (self.orbitlist != ''):
                shutil.copyfile(self.orbitlist,
                                'align_%s%s_%3.1f_named_orbit_list.txt' % 
                                (self.imgSuffix, self.filt, self.corrClean))
                            
            os.chdir(self.dirStart)
        except:
            os.chdir(self.dirStart)
            raise


def plotPosError(starlist, raw=False, suffix='', radius=4, magCutOff=15.0,
                 title=True):
    """
    Make three standard figures that show the data quality 
    from  a *_rms.lis file. 

    1. astrometric error as a function of magnitude.
    2. photometric error as a function of magnitude.
    3. histogram of number of stars vs. magnitude.

    Use raw=True to plot the individual stars in plots 1 and 2.
    """
    # Load up the starlist
    lis = Table.read(starlist, format='ascii')

    # Assume this is NIRC2 data.
    scale = 0.00995
    
    name = lis[lis.colnames[0]]
    mag = lis[lis.colnames[1]]
    x = lis[lis.colnames[3]]
    y = lis[lis.colnames[4]]
    xerr = lis[lis.colnames[5]]
    yerr = lis[lis.colnames[6]]
    snr = lis[lis.colnames[7]]
    corr = lis[lis.colnames[8]]

    merr = 1.086 / snr

    # Convert into arsec offset from field center
    # We determine the field center by assuming that stars
    # are detected all the way out the edge.
    xhalf = x.max() / 2.0
    yhalf = y.max() / 2.0
    x = (x - xhalf) * scale
    y = (y - yhalf) * scale
    xerr *= scale * 1000.0
    yerr *= scale * 1000.0

    r = np.hypot(x, y)
    err = (xerr + yerr) / 2.0

    magStep = 1.0
    radStep = 1.0
    magBins = np.arange(10.0, 20.0, magStep)
    radBins = np.arange(0.5, 9.5, radStep)
    
    errMag = np.zeros(len(magBins), float)
    errRad = np.zeros(len(radBins), float)
    merrMag = np.zeros(len(magBins), float)
    merrRad = np.zeros(len(radBins), float)

    ##########
    # Compute errors in magnitude bins
    ########## 
    #print '%4s  %s' % ('Mag', 'Err (mas)')
    for mm in range(len(magBins)):
        mMin = magBins[mm] - (magStep / 2.0)
        mMax = magBins[mm] + (magStep / 2.0)
        idx = (np.where((mag >= mMin) & (mag < mMax) & (r < radius)))[0]

        if (len(idx) > 0):
            errMag[mm] = np.median(err[idx])
            merrMag[mm] = np.median(merr[idx])
        
        #print '%4.1f  %5.2f' % (magBins[mm], errMag[mm])
        
                       
    ##########
    # Compute errors in radius bins
    ########## 
    for rr in range(len(radBins)):
        rMin = radBins[rr] - (radStep / 2.0)
        rMax = radBins[rr] + (radStep / 2.0)
        idx = (np.where((r >= rMin) & (r < rMax) & (mag < magCutOff)))[0]

        if (len(idx) > 0):
            errRad[rr] = np.median(err[idx])
            merrRad[rr] = np.median(err[idx])

    idx = (np.where((mag < magCutOff) & (r < radius)))[0]
    errMedian = np.median(err[idx])
    numInMedian = len(idx)

    ##########
    #
    # Plot astrometry errors
    #
    ##########
 
    # Remove figures if they exist -- have to do this
    # b/c sometimes the file won't be overwritten and
    # the program crashes saying 'Permission denied'
    if os.path.exists('plotPosError%s.png' % suffix):
        os.remove('plotPosError%s.png' % suffix)
    if os.path.exists('plotMagError%s.png' % suffix):
        os.remove('plotMagError%s.png' % suffix)
    if os.path.exists('plotNumStars%s.png' % suffix):
        os.remove('plotNumStars%s.png' % suffix)

    if os.path.exists('plotPosError%s.eps' % suffix):
        os.remove('plotPosError%s.eps' % suffix)
    if os.path.exists('plotMagError%s.eps' % suffix):
        os.remove('plotMagError%s.eps' % suffix)
    if os.path.exists('plotNumStars%s.eps' % suffix):
        os.remove('plotNumStars%s.eps' % suffix)

    py.figure(figsize=(6,6))
    py.clf()
    py.subplots_adjust(left=0.15, top=0.92, right=0.95)
    if (raw == True):
        idx = (np.where(r < radius))[0]
        py.semilogy(mag[idx], err[idx], 'k.')
        
    py.semilogy(magBins, errMag, 'g.-', lw=2)
    py.axis([8, 22, 1e-2, 30.0])
    py.xlabel('K Magnitude for r < %4.1f"' % radius, fontsize=16)
    py.ylabel('Positional Uncertainty (mas)', fontsize=16)
    if title == True:
        py.title(starlist)
    
    py.savefig('plotPosError%s.png' % suffix)
    #py.savefig('plotPosError%s.eps' % suffix)

    ##########
    #
    # Plot photometry errors
    #
    ##########
    py.clf()
    if (raw == True):
        idx = (np.where(r < radius))[0]
        py.plot(mag[idx], merr[idx], 'k.')
        
    py.plot(magBins, merrMag, 'g.-')
    py.axis([8, 22, 0, 0.15])
    py.xlabel('K Magnitude for r < %4.1f"' % radius)
    py.ylabel('Photo. Uncertainty (mag)')
    py.title(starlist)
    
    py.savefig('plotMagError%s.png' % suffix)
    #py.savefig('plotMagError%s.eps' % suffix)

    ##########
    # 
    # Plot histogram of number of stars detected
    #
    ##########
    py.clf()
    idx = (np.where(r < radius))[0]
    (n, bb, pp) = py.hist(mag[idx], bins=np.arange(9, 22, 0.5))
    py.xlabel('K Magnitude for r < %4.1f"' % radius)
    py.ylabel('Number of Stars')

    py.savefig('plotNumStars%s.png' % suffix)
    #py.savefig('plotNumStars%s.eps' % suffix)

    # Find the peak of the distribution
    maxHist = n.argmax()
    maxBin = bb[maxHist]


    ##########
    # 
    # Save relevant numbers to an output file.
    #
    ##########
    # Print out some summary information
    print('Number of detections: %4d' % len(mag))
    print('Median Pos Error (mas) for K < %2i, r < %4.1f (N=%i):  %5.2f' % \
          (magCutOff, radius, numInMedian, errMedian))
    print('Median Mag Error (mag) for K < %2i, r < %4.1f (N=%i):  %5.2f' % \
          (magCutOff, radius, numInMedian, np.median(merr[idx])))
    print('Turnover mag = %4.1f' % (maxBin))


    out = open('plotPosError%s.txt' % suffix, 'w')
    out.write('Number of detections: %4d\n' % len(mag))
    out.write('Median Pos Error (mas) for K < %2i, r < %4.1f (N=%i):  %5.2f\n' % \
          (magCutOff, radius, numInMedian, errMedian))
    out.write('Median Mag Error (mag) for K < %2i, r < %4.1f (N=%i:  %5.2f\n' % \
          (magCutOff, radius, numInMedian, np.median(merr[idx])))
    out.write('Turnover mag = %4.1f\n' % (maxBin))
    out.close()
    


       

