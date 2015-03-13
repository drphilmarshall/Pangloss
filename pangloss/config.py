# ===========================================================================

import os, glob

# ======================================================================

class Configuration(object):
    """
    NAME
        Configuration

    PURPOSE
        Structure of metadata supplied by the Pangloss user to define 
        the experiment they want to do.

    COMMENTS

    INITIALISATION
        configfile    Plain text file containing desired values of pars
    
    METHODS
        read(self): from configfile, get par names and values
        
        convert(self): strings to numbers, and expand paths
        
        prepare(self): set up workspace
        
        getLightconePickleName(self,flavor,pointing=None): 

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Collett & Marshall (Cambridge)
    """

    def __init__(self,configfile):
        self.file = configfile
        self.parameters = {}
        self.read()
        self.convert()
        self.prepare()
        return

    # ------------------------------------------------------------------
    # Read in values by keyword and populate the parameters dictionary.

    def read(self):
        thisfile = open(self.file)
        for line in thisfile:
            # Ignore empty lines and comments:
            if line[0:2] == '\n': continue
            if line[0] == '#': continue
            line.strip()
            line = line.split('#')[0]
            # Remove whitespace and interpret Name:Value pairs:
            line = ''.join(line.split())
            line = line.split(':')
            Name, Value = line[0], line[1]
            self.parameters[Name] = Value
        thisfile.close()
        return

    # ------------------------------------------------------------------
    # Convert string values into floats/ints where necessary, and expand
    # environment variables.

    def convert(self):

        # Some values need to be floats or integers:
        for key in self.parameters.keys():
            try:
                self.parameters[key] = float(self.parameters[key])
            except ValueError:
                pass
        intkeys = ['NCalibrationLightcones','NRealisations']
        for key in intkeys:
            self.parameters[key] = int(self.parameters[key])

        # Now sort out filenames etc:

        pathkeys = ['CalibrationCatalogs', 'CalibrationKappamaps',
                    'ObservedCatalog', 'CalibrationFolder', 'HMFfile']
        for key in pathkeys:
            paths = self.parameters[key]
            # Expand environment variables (eg $PANGLOSS_DIR)
            paths = os.path.expandvars(paths)
            # Expand wildcards - glob returns [] if no files found...
            found = glob.glob(paths)
            if len(found) > 0: paths = found
            # Make sure all paths are lists, for consistency:
            if len(paths[0]) == 1: paths = [paths]
            # Replace parameters:
            self.parameters[key] = paths

        # Calibration catalogs and kappa maps must come in pairs...
        if self.parameters['CalibrationKappamaps'] != None:
            assert len(self.parameters['CalibrationCatalogs']) == \
                  len(self.parameters['CalibrationKappamaps'])


        surveycoveragekeys=['PhotometricRadius','PhotometricDepth','SpectroscopicDepth','SpectroscopicRadius']
        for key in surveycoveragekeys:
            self.parameters[key]=self.parameters[key]\
                .split('[')[1].split(']')[0].strip().split(',')
            if self.parameters[key]!=['']:
                for i in range(len(self.parameters[key])):
                    self.parameters[key][i]=float(self.parameters[key][i])
            
        return

    # ------------------------------------------------------------------
    # Perform various other preparations.

    def prepare(self):

        # Make directories if necessary:
        folderkeys = ['CalibrationFolder']
        for key in folderkeys:
            folder = self.parameters[key]
            fail = os.system('mkdir -p '+folder[0])

        return

    # ------------------------------------------------------------------
    # Figure out pickle names:

    def getLightconePickleName(self,flavor,pointing=None):

        if flavor == 'real':
            # In this case, need the name of the obscat:
            x = self.parameters['ObservedCatalog'][0]
            return x.split('.')[0]+"_lightcone.pickle"
        
        elif flavor == 'simulated':
            # In this case, need the CALIB_DIR and pointing number:
            assert pointing != None
            CALIB_DIR = self.parameters['CalibrationFolder'][0]          
            x = "%s/pointing_%i" % (CALIB_DIR, pointing)
            return x+"_lightcone.pickle"

        elif flavor == 'simulated_borg':
            # In this case, need the CALIB_DIR and pointing number:
            assert pointing != None
            CALIB_DIR = self.parameters['CalibrationFolder'][0]

            # This is for multiple catalogs to help separate them
            EXP_NAME = self.parameters['ExperimentName']            
            x = "%s/%s_pointing_%i" % (CALIB_DIR, EXP_NAME, pointing)

            return x+"_lightcone.pickle"

        return


# ======================================================================

