# ===========================================================================

import pangloss
import os,cPickle
from astropy.table import Table, Column

# ======================================================================

"""
    NAME
        io

    PURPOSE
        Useful general functions to streamline file input and output.

    COMMENTS

    FUNCTIONS
        writePickle(contents,filename):

        readPickle(filename): returns contents of pickle

        readCatalog(filename,config): returns table, given column names
                                      in configuration config

        rm(filename): silent file removal

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Marshall & Collett (Oxford)
"""

#=========================================================================

def writePickle(contents,filename):
    F = open(filename,"wb")
    cPickle.dump(contents,F,protocol=2)
    F.close()
    return

def readPickle(filename):
    F = open(filename,"rb")
    contents = cPickle.load(F)
    F.close()
    return contents

# ----------------------------------------------------------------------------

def readCatalog(filename,config):
    table = Table.read(filename, format = 'ascii')

    try: table.rename_column(config.parameters['nRAName'],'nRA')
    except: pass
    try: table.rename_column(config.parameters['DecName'],'Dec')
    except: pass

    # Calibration catalogs:
    try: table.rename_column(config.parameters['CalibMhaloName'],'Mhalo_obs')
    except: pass
    try: table.rename_column(config.parameters['CalibRedshiftName'],'z_obs')
    except: pass

    # Observed catalog:
    try: table.rename_column(config.parameters['ObsMstarName'],'Mstar_obs')
    except: pass
    try: table.rename_column(config.parameters['ObsRedshiftName'],'z_obs')
    except: pass
    try:
        table['mag'] = table[config.parameters['MagName']]
    except:
        raise "Error in io.readCatalog: no mag column called "+config.parameters['MagName']

    return table

# ----------------------------------------------------------------------------
# Remove file, if it exists, stay quiet otherwise:

def rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    return

# ======================================================================
