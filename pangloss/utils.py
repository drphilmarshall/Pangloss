import pandas as pd
from collections import OrderedDict
import os,cPickle

"""
    NAME
        utils

    PURPOSE
        Useful general functions to streamline file input and output.

    COMMENTS

    FUNCTIONS
        writePickle(contents,filename):

        readPickle(filename): returns contents of pickle

        read_hilbert_catalog(filename,config): returns table, given column names
                                      in configuration config for a Hilbert catalog

        rm(filename): silent file removal

        number(s): extract integer or floating point from string

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Marshall & Collett (Oxford)
"""

def writePickle(contents,filename):
    F = open(filename, "wb")
    cPickle.dump(contents,F,protocol=2)
    F.close()
    return

def readPickle(filename):
    F = open(filename, "rb")
    contents = cPickle.load(F)
    F.close()
    return contents

def read_hilbert_catalog(filename,config):

    try: table = pd.read_table(filename)
    except:
        raise IOError("Cannot open %s\n" % filename)

    rename = OrderedDict()
    rename[config.parameters['nRAName']] = 'nRA'
    rename[config.parameters['DecName']] = 'Dec'
    rename[config.parameters['CalibMhaloName']] = 'Mhalo_obs'
    rename[config.parameters['CalibRedshiftName']] = 'z_obs'
    rename[config.parameters['ObsMstarName']] = 'Mstar_obs'
    rename[config.parameters['ObsRedshiftName']] = 'z_obs'

    for key in rename:
        try: table = table.rename(columns={key: rename[key]})
        except: pass

    # Add a RA column:
    try: table['RA'] = -table['nRA']
    except: pass
    try:
        table['mag'] = table[config.parameters['MagName']]
    except:
        raise NameError("Error in io.readCatalog: no mag column called %s\n" % config.parameters['MagName'])

    table['GalID'] = table['GalID'].astype(str)
    table = pd.DataFrame.copy(table[table['Mhalo_obs'] > 0])
    return table

# Remove file, if it exists, stay quiet otherwise:

def rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    return

def int_or_float(s):
    try:
        return int(s)
    except ValueError:
        return float(s)
