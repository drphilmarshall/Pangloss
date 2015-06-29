# ===========================================================================

import pangloss
import os,cPickle,astropy.io

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

# PJM: we need to switch to astropy tables...
# Here's how Richard McMahon uses them, admittedly when reading in FITS:

    '''
    here is an example; I like to write out the version number to help with debugging.
    
    import astropy
    print('astropy: ', astropy.version)
    from astropy.table import Table
    from astropy.io import fits
    
    infile_OM10="/home/rgm/soft/OM10/OM10/data/qso_mock.fits"
    print 'Read in with Astropy FITS table reader'
    table=Table.read(infile_OM10)
    table.pprint()
    print 'Numer of rows: ', len(table)
    print 'columns: ', table.columns
    
    print 'colnames: ', table.colnames
    print 'meta: ', table.meta
    
    # see http://astropy.readthedocs.org/en/latest/table/modify_table.html
    # to see how to add a column, e.g.
    # to insert before the first table column, do:
    
    table.add_column(aa, index=0)
    
    table.write('new.fits')
    
    '''
    table = astropy.io.Table(filename, type='ascii')

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
        mag = table[config.parameters['MagName']]
        table.add_column('mag',mag)
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
