# ===========================================================================

import pangloss

import os,cPickle,atpy

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
      Please cite: Collett et al 2013, arxiv/###

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

    table = atpy.Table(filename, type='ascii')

    try: table.rename_column(config.parameters['nRAName'],'nRA')
    except: pass
    try: table.rename_column(config.parameters['DecName'],'Dec')
    except: pass

    try: table.rename_column(config.parameters['MhaloName'],'Mhalo')
    except: pass
    try: table.rename_column(config.parameters['RedshiftName'],'z_spec')
    except: pass

    try: table.rename_column(config.parameters['MstarName'],'Mstar')
    except: pass
    try: table.rename_column(config.parameters['ObsRedshiftName'],'z_obs')
    except: pass
    try: table.rename_column(config.parameters['MagName'],'mag')
    except: pass
   
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
