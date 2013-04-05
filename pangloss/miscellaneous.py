# ===========================================================================

import numpy

# ======================================================================

"""
    NAME
        miscellaneous.py

    PURPOSE
        Bits and pieces to help make a nice output log.
        
    COMMENTS
            
    FUNCTIONS
        
    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Marshall & Collett (Oxford)
"""

#=========================================================================

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

dashedline = '--------------------------------------------------------------------------------'
doubledashedline = '================================================================================'

hello = '              Pangloss: reconstructing all the mass in the Universe             '
