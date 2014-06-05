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
arcmin2deg = (1.0/60.0)
rad2arcmin = 1.0/arcmin2rad
rad2deg = 180.0/numpy.pi
deg2rad = 1.0/rad2deg

squaredeg = (numpy.pi/180.0)**2./(4*numpy.pi) # solid angle of cone subtending 1 square degree

dashedline = '--------------------------------------------------------------------------------'
doubledashedline = '================================================================================'

hello = '              Pangloss: reconstructing all the mass in the Universe             '
