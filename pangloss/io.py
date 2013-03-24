# ===========================================================================

import cPickle

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

# ======================================================================
