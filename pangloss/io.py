# ===========================================================================

import cPickle

# ======================================================================

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
