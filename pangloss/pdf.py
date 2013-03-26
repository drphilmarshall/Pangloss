# ===========================================================================

import numpy

# ============================================================================

class PDF(object):
    """
    NAME
        PDF

    PURPOSE
        Define a probability distribution for some named parameters, and
        store samples drawn from it.

    COMMENTS
        The function itself is defined elsewhere - this class is just a 
        data structure.

    INITIALISATION
        parameters     List of parameter names 
        
    METHODS
        append(self,sample): add a sample (or numpy array of samples) to the ensemble
    
    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, arxiv/###

    HISTORY
      2013-03-21  Marshall & Collett (Oxford)
    """

# ----------------------------------------------------------------------------

    def __init__(self,parameters):
        
        self.name = 'Probability Density Function'
        if type(parameters) != list: parameters = [parameters]
        self.parameters = parameters
        self.Ndim = len(parameters)
        self.samples = numpy.empty((0,self.Ndim))
        self.truth = numpy.empty(self.Ndim)
        self.parstring=", ".join(self.parameters)
        
        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Probability density function'

# ----------------------------------------------------------------------------
# Add one sample to the ensemble:

    def append(self,sample):
        assert len(sample) == self.Ndim
        self.samples = numpy.append(self.samples,[sample],axis=0)
        return 

# ----------------------------------------------------------------------------
# accesscommands
    def call(self,key):
        assert key in self.parameters, "not a valid key. Keys are %s"%self.parstring
        for i in range(len(self.parameters)):
            if key == self.parameters[i]:
                return self.samples[:,i]



#=============================================================================

if __name__ == '__main__':
    p = PDF(['a','b'])
    x1 = numpy.array([1,2])
    x2 = numpy.array([3,3])
    p.append(x1)
    p.append(x2)
    print p.samples

#=============================================================================
