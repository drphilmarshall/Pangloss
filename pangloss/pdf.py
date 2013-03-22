# ===========================================================================

import numpy

# ============================================================================

class PDF(object):
    """
    TEST DOCSTRING
    """

# ----------------------------------------------------------------------------

    def __init__(self,parameters):
        
        self.name = 'Probability Density Function'
        if type(parameters) != list: parameters = [parameters]
        self.parameters = parameters
        self.Ndim = len(parameters)
        self.samples = numpy.empty((0,self.Ndim))
        self.truth = numpy.empty(self.Ndim)
        
        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Probability density function'

# ----------------------------------------------------------------------------
# Add one sample to the ensemble:

    def append(self,sample):
    
        assert len(sample) == self.Ndim
        assert sample.shape == (self.Ndim,)
        self.samples = numpy.append(self.samples,[sample],axis=0)
    
        return 

#=============================================================================

if __name__ == '__main__':
    p = PDF(['a','b'])
    x1 = numpy.array([1,2])
    x2 = numpy.array([3,3])
    p.append(x1)
    p.append(x2)
    print p.samples

#=============================================================================
