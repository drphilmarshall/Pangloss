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
    def getParameter(self,key):
        assert key in self.parameters, "not a valid key. Keys are %s"%self.parstring
        for i in range(len(self.parameters)):
            if key == self.parameters[i]:
                return self.samples[:,i]

# ----------------------------------------------------------------------------
# plotting commands
    def plot(self,key1,key2=None,weightkey="weight",bins=None,title=None):
        import pylab as plt
        assert key1 in self.parameters, "not a valid key. Keys are %s"%self.parstring
        if key2!=None:assert key2 in self.parameters, "not a valid key. Keys are %s"%self.parstring

        if key2=="weight": 
            print "This code automatically uses the weights to form a histogram, you can state this explicitly using weightkey=weight if you have several weight columns"
            key2==None
        
        reformatnames={}
        reformatnames["kappa_ext"]="$\kappa_{\mathrm{ext}}$"
        reformatnames["kappa_h_median"]="$\widetilde{\kappa}_{\mathrm{halos}}$"
        key1name=key1[:]
        if key1name in reformatnames.keys():
            key1name = reformatnames[key1name]

        if key2 != None:
            key2name=key2[:]
            if key2name in reformatnames.keys():
                key2name = reformatnames[key2name]

        if key2 == None:
            par1=self.getParameter(key1)
            #make a histogram
            
            #check to see if samples are weighted:
            print "pdf: looking for a weight key"
            if weightkey!=None and weightkey in self.parameters: 
                weights=self.getParameter(weightkey)
                print "pdf: plotting weighted histogram..."
            else: 
                weights= numpy.ones(len(par1))*1.0
                print "pdf: plotting unweighted histogram..."
            weights/=numpy.sum(weights)
            
            if bins==None:
                bins=numpy.linspace(par1.min(),par1.max(),len(par1)*0.05)
        
            plt.figure()
            plt.hist(par1,weights=weights,bins=bins)
            plt.xlabel(key1name)
            plt.ylabel("P(%s|$\mathcal{D}$)"%key1name)
            if title != None: plt.title(title)


        if key2 != None:
            assert weightkey==None, "I don't know about weighted 2D plots yet."
            print "pdf: plotting a 2D scatter plot ..."
            #make a 2d scatterplot (upgrade to a cornerplot in future!)
            plt.figure()
            plt.scatter(self.getParameter(key1),self.getParameter(key2),\
                            c='k',s=1,edgecolors=None)
            plt.xlabel(key1name)
            plt.ylabel(key2name)
            if title != None: plt.title(title)

        plt.show(block=True)
        return None

#=============================================================================

if __name__ == '__main__':
    p = PDF(['a','b'])
    x1 = numpy.array([1,2])
    x2 = numpy.array([3,3])
    p.append(x1)
    p.append(x2)
    print p.samples

#=============================================================================
