# ===========================================================================

import pangloss

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
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-21  Marshall & Collett (Oxford)
      2014-04-09  modified for BoRG, C Mason (UCSB)
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
# Extract samples in one parameter:
   
    def errmsg(self):
        return "not a valid parameter name. These are %s"%self.parstring

# ----------------------------------------------------------------------------
# Extract samples in one parameter:
   
    def getParameter(self,key):
        assert key in self.parameters, self.errmsg()
        for i in range(len(self.parameters)):
            if key == self.parameters[i]: return self.samples[:,i]

# ----------------------------------------------------------------------------
# Plot rough 1D histogram or 2D scatter plot:

    def plot(self,key1,key2=None,weight="weight",output=None,bins=None,title=None):

        import pylab as plt
        
        assert key1 in self.parameters, self.errmsg()
        if key2!=None:
            assert key2 in self.parameters, self.errmsg()
        if key2=="weight": 
            print "This code automatically uses the weights to form a histogram, you can state this explicitly using weightkey=weight if you have several weight columns"
            key2==None
        
        reformatnames={}
        
        reformatnames["mu"]="$\mu$"

        reformatnames["mu_cone"]="$\mu_{lc}$"
        reformatnames["kappa_cone"]="$\kappa_{lc}$"
        
        reformatnames["mu_tot"]="$\mu_{\mathrm{tot}}$"
        reformatnames["kappa_tot"]="$\kappa_{\mathrm{tot}}$"
        
        reformatnames["kappa_ext"]="$\kappa_{\mathrm{ext}}$"
        reformatnames["Kappah_median"]="$\widetilde{\kappa}_{\mathrm{halos}}$"

        key1name=key1[:]
        if key1name in reformatnames.keys():
            key1name = reformatnames[key1name]

        if key2 != None:
            key2name=key2[:]
            if key2name in reformatnames.keys():
                key2name = reformatnames[key2name]

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # 1D histogram:

        if key2 == None:

            par1 = self.getParameter(key1)
            
            # Check to see if samples are weighted:
            if weight!=None and weight in self.parameters: 
                weights = self.getParameter(weight)
                # print "pdf: plotting weighted histogram..."
            else: 
                weights = numpy.ones(len(par1))*1.0
                # print "pdf: plotting unweighted histogram..."
            weights /= numpy.sum(weights)

            if bins==None:
                nb=len(par1)*0.05
                if nb<20: nb=20
                bins=numpy.linspace(par1.min(),par1.max(),20)
            
            plt.figure()
            plt.hist(par1,weights=weights,bins=bins)
            plt.xlabel(key1name)
            plt.ylabel("P(%s|$\mathcal{D}$)"%key1name)
            if title != None: plt.title(title)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # 2D scatter plot:

        elif key2 != None:
            
            assert weight==None, "I don't know about weighted 2D plots yet."
            
            # print "pdf: plotting a 2D scatter plot ..."

            plt.figure()
            plt.scatter(self.getParameter(key1),self.getParameter(key2),\
                            c='k',s=2, edgecolors=None)
            plt.xlabel(key1name)
            plt.ylabel(key2name)
            if title != None: plt.title(title)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if output == None:
            plt.show(block=True)
        else:
            pangloss.rm(output)
            plt.savefig(output,dpi=300)
                
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
