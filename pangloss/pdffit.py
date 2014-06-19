import numpy as np
import pylab as plt
import numpy.polynomial.hermite as H
from math import pi

class GaussHermiteModel:
    def __init__(self, n, mean, sd, uncertainty=1.):
        # Initial mean and sd are directly from the data
        self.initial_mean = mean
        self.initial_sd = sd
        self.initial_guess = np.array([mean, sd, 1. ,0., 0., 0., 0.])
        self.n = n
        self.uncertainty = uncertainty
    
    def B(self, x, order):
        """"
        Gauss-Hermite basis functions. Need to divide by function by sigma
        ----------------------------------------------------------------------- 
        Input: x = mu - <mu> / sigma
               order < 5
        """
        herm_poly = [1., 2.*x, 4*(x**2.) - 2., 8*(x**3.) - 12*x, 16*(x**4.) - 48*(x**2) + 12.]
        prefactor = (2.**float(order) * np.sqrt(pi) * np.factorial(order))**0.5
        return prefactor*herm_poly[order]*np.exp(-(x**2.)/2.)
    
    def fit_matrix(self, data, parameters):
        """
        A is matrix of basis functions (known - calculated from mean, sd)     
        """
        mean, sd = parameters
        self.x = (data - mean)/sd
        
        self.A = []
        for i in range(5): self.A.append(self.B(self.x, i))
        
    def predicted_values(self, Cn):
        """
        Fit for yi = Aij*Cj + ei, returns fitted values yfi
        -----------------------------------------------------------------------
        A is matrix of basis functions (known - calculated from mean, sd)
        Cn is vector of coefficients (unknown - we want to find this)        
        """
        self.Cn = Cn
        return np.dot(self.A, Cn)    
    
    def uncertainties(self, parameters):
        return self.uncertainty*np.ones(self.n)
    
    def generate(self, parameters):
        return (self.predicted_values(parameters) +
            np.random.normal(scale=self.uncertainty,size=self.n))   
            
            
            

# ===============================================  

def fit_linear_least_squares(model, data):
    mean = np.mean(data)
    sd = np.std(data)
    
    n_params = len(model.initial_guess)
    n_data = len(data)
    assert len(model.predicted_values(model.initial_guess))==n_data    
    
    coefficient_matrix = np.zeros((n_data,n_params))
    for i in range(n_params):
        params = np.zeros(n_params)
        params[i] = 1
        coefficient_matrix[:,i] = model.predicted_values(params)
        
    U,S,V = np.linalg.svd(arr, full_matrices=False)

    # v should be sorted.
    # this solution should be equivalent to v[1,0] / -v[1,1]
    # but I'm using this: http://stackoverflow.com/questions/5879986/pseudo-inverse-of-sparse-matrix-in-python
    M = V[-1,0]/-V[-1,-1]
    
    x, residues, rank, s = np.linalg.lstsq(coefficient_matrix, data)
    
    return x

# ===============================================   
     
def demo_linear_least_squares():
    true_parameters = np.array([2.,-1.])
    n = 10
    uncertainty = 0.1
    
    model = GaussHermiteModel(n,uncertainty)
    
    xs = np.arange(n)/float(n)
    data = model.generate(true_parameters)
    
    fit_parameters = fit_linear_least_squares(model, data)
    
    
    plt.figure()
    plt.errorbar(xs, data, uncertainty, label="data", fmt="+")
    plt.plot(xs,model.fitting_function(true_parameters), label="true value")
    plt.plot(xs,model.fitting_function(fit_parameters), label="fitted value")
    plt.xlim(0,1)
    plt.title("Simple linear least-squares fitting")
    plt.legend(loc="best")
    plt.savefig("linear-least-squares-1.png")
        
        
if __name__=='__main__':
  np.random.seed(0)
  demo_linear_least_squares()
  plt.show()            