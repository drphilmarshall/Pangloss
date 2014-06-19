import numpy as np
import pylab as plt

class SinCosModel:

    def __init__(self, n, uncertainty=1.):
        self.initial_guess = np.array([1.,0.])
        self.n = n
        self.uncertainty = uncertainty
    
    def fitting_function(self, parameters):
        a, b = parameters
        xs = np.arange(self.n)/float(self.n)
        return a*np.sin(2*np.pi*xs) + b*np.cos(2*np.pi*xs)
    
    def uncertainties(self, parameters):
        return self.uncertainty*np.ones(self.n)
    
    def generate(self, parameters):
        return (self.predicted_values(parameters) +
            np.random.normal(scale=self.uncertainty,size=self.n))
          
 

# ===============================================  

def fit_linear_least_squares(model, data):
    n_params = len(model.initial_guess)
    n_data = len(data)
    assert len(model.predicted_values(model.initial_guess))==n_data    
    
    coefficient_matrix = np.zeros((n_data,n_params))
    for i in range(n_params):
        params = np.zeros(n_params)
        params[i] = 1
        coefficient_matrix[:,i] = model.predicted_values(params)
        
    x, residues, rank, s = np.linalg.lstsq(coefficient_matrix, data)
    
    return x

# ===============================================   
     
def demo_linear_least_squares():
    true_parameters = np.array([2.,-1.])
    n = 10
    uncertainty = 0.1
    
    model = SinCosModel(n,uncertainty)
    
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