import numpy as np
import pylab as plt
import scipy
import numpy.polynomial.hermite as H
from math import pi

class GaussHermiteModel:
    def __init__(self, mean, sd, uncertainty=0.0):
        # Initial mean and sd are directly from the data
        self.initial_mean = mean
        self.initial_sd = sd
        self.initial_parameters = np.array([1. ,0., 0., 0., 0.])
        
    def B(self, x, order):
        """"
        Gauss-Hermite basis functions. Need to divide by sqrt(sigma) in the matrix
        ----------------------------------------------------------------------- 
        Input: x = mu - <mu> / sigma
               order < 5
        """
        herm_poly = [1., 2.*x, 4*(x**2.) - 2., 8*(x**3.) - 12*x, 16*(x**4.) - 48*(x**2) + 12.]
        prefactor = (2.**float(order) * np.sqrt(pi) * np.math.factorial(order))**-0.5
        return prefactor*herm_poly[order]*np.exp(-(x**2.)/2.)
    
    def matrix_row(self, datapoint, mean, sd):
        """
        Fill each row with basis functions (known - calculated from mean, sd)     
        """
        x = (datapoint - mean)/sd        
        row = np.array([self.B(x, 0), self.B(x, 1), self.B(x, 2), self.B(x, 3), self.B(x, 4)])
        self.row = row/sd
#        self.row = row/np.sqrt(sd)
        
        return self.row
    
    def fill_matrix(self, data, mean, sd):
        """
        A is matrix of basis functions (known - calculated from mean, sd)     
        Fill each row of the matrix with the basis using for a given data point
        """
        n_points = len(data)
        n_params = len(self.initial_parameters)
        
        self.A = np.zeros((n_points,n_params))
        for i in range(n_points):
            self.A[i,:] = self.matrix_row(data[i], mean, sd) 

        return self.A
        
    def predicted_values(self, A, Cn):
        """
        Fit for yi = Aij*Cj + ei, returns fitted values yfi
        -----------------------------------------------------------------------
        A is matrix of basis functions (known - calculated from mean, sd)
        Cn is vector of coefficients (unknown - we want to find this)        
        """
        self.Cn = Cn
        return np.dot(A, Cn)     
                             

# ===============================================  

def fit_linear_least_squares(model, measured, points):
    mean = np.mean(measured)
    sd = np.std(measured) 
    
    # Fill each row of the matrix with the basis using for a given data point
    # (data is the mu values)
    coefficient_matrix = model.fill_matrix(points, mean, sd)

    # Calculate the Least Squares solution from the matrix we just made
    C, residues, rank, s = np.linalg.lstsq(coefficient_matrix, measured)
#    C, residues = scipy.optimize.nnls(coefficient_matrix, measured)
    print 'singular values:',s
    return C, coefficient_matrix, residues

# ===============================================   
     
def demo_linear_least_squares():
    """
    Demo with a true gaussian
    """
    pdf = np.genfromtxt("../../BORG/LensingModifications/pangloss/figs/borg/all_PofMu.txt")            
    points = pdf[:,0]
    measure = pdf[:,1]
    mean, sd = np.mean(pdf), np.std(pdf)
    
    # Set up the true Cn
    # Make a Gaussion over the range 0-3 with mean=1, sd=0.2
    '''
    true_parameters = np.array([1., 0.0, 0., 0., 0.])
    points = np.linspace(-3.,3.,1000)
    mean, sd = 1.0, 0.5
    '''
    print mean, sd
    
    # Initialise the model with given mean and sd
    model = GaussHermiteModel(mean, sd) 
    
    # Build the measured values      
#    measure_pure = np.exp(-(points - mean)**2./(2.*(sd**2.)))/(np.sqrt(2.*pi)*sd)
#    measure_matrix = model.fill_matrix(points, mean, sd)
    #MTM = np.dot(measure_matrix.T, measure_matrix)
    #detMTM = np.linalg.det(MTM)
    #print detMTM
#    measure = model.predicted_values(measure_matrix, true_parameters)
                                         
    # Find the Cn and matrix by fitting the model to the 'measured' data and the data points    
    fit_parameters, matrix, residues = fit_linear_least_squares(model, measure, points)
    print 'residues:',residues  
    print 'parameters:',fit_parameters  
    
    plt.figure()
    plt.plot(points, measure, label="true values")
#    plt.plot(points, model.predicted_values(matrix, true_parameters), label="true values via matrix")
    plt.plot(points, model.predicted_values(matrix, fit_parameters), label="fitted values")
    plt.title("Simple linear least-squares fitting")
    plt.legend(loc="best")
    plt.savefig("linear-least-squares-hermitegauss1.png")
        
        
if __name__=='__main__':
  np.random.seed(0)
  demo_linear_least_squares()
  plt.show()            