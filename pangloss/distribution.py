from numpy.random import permutation, randint
from numpy import array

class PermuteTrueHaloMassAssignment(object):
    """
    Halo mass distribution that permutes the true halo mass values of the galaxies.
    """
    def draw(self, galaxies):
        return permutation(galaxies['Mhalo_obs'])

class TrueHaloMassDistribution(object):
    """
    Halo mass distribution that draws from the true halo mass distribution of the galaxies.
    """
    def draw(self, galaxies):
        n = len(galaxies)
        return array(galaxies['Mhalo_obs'])[randint(0, n, n)]