import pylab, atpy
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
import lightcone

D = distances.Distance()
D.h = 0.7
arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"

# First 99 lines of the same, for testing:
#datafile = "test.txt"
    
cat = atpy.Table(datafile, type='ascii')
print "Read in master table, length",len(cat)

def MCrelation(M200,MCerror=False):
    c_200 = 4.67*(M200/(10**14))**0.11 #Neto et al. equation 5
    return c_200

def Hsquared(z):
    H0 =D.h*3.241*10**-18
    Hsq=(H0**2.)*(D.OMEGA_M*(1.+z)**3.+(1.-D.OMEGA_M)) #Flat LambdaCDM only at this stage
    return Hsq

def rho_crit_univ(z):   #critical density of the universe at z
    ro= (2.642*10**46)*Hsquared(z) #units of solar mass per cubic megaparsec, H(z) must be in units of per second.
    return ro 

def delta_c(c):
    return (200./3)*(c**3)/(numpy.log(1+c)-c/(1+c))


M200 = cat['M_Subhalo[M_sol/h]']
c200 = MCrelation(M200)
cat.add_column('NetoC', c200)
rhocrit=rho_crit_univ(cat['z_spec'] )
r200 = (3*M200/(800*3.14159*rhocrit))**(1./3)
cat.add_column('r200TRUE', r200)

rs = r200/c200
rhos = delta_c(c200)*rhocrit 

truncationscale=[1,2,3,5,8,10]
for i in range(len(truncationscale)):
    R_trunc=truncationscale[i]*r200
    
    mass=4*3.14159*rhos*(rs**3)  *     \
        (  numpy.log(1+(R_trunc)/rs) - \
               R_trunc/(rs+R_trunc)    \
        )
    cat.add_column('Mtrunc%iR200'%truncationscale[i], mass)

filename="testtable.table"
print "Table printed to", filename
cat.write(filename, type='ascii')

nplanes=200
zp,dz=numpy.linspace(0,self.zs,nplanes,endpoint=False,retstep=True)

