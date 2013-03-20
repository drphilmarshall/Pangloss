import numpy
import distances
from scipy import interpolate,optimize
import cPickle

import Relations as Rel
#import pylab as plt

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

vb = False

# ============================================================================
# Attempt at good Object Orientation. Made code slower :-(.   
class Lens_Plane():
    def __init__(self, z,zl,zs,Da_l,Da_s,Da_ls):
        if z==0: z+=1e-99
        D=distances.Distance()
        self.name = "A sheet at redshift %.3f"%z
        self.z=z
        self.Da_p=D.Da(0,z)
        self.rho_crit=D.rho_crit_univ(z)
        self.Da_l,self.Da_s,self.Da_ls=Da_l,Da_s,Da_ls
        self.Da_ps = D.Da(z,zs)
        self.Da_pl = D.Da(z,zl)
        self.sigma_crit=(1.663*10**18)*(self.Da_s/(self.Da_p*self.Da_ps))
        if self.z>zl: #1 is lens, 2 is perturber
            D1s = self.Da_ls
            D2  = self.Da_p
        else: #1 is perturber, 2 is lens
            D1s = self.Da_ps
            D2  = self.Da_l
        D12= self.Da_pl
        self.beta=(D12*self.Da_s)/(D2*D1s)

# ============================================================================
class Grid():
    def __init__(self,zl,zs,nplanes=100,cosmo=[0.25,0.75,0.73]): 
        if zl > zs: zl,zs=zs,zl
        D=distances.Distance()
        self.name = '1D Redshift grid of lens planes, each containing precalculated quantities '
        self.zmax=zs*1.0
        self.zs=zs*1.0
        self.zltrue=zl
        self.nplanes=nplanes
        self.cosmo=cosmo
        #These are the plane redshifts:
        self.redshifts,self.dz=self.redshiftbins=numpy.linspace(0.0,self.zmax,self.nplanes,endpoint=True,retstep=True)
        self.redshifts+=(self.dz/2.)
        # SNAP LENS TO GRID:
        self.zl=self.snap([zl])[0][0]
        self.Da_l=D.Da(zl)
        self.Da_s=D.Da(zs)
        self.Da_ls=D.Da(zl,zs)
        self.plane={}

        #nice object orientation: can't call easily
        #for p in range(nplanes):
        #    z=self.redshifts[p]
        #    self.plane[p] = Lens_Plane(z,zl,zs,self.Da_l,self.Da_s,self.Da_ls)

        #fast and nasty python :-(
        self.Da_p=numpy.zeros(len(self.redshifts))
        self.rho_crit=numpy.zeros(len(self.redshifts))
        self.Da_ps=numpy.zeros(len(self.redshifts))
        self.Da_pl=numpy.zeros(len(self.redshifts))
        self.sigma_crit=numpy.zeros(len(self.redshifts))
        self.beta=numpy.zeros(len(self.redshifts))
        for i in range(len(self.redshifts)):
            z=self.redshifts[i]
            self.Da_p[i]=D.Da(0,z)
            self.rho_crit[i]=D.rho_crit_univ(z)
            self.Da_ps[i] = D.Da(z,zs)
            self.Da_pl[i] = D.Da(z,zl)
            self.sigma_crit[i]=(1.663*10**18)*(self.Da_s/(self.Da_p[i]*self.Da_ps[i]))
            if z>zl: #1 is lens, 2 is perturber
                D1s = self.Da_ls
                D2  = self.Da_p[i]
            else: #1 is perturber, 2 is lens
                D1s = self.Da_ps[i]
                D2  = self.Da_l
            D12= self.Da_pl[i]
            self.beta[i]=(D12*self.Da_s)/(D2*D1s)

    def snap(self,z):
        snapped_p=numpy.digitize(z,self.redshifts-(self.dz)/2)-1
        snapped_p[snapped_p<0]=0 #catalogue has some blue shifted objects!
        snapped_z=self.redshifts[snapped_p]
        return snapped_z,snapped_p

    def __str__(self):
        return '1-D Grid of %i planes seperated in redshift by dz= %f' % (self.nplanes,self.dz)

# ---------------------------------------------------------------------------
if __name__ == '__main__':
    nplanes=100
    g=Grid(0.6,1.4,nplanes=nplanes)
    TEST=open("/data/tcollett/Pangloss/grid%i.grid"%nplanes,"wb")
    cPickle.dump(g,TEST,2)
#TEST.close()

"""
FILE=open("test.pkl","r")
k=cPickle.load(FILE)
print k.plane[3].beta
nn=1000
reg=numpy.zeros(nn)
for i in range(nn):
    reg[i]= numpy.log10(k.plane[3].drawMhalo(k.plane[3].drawMstar([1e11])))
plt.hist(reg)
plt.show()
print numpy.log10(k.plane[3].drawMstar([10**10.1]))
"""
