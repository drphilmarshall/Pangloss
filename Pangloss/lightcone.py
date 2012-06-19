'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

description:
------------

Given a catalog, drill out a narrow pencil beam defined by position
vector  xc and some radius.

# NB. Anything marked '###' is still unfinished/preliminary and 
probably wrong! #### are finished functions that haven't been 
bug-checked.

to-do:
------
- finish basic functions
- carry out curve of growth test
- decide whether to continue

issues:
-------
- stellar mass is not included yet
- so lines of sight passing close to halos will not be accurate
- Born approximation will break down in this case anyway
- Keeton approximation may well be way off, need better approx from 
    Roger, Sherry and Stefan.
- lightcone and lenslightcone should be separate classes, where llc 
    inherits from lc



In need of array optimization
-----------------------------
- The LensingProfiles class
- The beta function.
- Mstar to MBehrooziHalo



The warning: divide by zero comes from the grid class, but it isn't causing any errors.
'''

# ======================================================================

import pylab
import matplotlib.pyplot as plt
import numpy, numpy.random as rnd, atpy
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
import LensingProfiles as LP
import LensingFunc as LF
import grid as GRID
import Relations as Rel
import cPickle
#import time
#t0=time.clock()    

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

vb = False

# ============================================================================
 
class lightcone(object):
    def __init__(self,catalog,position,radius,magnitudecut=99,band="r"):
        self.name = 'Lightcone through the observed Universe'
        self.catalog = catalog
        self.xmax = self.catalog['pos_0[rad]'].max()
        self.xmin = self.catalog['pos_0[rad]'].min()
        self.ymax = self.catalog['pos_1[rad]'].max()
        self.ymin = self.catalog['pos_1[rad]'].min() 
        self.rmax = radius
        self.xc = [position[0],position[1]]

        dx = self.rmax*arcmin2rad
        self.galaxies = self.catalog.where((self.catalog['pos_0[rad]'] > (self.xc[0]-dx)) & \
                                           (self.catalog['pos_0[rad]'] < (self.xc[0]+dx)) & \
                                           (self.catalog['pos_1[rad]'] > (self.xc[1]-dx)) & \
                                           (self.catalog['pos_1[rad]'] < (self.xc[1]+dx))   )

 
        x = (self.galaxies['pos_0[rad]'] - self.xc[0])*rad2arcmin
        y = (self.galaxies['pos_1[rad]'] - self.xc[1])*rad2arcmin
        r = numpy.sqrt(x*x + y*y)
        self.galaxies.add_column('x',x)
        self.galaxies.add_column('y',y)
        self.galaxies.add_column('r',r)
        self.galaxies = self.galaxies.where(self.galaxies.r < self.rmax)
        F814 = (self.galaxies.mag_SDSS_i + self.galaxies.mag_SDSS_z)/2. #approximate F814 colour.
        self.galaxies.add_column("mag_F814W",F814)
        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W"
        else:
            col = "mag_%s" % band
        self.galaxies=self.galaxies.where(self.galaxies["%s"%col] < magnitudecut)

        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Lightcone of radius %.2f arcmin, centred on (%.3f,%.3f) rad' % (self.rmax,self.xc[0],self.xc[1])

# ----------------------------------------------------------------------------
   #Tell me the number of galaxies within a certain radius, that pass a certain magnitude cut.
    def N_radius_cat(self,radius,cut=[18.5,24.5], band="F814W", radius_unit="arcsec"):
        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W"
        else:
            col = "mag_%s" % band
        if radius < 10: print "Warning: Default units for N_radius are arcsec!"
        if radius_unit == "arcsec":
            radius = radius/60.
        if col != "warning":
            self.N_cut=self.galaxies.where((self.galaxies.r < radius)  & \
                                              (self.galaxies["%s"%col] < cut[1])& \
                                          (self.galaxies["%s"%col] > cut[0]))

            return self.galaxies.where((self.galaxies.r < radius) & \
                                          (self.galaxies["%s"%col] < cut[1])& \
                                          (self.galaxies["%s"%col] > cut[0]))

    def N_radius(self,radius,cut=[18.5,24.5],band="F814W", radius_unit="arcsec"):
        Ntable=self.N_radius_cat(radius,cut,band, radius_unit)
       #print len(Ntable.r)
       #plt.hist(Ntable.r)
       #plt.show()
        return len(Ntable.r)

# ============================================================================
# ============================================================================

class lens_lightcone(lightcone):
    def __init__(self,catalog,position,radius,zl,zs,nplanes=50,cosmo=[0.25,0.75,0.73],magnitudecut=99,band="r",grid=None):
        self.name = 'Snapped Lens Lightcone through the observed Universe'
        lightcone.__init__(self,catalog,radius,position,magnitudecut=magnitudecut,band=band)
        if grid==None:
            self.grid=GRID.lensgrid(zl,zs,nplanes=nplanes,cosmo=cosmo)
            self.grid.populatelensgrid()
        else: self.grid=grid


        #correct M_subhalos that == 0
        #Msub=self.galaxies['M_Subhalo[M_sol/h]']
        #Mhal=self.galaxies['M_Halo[M_sol/h]']
        #Msub[Msub==0]=Mhal[Msub==0]
        #self.galaxies.remove_columns("M_Subhalo[M_sol/h]")
        #self.galaxies.add_column("M_Subhalo[M_sol/h]",Msub)


        self.galaxies=self.galaxies.where(self.galaxies.z_spec < zs)

        zsnapped,psnapped=self.grid.snap(self.galaxies.z_spec)

        self.galaxies.add_column('zsnapped',zsnapped)
        self.galaxies.add_column('psnapped',psnapped)


        zd = self.galaxies.zsnapped
        p  = self.galaxies.psnapped

        self.zl = self.grid.snap([zl])
        self.zs = self.grid.snap([zs])

        # Grab relevant quantities from the grid
        self.galaxies.add_column('Da_d', self.grid.Da_p[p])
        self.galaxies.add_column('Da_ds', self.grid.Da_ps[p])
        self.galaxies.add_column('Da_dl', self.grid.Da_pl[p])
        self.galaxies.add_column('beta', self.grid.beta_p[p])
        self.galaxies.add_column('rho_crit',  self.grid.rho_crit_p[p])    # Mpc
        self.galaxies.add_column('sigma_crit', self.grid.sigma_crit_p[p]) #solarmasses per Mpc^2

        #calculate rphys
        self.galaxies.add_column('rphys',self.galaxies.Da_d*self.galaxies.r*arcmin2rad) #units 


        return None

# ----------------------------------------------------------------------------

    def make_kappa_contributions(self,BehrooziHalos=False,hardcut="Rvir",truncationscale=5,scaling="add",errors=True,BehrooziSpline=None,eBer=1e-99,centralsonly=False, Mh2Mh=False, Mstar2Mh=False, perfectsatellites=False): 

        if centralsonly==True:
            self.galaxies=self.galaxies.where(self.galaxies.Type==0)
            #self.galaxies=self.galaxies.where(self.galaxies.Type!=2)





        self.galaxies['M_Subhalo[M_sol/h]'][self.galaxies['M_Subhalo[M_sol/h]']==0.0]=1.0

        if Mh2Mh:
                if BehrooziSpline==None:
                    HALOSTARlowz,STARHALOlowz,HALOSTARhighz,STARHALOhighz=Rel.Behroozi_Spline()
                else:
                    HALOSTARlowz,STARHALOlowz,HALOSTARhighz,STARHALOhighz=BehrooziSpline[0],BehrooziSpline[1],BehrooziSpline[2],BehrooziSpline[3],
                M200 = Rel.Mhalo_to_Mhalo(self.galaxies['M_Subhalo[M_sol/h]'],self.galaxies['z_spec'],HALOSTARlowz,STARHALOlowz,HALOSTARhighz,STARHALOhighz,eBer=eBer)

        elif Mstar2Mh:
                Mstar=10**(numpy.log10(self.galaxies['M_Stellar[M_sol/h]'])+rnd.normal(0,eBer,len(self.galaxies.z_spec)))
                M200=Rel.Mstar_to_M200(Mstar,self.galaxies['z_spec'],scatter=False,)

        else: M200=self.galaxies['M_Subhalo[M_sol/h]']


        if perfectsatellites==True:
            M200[self.galaxies.Type==1]=self.galaxies['M_Subhalo[M_sol/h]'][self.galaxies.Type==1]
            M200[self.galaxies.Type==2]=self.galaxies['M_Subhalo[M_sol/h]'][self.galaxies.Type==2]




        c200 = Rel.MCrelation(M200,MCerror=errors)
        r200 = (3*M200/(800*3.14159*self.galaxies.rho_crit))**(1./3)

        #print (M200/self.galaxies['M_Subhalo[M_sol/h]']).max()


        r_s = r200/c200
        rho_s = LP.delta_c(c200)*self.galaxies.rho_crit
        kappa_s = rho_s * r_s /self.galaxies.sigma_crit
    
        x=self.galaxies.rphys/r_s
        
        if hardcut == "Rvir" or hardcut=="RVir" or hardcut == "r_vir" or hardcut == "rvir":
            R_trunc=truncationscale*r200
        elif hardcut == "rs" or hardcut=="Rs" or hardcut == "R_s" or hardcut == "r_s":
            R_trunc=truncationscale*r_s
        elif hardcut== "subhalo":
            R_trunc=3*r200
            R_trunc[self.galaxies.Type == 1]=(r200[self.galaxies.Type == 1])*truncationscale
        else: print "what hardcut did you mean?"


        

        mass=4*3.14159*rho_s*(r_s**3)  *     \
            (  numpy.log(1+(R_trunc)/r_s) - \
                   R_trunc/(r_s+R_trunc)    \
            )

        #print "boo"
        self.galaxies.add_column('Mtrunc', mass)
        #print "hoo"

        kappaNFW=kappa_s*1.0
        shearNFW=kappa_s*1.0
        for i in range(len(x)):
               #treat as NFW if within truncation radius:
               if self.galaxies.rphys[i]<R_trunc[i]: 
                   kappaNFW[i]*=LP.Ffunc([x[i]])
                   shearNFW[i]*=LP.Gfunc([x[i]])
               #treat as point mass if outside truncation radius:
               else:
                   kappaNFW[i]*=0.0
                   shearNFW[i]=(((mass[i])\
                                    /(3.14159*( self.galaxies.rphys[i])**2))\
                                    /self.galaxies.sigma_crit[i])
                                    
        
        #-------------------------------------------------------
        # Now computer starlight lensing component.

        #Mstar=self.galaxies['M_Stellar[M_sol/h]']
        #reff=Rel.reffFromMass(Mstar) # in kiloparsecs. # Remember self.galaxies.rphys[i] is in Mpc so
        #reff=reff/1000 # now in megaparsecs :-)
        #kappaStar,shearStar=Mstar*LP.sersic(self.galaxies.rphys,reff)/self.galaxies.sigma_crit
        #should ^ these be times Mstar or not? TC thinks so.
        ###
        ####
        # This is disabled currently - it gives silly answers for a minority of LoS that get very cose to the galaxies (SL lines where keeton approximation breaks.)

        #-------------------------------------------------------
        #compute a fudge component (nonNFW additional component to the halos), just to see what happens
        #kappaFudge=kappaNFW*0#.5

        #-------------------------------------------------------
        kappa = kappaNFW #+ kappaStar
        shear = shearNFW #+ shearStar

        self.galaxies.add_column('kappa',kappa)
        self.galaxies.add_column('gamma',shear)

        kappa_Scaled= LF.KappaScale_beta(self.galaxies.beta,self.galaxies.kappa,self.galaxies.gamma,scaling=scaling)
        self.galaxies.add_column('kappa_Scaled',kappa_Scaled)
        
        self.kappa_Scaled_total=numpy.sum(kappa_Scaled)
        
        kappa_fail=kappa_Scaled[kappa_Scaled<0]
        #if kappa_fail!=[]:print kappa_fail
        
        return None

# ----------------------------------------------------------------------------

   # 2-panel plot, showing view from Earth and also line of sight section:
    def plot(self,starlight=False,dmglow=False,kappa_indiv=False,kappa_keeton=False,observed_light=True):
       scale=numpy.max([ numpy.absolute((numpy.min(self.galaxies.kappa_Scaled))), numpy.max(self.galaxies.kappa_Scaled)])/200
       scale2= (numpy.max(self.galaxies.kappa))/200
       # Galaxy positions:
       ax1=plt.subplot(2,1,1, aspect ='equal')
       plt.subplot(2,1,1,aspect='equal')
       empty = True

       # plot different properties depending on options:
       if dmglow==True:
         #plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Subhalo[M_sol/h]']))/3),edgecolor = 'none' )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2),edgecolor = 'none' )     
         empty = False
       if kappa_keeton:
         for galaxy in self.galaxies:
            kappa = galaxy['kappa_Scaled']
            if kappa <0:
               plt.scatter(galaxy['x'], galaxy['y'], c='b', marker='o',s=-kappa/scale , edgecolor = 'none')    
            else:
               plt.scatter(galaxy['x'], galaxy['y'],  c='r', marker='o',s=kappa/scale , edgecolor = 'none')
         empty = False
       if kappa_indiv: 
         plt.scatter(self.galaxies.x, self.galaxies.y, c='none', marker='o',s=(self.galaxies.kappa)/scale2)     
         empty = False
       if observed_light==True:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=2**(25-(self.galaxies['mag_SDSS_r'])),edgecolor = 'none' )     
         empty = False


       if empty:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=1)      

       # Lightcone boundary and centroid:
       circ=pylab.Circle(self.xc,radius=self.rmax,fill=False,linestyle='dotted')
       ax1.add_patch(circ)
       ax1.plot([self.xc[0]],[self.xc[1]], c='k', marker='+',markersize=10)


       # Labels:
       plt.xlabel('x / arcmin')
       plt.ylabel('y / arcmin')

       #axis limits
       ax1.axis([self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1])

       #plt.title('%s' % self)


#      -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


       # Phil Marshall's subplot  View along redshift axis: 

       plt.subplot(2,1,2)
       empty = True
       

       if dmglow==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='k', marker='o',s=(((self.galaxies['M_Subhalo[M_sol/h]']))*3e-12),edgecolor = 'none' )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2),edgecolor = 'none' )     
         empty = False
       if kappa_keeton:
         for galaxy in self.galaxies:
            kappa = galaxy['kappa_Scaled']
            if kappa <0:
               plt.scatter(galaxy['z_spec'], galaxy['y'], c='b', marker='o',s=-kappa/scale , edgecolor = 'none')    
            else:
               plt.scatter(galaxy['z_spec'], galaxy['y'],  c='r', marker='o',s=kappa/scale , edgecolor = 'none')
         empty = False
       if kappa_indiv: 
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='none', marker='o',s=(self.galaxies.kappa)/scale2)     
         empty = False
       if observed_light==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='y', marker='o',s=2**(25-(self.galaxies['mag_SDSS_r'])),edgecolor = 'none' )     
         empty = False


       if empty:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='k', marker='o',s=1)


       plt.xlabel('Redshift')
       plt.ylabel('y / arcmin')
       zmax = self.galaxies['z_spec'].max()


       plt.axis([0,self.zs+0.1,0,self.rmax+0.1])

       plt.xlabel('redshift z')
#        plt.ylabel('distance from LoS / arcmin')
       plt.ylabel('y / arcmin')
       
       zmax = max(self.galaxies['z_spec'].max(),self.zs)
#        print self.zs,self.galaxies['z_spec'].max()
#        print zmax
       plt.axis([0,zmax+0.1,-self.rmax-0.1,self.rmax+0.1])

       # Add lines marking source and lens plane, and optical axis:
       plt.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dotted',label='bla')
       plt.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='dotted')
       plt.axhline(y=0.0, xmin=0.0, xmax=zmax, color='black', ls='dashed')

 
       return None


# # ----------------------------------------------------------------------------

    def curve_of_growth(self,ordering="contribution",starlight=False,dmglow=False,kappa_indiv=False,kappa_keeton=False,observed_light=True):
       plt.clf()
       scale=numpy.max([ numpy.absolute((numpy.min(self.galaxies.kappa_Scaled))), \
                            numpy.max(self.galaxies.kappa_Scaled)])/200
       scale2= (numpy.max(self.galaxies.kappa))/200

       zeropoint=-self.kappa_expected()
       zero=numpy.ones(1)*zeropoint


       plt.subplot(2,1,1)
       #set up things to plot:
       if ordering=="contribution":
          args=numpy.argsort(-numpy.absolute(self.galaxies.kappa_Scaled))
       if ordering=="distance":
          args=numpy.argsort((self.galaxies.r))
          dist=numpy.concatenate((zero,numpy.take(self.galaxies.r,args)))
       if ordering=="r_mag":
          args=numpy.argsort((self.galaxies['mag_SDSS_r']))
          rmag=numpy.concatenate((zero,numpy.take(self.galaxies['mag_SDSS_r'],args)))
       ordered_kappas=numpy.take(self.galaxies.kappa_Scaled,args)


       cat=numpy.concatenate((zero,(ordered_kappas)))
       cumtot=numpy.cumsum(cat)
       n=range(len(args)+1)
      
       #work out limits for plot
       #maximum=ordered_kappas.max()      
       #minimum=ordered_kappas.min()
       #if numpy.absolute(maximum) < numpy.absolute(minimum):
       #   max = minimum
       #   min = numpy.max([maximum,0])
       #else: 
       #   max = maximum
       #   min = numpy.min([minimum,0])

       #choose first m points
       #n1=n[:20]
       #cumtot1=cumtot[:20]
       #corresponding_r1=corresponding_r[:20]

       if ordering=="contribution":
          plt.plot(n,cumtot)
          #plt.xlim([0,200])
          plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
          plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')
       if ordering=="distance":
          plt.plot(dist,cumtot)
          #plt.xlim([0,200])
          plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
          plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')
          plt.xlabel('LOS distance (arcmin)')
       if ordering=="r_mag":
          plt.plot(rmag,cumtot)
          #plt.xlim([0,200])
          plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
          plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')
          plt.xlim([16,26])
       plt.ylabel('$\kappa_{ext}$ (cumulative)')


#      -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


#      #Redshift axis vs LOS distance


       plt.subplot(2,1,2)
       empty = True

       self.galaxies=self.galaxies.where(self.galaxies.r<8)
       if kappa_keeton:
         for galaxy in self.galaxies:
            kappa = galaxy['kappa_keeton']
            if kappa <0:
               plt.scatter(galaxy['z_spec'], galaxy['r'], c='b', marker='o',s=-kappa/scale , edgecolor = 'none')    
            else:
               plt.scatter(galaxy['z_spec'], galaxy['r'],  c='r', marker='o',s=kappa/scale , edgecolor = 'none')
         empty = False
       if kappa_indiv: 
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='none', marker='o',s=(self.galaxies.kappa)/scale2)     
         empty = False

       if dmglow==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='k', marker='o',s=((numpy.log(self.galaxies['M_Subhalo[M_sol/h]']))/3),edgecolor = 'none' )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='g', marker='o',s=2**(numpy.log10(self.galaxies['M_Stellar[M_sol/h]'])-7),edgecolor = 'none' )     
         empty = False

       if observed_light==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='y', marker='o',s=2**(25-(self.galaxies['mag_SDSS_r'])),edgecolor = 'none' )     
         empty = False


       if empty:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='k', marker='o',s=1)


       plt.xlabel('Redshift')
       plt.ylabel('LoS distance / arcmin')
       zmax = self.galaxies['z_spec'].max()


       plt.axis([0,self.zs+0.1,0,self.rmax+0.1])

       #add lines for source and lens plane
       plt.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dashed',label='bla')
       plt.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='solid')       # plt.title('%s' % self)

# ============================================================================

# TESTS:

def test1(catalog):

    plt.clf()
    rmax = 3    
    zs=1.4 

    xpos = -0.008
    ypos = -0.008
    zl =  1.0 
    xc = []
    #xc = [xpos,ypos,zl] #leave as [] to select a lens at random
    rmax = 5

    xmax = catalog['pos_0[rad]'].max()
    xmin = catalog['pos_0[rad]'].min()
    ymax = catalog['pos_1[rad]'].max()
    ymin = catalog['pos_1[rad]'].min()

    x = rnd.uniform(xmin+rmax*arcmin2rad,xmax-rmax*arcmin2rad)
    y = rnd.uniform(ymin+rmax*arcmin2rad,ymax-rmax*arcmin2rad)

    xc = [x,y]

    lc = lens_lightcone(catalog,rmax,xc,zl,zs)

    lc.make_kappa_contributions(hardcut="RVir",truncationscale=3,scaling="add",errors=False)


    print "Total external convergence =",numpy.sum(lc.galaxies.kappa_Scaled)
    # Now make an illustrative plot:
    
    if rmax < 8:
       print "Plotting objects in lightcone..."
       lc.plot(starlight=False,dmglow=False, kappa_indiv=False, kappa_keeton=True,observed_light=True)
       

       pngfile = 'test.png'
       plt.savefig(pngfile)
       print "Plot saved in",pngfile

    print "Plotting curve of growth..."
    plt.clf()
    lc.curve_of_growth(ordering="distance")
    pngfile = 'curve_of_growth.png'
    plt.savefig(pngfile)
    print "Plot saved in",pngfile
    

    return

#-------------------------------------------------------------------------

def test2(catalog): 

    rmax = 5

    xmax = catalog['pos_0[rad]'].max()
    xmin = catalog['pos_0[rad]'].min()
    ymax = catalog['pos_1[rad]'].max()
    ymin = catalog['pos_1[rad]'].min()

    x = rnd.uniform(xmin+rmax*arcmin2rad,xmax-rmax*arcmin2rad)
    y = rnd.uniform(ymin+rmax*arcmin2rad,ymax-rmax*arcmin2rad)


 
    xc = [x,y]
    #print xc
    xc=[-0.007537950324046427, -0.006001515672089354]

    zl=0.6
    zs=1.4
    lc = lens_lightcone(catalog,rmax,xc,zl,zs)
    lc.make_kappa_contributions(errors=False)

    print numpy.max(lc.galaxies.kappa_Scaled[lc.galaxies.Type==2])

    return

#-------------------------------------------------------------------------

def test3(catalog):
    rmax = 5
    truncationscale=10
    ncones=12
    zl=0.6
    zs=1.4

    xmax = catalog['pos_0[rad]'].max()
    xmin = catalog['pos_0[rad]'].min()
    ymax = catalog['pos_1[rad]'].max()
    ymin = catalog['pos_1[rad]'].min()

    H0=numpy.linspace(0.65,0.83,30,endpoint=True)
    Om=numpy.linspace(0.2,0.4,21,endpoint=True)

    print Om

    import smooth
    datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"
    data=[datafile]
    kappa_empty_true= smooth.smooth(zl,zs,data,truncationscale=truncationscale,hardcut="Rvir",nplanes=50,cosmo=[0.25,0.75,0.73])

    x = rnd.uniform(xmin+rmax*arcmin2rad,xmax-rmax*arcmin2rad,ncones)
    y = rnd.uniform(ymin+rmax*arcmin2rad,ymax-rmax*arcmin2rad,ncones)

    matter=False
    hubble=True
    
    if matter:
        cosmologies=Om
    elif hubble:
        cosmologies=H0
    else:
        print "what should I plot?"


    kappacos=numpy.zeros((len(cosmologies),ncones))
    kappacostruth=numpy.zeros(ncones)
    
    for i in range(len(cosmologies)):
        if matter:
            cosmo=[Om[i],1-Om[i],0.73]
        elif hubble:
            cosmo=[0.25,0.75,H0[i]]
        kappa_empty = smooth.smooth(zl,zs,data,truncationscale=truncationscale,hardcut="Rvir",nplanes=50,cosmo=cosmo)

        for j in range(ncones):
            xc = [x[j],y[j]]

            lc = lens_lightcone(catalog,rmax,xc,zl,zs,cosmo=[0.25,0.75,0.73])
            lc.make_kappa_contributions(truncationscale=truncationscale)
            kappacostruth[j]=numpy.sum(lc.galaxies.kappa_keeton)-kappa_empty_true

            lc = lens_lightcone(catalog,rmax,xc,zl,zs,cosmo=cosmo)
            lc.make_kappa_contributions(truncationscale=truncationscale)
            kappacos[i,j]=numpy.sum(lc.galaxies.kappa_keeton)-kappa_empty-kappacostruth[j]


    order=numpy.argsort(kappacostruth)

    ktrue=numpy.take(kappacostruth,order)
    kcos=numpy.take(kappacos,order,axis=1)




    #now make a plot
    plt.figure(figsize=(10,6))
    for j in range(ncones):
            if j < 7: plt.plot(cosmologies,kcos[:,j],label="%.3f"%ktrue[j])
            elif j < 14: plt.plot(cosmologies,kcos[:,j],ls= 'dashed',label="%.3f"%ktrue[j])
            plt.xlim([cosmologies[0],cosmologies[-1]])
    plt.ylabel("$\Delta \kappa_{\mathrm{Keeton}}$")
    if matter:
        plt.xlabel("$\Omega_{\mathrm{M}}$")
        plt.legend(title = "$\kappa_{\mathrm{Keeton}}$$(\Omega_{\mathrm{M}}=0.25)$")
    elif hubble:
        plt.xlabel("$h$")
        plt.legend(title = "$\kappa_{\mathrm{Keeton}}$$(h=0.73)$")
    plt.savefig("kappavsh.png")
    plt.show()


    #print "Computing Keeton (2003) convergence at optical axis due to each halo..."
    #lc.make_kappa_contributions()

    #print "Total external convergence =",numpy.sum(lc.galaxies.kappa_keeton)
    # Now make an illustrative plot:
    


# ============================================================================

if __name__ == '__main__':

# Simulated (Millenium) lightcone data from Hilbert:
    datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"

# First 99 lines of the same, for testing:
#    datafile = "test.txt"
    
    master = atpy.Table(datafile, type='ascii')
    print "Read in master table, length",len(master)
    
    test1(master)

# ============================================================================





"""
kappaN=numpy.ones(n)
z_phot=numpy.ones(n)



n=numpy.size(kappaN)
KappaN=numpy.zeros(n)
KappaNtrue=numpy.zeros(n)

#{ Test 1: compare one line of sight with randomised kappas and redshifts
for i in range(n):
    kappaI=kappaN[i]*(rnd.randn()) #*?
    zN=z_phot[i]*(1+rnd.randn()*.05) #?

    KappaN[i]=kappaI*beta(zl,zN,zs)
    KappaNtrue[i]=kappaN[i]*beta(zl,z_phot[i],zs)



Kappa=KappaN.sum()
KappaTrue=KappaNtrue.sum()

Kappasum=KappaN.cumsum()
KappaTruesum=KappaNtrue.cumsum()
y=range(n)
plt.plot(KappaTruesum,y,'b',linestyle='--')
plt.plot(KappaTruesum,y,'b')
plt.show()
#}



#{ Test2: calculate lots of kappas, and see how we do on average:

kappaTrue=0.0
kappa=0.0
m=1000
Kappa=numpy.zeros(m)
for i in range(n):
    KappaNtrue[i]=kappaN[i]*beta(zl,z_phot[i],zs)
    for j in range(m):
        print m
        zN=z_phot[i]*rnd.randn()*.05 #?
        kappaI=kappaN[i]*rnd.randn()*1 #?
        Kappa[m]+=kappaI*beta(zl,zN,zs)

Kappatrue=KappaNtrue.sum()    

plt.hist(Kappa)
plt.axvline(x=Kappatrue, ymin=0, ymax=1,color='black', ls='dotted')
plt.show()
#}
"""
