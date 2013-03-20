import lightcone,kappamap
import cPickle
import atpy
import numpy.random as rnd
import numpy
#==========================================================================
arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad
#==========================================================================
def MakeCones(N,savefolder,radius=5,catalogue=[],kappa=[],gamma1=[],gamma2=[]):
    if len(catalogue)!=len(kappa): 
        print "the catalogues and the ray-traced results need to match up!"
	exit()
    Ncones=N/len(catalogue)
    l=0
    for i in range(len(catalogue)):
	print "reading in catalogue"
	master = atpy.Table(catalogue[i], type='ascii')

        xmax = master['pos_0[rad]'].max()
        xmin = master['pos_0[rad]'].min()
        ymax = master['pos_1[rad]'].max()
        ymin = master['pos_1[rad]'].min()

        
        print "reading in kappamap"
        kappafile = kappa[i]
        MSconvergence = kappamap.Kappamap(kappafile)
        gammafile1 = gamma1[i]
        MSgamma1 = kappamap.Kappamap(gammafile1)
        gammafile2 = gamma2[i]
        MSgamma2 = kappamap.Kappamap(gammafile2)

        x = rnd.uniform(xmin+Rcone*arcmin2rad,xmax-Rcone*arcmin2rad,Ncones)
        y = rnd.uniform(ymin+Rcone*arcmin2rad,ymax-Rcone*arcmin2rad,Ncones)

        print "Writing out cones"
        for k in range(Ncones):
            if k % 200 == 0 and k !=0: print ("saving cone %i of %i" %(k,Ncones))
            xc = [x[k],y[k]]
            lc= lightcone.Lightcone(master,xc,Rcone)
            lc.kappa_hilbert=MSconvergence.at(x[k],y[k],coordinate_system='physical')
            lc.gamma1_hilbert=MSgamma1.at(x[k],y[k],coordinate_system='physical')
            lc.gamma2_hilbert=MSgamma2.at(x[k],y[k],coordinate_system='physical')
            CONE=open("%s/cone%i.cone"%(savefolder,l),"wb")
            cPickle.dump(lc,CONE,protocol=2)
            CONE.close()
            l+=1
            
