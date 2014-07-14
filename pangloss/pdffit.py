#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ===========================================================================
# Fit BoRG pdfs with 2 gaussians

# http://nbviewer.ipython.org/github/CamDavidsonPilon/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/blob/master/Chapter3_MCMC/IntroMCMC.ipynb

# ===========================================================================

__author__ = 'cmason'

# ===========================================================================
# Import useful packages

import numpy as np
import pymc
import pylab as plt
import scipy.stats as stats
from astropy.io import ascii
import os

from math import pi

# ============================================================================

def fit_gauss(pdf, center1, center2, sample1, burnin, sample2):
    
    """
    We're going to fit with 2 gaussians - one broad, with mean ~ 1, 
    one narrow, with mean <1
    
    We do this by assuming two clusters of the data:
        1) For each data point, it belongs in cluster 1 with probability p
            or in cluster 2 with probability 1-p (cluster number i=0,1)
        2) Draw a random variate from a Normal distribution with parameters 
            μi and σi where i was chosen in step 1
        3) Repeat
    """

    # ===========================================================================
    # MCMC Fit

    # ----------------------------------------------------------------------------
    # We don't know p but it is between 0 and 1. Assume it is uniformly distributed
    p = pymc.Uniform("p", 0, 1)
    
    # Assign data points to clusters using Categorical: value=k-length array of probs
    assignment = pymc.Categorical("assignment", [p, 1 - p], size=borg_pdf.shape[0])
    print "prior assignment, with p = %.2f:" % p.value
    print assignment.value[:10], "..."
    
    # ----------------------------------------------------------------------------
    # tau is 'Precision'. Give the range of standard deviations
    taus = 1.0 / pymc.Uniform("stds", 0, 1 , size=2) ** 2
    #taus = 1.0 / pymc.Normal("stds", [0.2, 0.8], [4.0, 1.0], size=2) ** 2
    
    # Priors: predict centres and taus (1/std^2)
    centers = pymc.Normal("centers", [center1, center2], [25.0, 0.25], size=2)
    
    """
    The below deterministic functions map an assignment, in this case 0 or 1,
    to a set of parameters, located in the (1,2) arrays `taus` and `centers`.
    """
    
    @pymc.deterministic
    def center_i(assignment=assignment, centers=centers):
        return centers[assignment]
    
    @pymc.deterministic
    def tau_i(assignment=assignment, taus=taus):
        return taus[assignment]
    
    print "Random assignments: ", assignment.value[:4], "..."
    print "Assigned center: ", center_i.value[:4], "..."
    print "Assigned precision: ", tau_i.value[:4], "..."
    
    # and to combine it with the observations:
    observations = pymc.Normal("obs", center_i, tau_i, value=borg_pdf, observed=True)
    
    # below we create a model class
    model = pymc.Model([p, assignment, taus, centers])
    
    # Sample the space with 50000 steps
    map_ = pymc.MAP( model )
    map_.fit()
    print "means MAP:",centers.value,"std MAP:",1./np.sqrt(taus.value)

    mcmc = pymc.MCMC(model)
    mcmc.sample(sample1, burnin)
    sample1_ = sample1 - burnin
    
    # ----------------------------------------------------------------------------
    # Plot the traces of the parameters - i.e. do they converge?
    center_trace = mcmc.trace("centers")[:]
    std_trace = mcmc.trace("stds")[:]    
    p_trace = mcmc.trace("p")[:]

    colors = ["#A60628", "#348ABD"]

    # plt.figure(figsize=(12.5, 9))

    # lw = 1

    
    # plt.subplot(311)
    # plt.plot(center_trace[:, 0], label="trace of center 0", c=colors[0], lw=lw)
    # plt.plot(center_trace[:, 1], label="trace of center 1", c=colors[1], lw=lw)
    # plt.title("Traces of unknown parameters")
    # leg = plt.legend(loc="upper right")
    # leg.get_frame().set_alpha(0.7)
    
    # plt.subplot(312)
    # plt.plot(std_trace[:, 0], label="trace of standard deviation of cluster 0",
    #     c=colors[0], lw=lw)
    # plt.plot(std_trace[:, 1], label="trace of standard deviation of cluster 1",
    #     c=colors[1], lw=lw)
    # plt.legend(loc="upper left")
    
    # plt.subplot(313)
    # plt.plot(p_trace, label="$p$: frequency of assignment to cluster 0",
    #     color="#467821", lw=lw)
    # plt.xlabel("Steps")
    # plt.ylim(0, 1)
    # plt.legend()
    
    # ----------------------------------------------------------------------------   
    # sample again
    mcmc.sample(sample2)
    
    # plt.figure(figsize=(12.5, 4))
    # center_trace = mcmc.trace("centers", chain=1)[:]
    # prev_center_trace = mcmc.trace("centers", chain=0)[:]
    
    # x = np.arange(sample1_)
    # plt.plot(x, prev_center_trace[:, 0], label="previous trace of center 0",
    #     lw=lw, alpha=0.4, c=colors[1])
    # plt.plot(x, prev_center_trace[:, 1], label="previous trace of center 1",
    #     lw=lw, alpha=0.4, c=colors[0])
    
    # x = np.arange(sample1_, sample1_ + sample2)
    # plt.plot(x, center_trace[:, 0], label="new trace of center 0", lw=lw, c="#348ABD")
    # plt.plot(x, center_trace[:, 1], label="new trace of center 1", lw=lw, c="#A60628")
    
    # plt.title("Traces of unknown center parameters")
    # leg = plt.legend(loc="upper right")
    # leg.get_frame().set_alpha(0.8)
    # plt.xlabel("Steps")


    # ----------------------------------------------------------------------------   
    # # Plot the posteriors for each cluster
    # plt.figure()
    # _i = [1, 2, 3, 0]
    # for i in range(2):
    #     plt.subplot(2, 2, _i[2 * i])
    #     plt.title("Posterior of center of cluster %d" % i)
    #     plt.hist(center_trace[:, i], color=colors[i], bins=30,
    #             histtype="stepfilled")
    
    #     plt.subplot(2, 2, _i[2 * i + 1])
    #     plt.title("Posterior of standard deviation of cluster %d" % i)
    #     plt.hist(std_trace[:, i], color=colors[i], bins=30,
    #             histtype="stepfilled")
    #     # plt.autoscale(tight=True)
    
 


    # p, std are uniform - so get their means
    # centers are normal, so want most probable value

    posterior_center_means = center_trace.mean(axis=0)
#    posterior_center_mostprob = [np.mean(center_trace[:,0]), np.mean(center_trace[:,1])]
#    posterior_std_means = [np.mean(std_trace[:,0]), np.mean(std_trace[:,1])]
    posterior_std_means = std_trace.mean(axis=0)
    posterior_p_mean = mcmc.trace("p")[:].mean()

    del model
    del mcmc
    del map_
    del taus, assignment, observations, centers

    return posterior_p_mean, posterior_center_means[0], posterior_std_means[0], posterior_center_means[1], posterior_std_means[1]

# ===========================================================================
# Get the BoRG pdfs

borg_field = np.genfromtxt('../BORG/LensingModifications/pangloss/data/borg_overdensity.txt', comments = '#', skip_header=1, usecols=0, dtype='S30')
borg_pdf_table = []

# doagain = ['borg_0240-1857','borg_0952+5304', 'borg_1031+5052', 'borg_1059+0519', 
#             'borg_1119+4026', 'borg_1301+0000', 'borg_1358+4326', '1408+5503']
doagain = ['borg_1059+0519']


# 12 at a time
for i in range(1):
    field = "All BoRG fields"
    pdf_file = 

    # field = doagain[i] #borg_field[i+60]
    # pdf_file = "../BORG/LensingModifications/pangloss/figs/borg/"+field+"_PofMu.txt"

    borg_pdf = np.genfromtxt(pdf_file, comments = '#')

    mask = np.where(borg_pdf >= 0.)
    borg_pdf = borg_pdf[mask]

    mode = stats.mode(borg_pdf)[0]
    # ----------------------------------------------------------------------------   
    # Do the fit
    p, mean1, std1, mean2, std2 = fit_gauss(pdf=borg_pdf, center1 = mode, center2 = 1.2, sample1=50000, burnin=40000, sample2=200000)


    # Save the parameters to a file

    table = [field, p, mean1, std1, mean2, std2 ]
    borg_pdf_table.append(table)
    
    # plt.clf()

    # del table

borg_pdf_table = np.array(borg_pdf_table)
ascii.write(borg_pdf_table, 'borg_pdf_table_again2.txt', names=['Field', 'p', 'mean1', 'std1', 'mean2', 'std2'])




