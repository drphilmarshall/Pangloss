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

# ===========================================================================
# Get the BoRG pdfs

dropout_fields = np.genfromtxt('../BORG/LensingModifications/python/data/borg_dropouts.txt', comments = '#', delimiter="\t", usecols=0, dtype='S30')
borg_pdf_table = []

for i in range(1):
    pdf_file = "../BORG/LensingModifications/pangloss/figs/borg/"+dropout_fields[i]+"_PofMu.txt"
    borg_pdf = np.genfromtxt(pdf_file, comments = '#')

    mask = np.where(borg_pdf >= -2.)
    borg_pdf = borg_pdf[mask]

    # ===========================================================================
    # MCMC Fit
    
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
    # ----------------------------------------------------------------------------
    n = 3

    # We don't know p but it is between 0 and 1. Assume it is uniformly distributed
    p = pymc.Dirichlet("p", [1,1,1])
    
    # Assign data points to clusters using Categorical: value=k-length array of probs
    assignment = pymc.Categorical("assignment", p=p, size=borg_pdf.shape[0])
    print "prior assignment, with p = %.2f:" % p[0].value
    print assignment.value[:10], "..."
    
    # ----------------------------------------------------------------------------
    # tau is 'Precision'. Give the range of standard deviations
    taus = 1.0 / pymc.Uniform("stds", 0., 5., size=3) ** 2
    
    # Priors: predict centres and taus (1/std^2)
    centers = pymc.Normal("centers", [0.8, 1.2, 1.8], [25.0, 2., 0.25], size=3)
    
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
    mcmc = pymc.MCMC(model)
    mcmc.sample(50000)
    
    # ----------------------------------------------------------------------------
    # Plot the traces of the parameters - i.e. do they converge?
    plt.figure(figsize=(12.5, 9))

    lw = 1
    center_trace = mcmc.trace("centers")[:]
    std_trace = mcmc.trace("stds")[:]    
    p_trace = mcmc.trace("p")[:]

    # for pretty colors later in the book.    
    colors = ["#A60628", "#348ABD", "orange"]
    
    plt.subplot(311)
    plt.plot(center_trace[:, 0], label="trace of center 0", c=colors[0], lw=lw)
    plt.plot(center_trace[:, 1], label="trace of center 1", c=colors[1], lw=lw)
    plt.plot(center_trace[:, 2], label="trace of center 2", c=colors[2], lw=lw)
    plt.title("Traces of unknown parameters")
    leg = plt.legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)
    
    plt.subplot(312)
    plt.plot(std_trace[:, 0], label="trace of standard deviation of cluster 0",
        c=colors[0], lw=lw)
    plt.plot(std_trace[:, 1], label="trace of standard deviation of cluster 1",
        c=colors[1], lw=lw)
    plt.plot(std_trace[:, 2], label="trace of standard deviation of cluster 2",
        c=colors[2], lw=lw)
    plt.legend(loc="upper left")
    
    plt.subplot(313)
    plt.plot(p_trace[:, 0], label="$p$: frequency of assignment to cluster 0",
        c=colors[0], lw=lw)
    plt.plot(p_trace[:, 1], label="$p$: frequency of assignment to cluster 1",
        c=colors[1], lw=lw)
    plt.plot(1-p_trace[:, 0]-p_trace[:, 1], label="$p$: frequency of assignment to cluster 2",
        c=colors[2], lw=lw)
    plt.xlabel("Steps")
    plt.ylim(0, 1)
    plt.legend()
    
    # ----------------------------------------------------------------------------   
    # sample again
    mcmc.sample(100000)
    
    plt.figure(figsize=(12.5, 4))
    center_trace = mcmc.trace("centers", chain=1)[:]
    prev_center_trace = mcmc.trace("centers", chain=0)[:]
    
    x = np.arange(50000)
    plt.plot(x, prev_center_trace[:, 0], label="previous trace of center 0",
        lw=lw, alpha=0.4, c=colors[0])
    plt.plot(x, prev_center_trace[:, 1], label="previous trace of center 1",
        lw=lw, alpha=0.4, c=colors[1])
    plt.plot(x, prev_center_trace[:, 2], label="previous trace of center 2",
        lw=lw, alpha=0.4, c=colors[2])

    x = np.arange(50000, 150000)
    plt.plot(x, center_trace[:, 0], label="new trace of center 0", lw=lw, c=colors[0])
    plt.plot(x, center_trace[:, 1], label="new trace of center 1", lw=lw, c=colors[1])
    plt.plot(x, center_trace[:, 2], label="new trace of center 2", lw=lw, c=colors[2])
    
    plt.title("Traces of unknown center parameters")
    leg = plt.legend(loc="upper right")
    leg.get_frame().set_alpha(0.8)
    plt.xlabel("Steps")


    # ----------------------------------------------------------------------------   
    # Plot the posteriors for each cluster
    plt.figure()
    _i = [1, 2, 3, 4, 5, 6, 0]
    for i in range(3):
        plt.subplot(2, 3, _i[2 * i])
        plt.title("Posterior of center of cluster %d" % i)
        plt.hist(center_trace[:, i], color=colors[i], bins=30,
                histtype="stepfilled")
    
        plt.subplot(2, 3, _i[2 * i + 1])
        plt.title("Posterior of standard deviation of cluster %d" % i)
        plt.hist(std_trace[:, i], color=colors[i], bins=30,
                histtype="stepfilled")
        # plt.autoscale(tight=True)
    
    plt.tight_layout()
 

    # ----------------------------------------------------------------------------   
    # Plot the final pdfs

    plt.figure()
    x = np.linspace(borg_pdf.min(), borg_pdf.max(), 500)
    
    # Plot kde   
    kde = stats.kde.gaussian_kde(borg_pdf)
    plt.plot(x, kde(x), color='k', linestyle='dashed', label="Gaussian kde")
    
    # Plot histogram
    plt.hist(borg_pdf, bins=20, histtype="step", normed=True, color="k",
             lw=2, label="histogram of data")
    
    # p, std are uniform - so get their means
    # centers are normal, so want most probable value

    # Plot the 2 distributions MEANS
    posterior_center_means = center_trace.mean(axis=0)
 #   posterior_center_mostprob = [np.mean(center_trace[:,0]), np.mean(center_trace[:,1])]
    posterior_std_means = std_trace.mean(axis=0)
    posterior_p_mean = p_trace.mean(axis=0)

    # Cluster 1
    y1 = posterior_p_mean[0] * stats.norm.pdf(x, loc=posterior_center_means[0],
                                    scale=posterior_std_means[0])
    plt.plot(x, y1, color=colors[0], label="Cluster 0", lw=3)
    plt.fill_between(x, y1, color=colors[0], alpha=0.3)
   
    # Cluster 2    
    y2 = posterior_p_mean[1] * stats.norm.pdf(x, loc=posterior_center_means[1],
                                        scale=posterior_std_means[1])
    plt.plot(x, y2, color=colors[1], label="Cluster 1", lw=3)
    plt.fill_between(x, y2, color=colors[1], alpha=0.3)

    # Cluster 3    
    y3 = (1-posterior_p_mean[0]-posterior_p_mean[1]) * stats.norm.pdf(x, loc=posterior_center_means[2],
                                        scale=posterior_std_means[2])
    plt.plot(x, y2, color=colors[2], label="Cluster 2", lw=3)
    plt.fill_between(x, y2, color=colors[2], alpha=0.3)

    # Joint pdf
    plt.plot(x, y1+y2+y3, label="Joint pdf", lw=3)                

    table = [dropout_fields[i], posterior_p_mean, posterior_center_means[0], posterior_std_means[0],
             posterior_center_means[1], posterior_std_means[1], posterior_center_means[2], posterior_std_means[2]]

    borg_pdf_table.append(table)

# =======================================================================
    # Cluster 1
   # y1 = 0.8 * stats.norm.pdf(x, loc=0.75,
   #                                 scale=0.1)
   # plt.plot(x, y1, label="Cluster 0", lw=3)
   # plt.fill_between(x, y1, color=colors[1], alpha=0.3)
   #
   # # Cluster 2    
   # y2 = (1 - posterior_p_mean) * stats.norm.pdf(x, loc=1.0,
   #                                     scale=0.25)
   # plt.plot(x, y2, label="Cluster 1", lw=3)
   # plt.fill_between(x, y2, color=colors[0], alpha=0.3)    
   #     
    
    plt.legend(loc="upper right")

    plt.title(dropout_fields[i]+": Visualizing Clusters using posterior-mean parameters")
    
 #   plt.xlim(0.0, 3.0)
    
    plt.tight_layout()
    savedfile = dropout_fields[i]+"_PofMu_MCMC_3mean.pdf"
    plt.savefig(savedfile,dpi=300)
    print "Plot saved as "+os.getcwd()+"/"+savedfile

    plt.show()

borg_pdf_table = np.array(borg_pdf_table)
ascii.write(borg_pdf_table, 'borg_pdf_table.txt', names=['Field', 'p', 'mean1', 'std1', 'mean2', 'std2', 'mean3', 'std3'])
