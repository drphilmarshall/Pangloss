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

from math import pi

# ===========================================================================
# Get the BoRG pdfs

dropout_fields = np.genfromtxt('../../BORG/LensingModifications/python/data/borg_dropouts.txt', comments = '#', delimiter="\t", usecols=0, dtype='S30')

for i in range(1):
    pdf_file = "../../BORG/LensingModifications/pangloss/figs/borg/"+dropout_fields[i]+"_PofMu.txt"
    borg_pdf = np.genfromtxt(pdf_file, comments = '#')

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
    # We don't know p but it is between 0 and 1. Assume it is uniformly distributed
    p = pymc.Uniform("p", 0, 1)
    
    # Assign data points to clusters using Categorical: value=k-length array of probs
    assignment = pymc.Categorical("assignment", [p, 1 - p], size=borg_pdf.shape[0])
    print "prior assignment, with p = %.2f:" % p.value
    print assignment.value[:10], "..."
    
    # ----------------------------------------------------------------------------
    # tau is 'Precision'. Give the range of standard deviations
    taus = 1.0 / pymc.Uniform("stds", 0., 1.5, size=2) ** 2
    
    # Priors: predict centres and taus (1/std^2)
    centers = pymc.Normal("centers", [0.75, 1.0], [100.0, 16.0], size=2)
    
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
    plt.subplot(311)
    lw = 1
    center_trace = mcmc.trace("centers")[:]
    std_trace = mcmc.trace("stds")[:]    
    p_trace = mcmc.trace("p")[:]
    std_trace = mcmc.trace("stds")[:]

    # for pretty colors later in the book.    
    if center_trace[-1, 0] > center_trace[-1, 1]: colors = ["#348ABD", "#A60628"]
    else: colors = ["#A60628", "#348ABD"]
    
    plt.plot(center_trace[:, 0], label="trace of center 0", c=colors[0], lw=lw)
    plt.plot(center_trace[:, 1], label="trace of center 1", c=colors[1], lw=lw)
    plt.title("Traces of unknown parameters")
    leg = plt.legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)
    
    plt.subplot(312)

    plt.plot(std_trace[:, 0], label="trace of standard deviation of cluster 0",
        c=colors[0], lw=lw)
    plt.plot(std_trace[:, 1], label="trace of standard deviation of cluster 1",
        c=colors[1], lw=lw)
    plt.legend(loc="upper left")
    
    plt.subplot(313)
    plt.plot(p_trace, label="$p$: frequency of assignment to cluster 0",
        color="#467821", lw=lw)
    plt.xlabel("Steps")
    plt.ylim(0, 1)
    plt.legend()
    
   
    # sample again
    mcmc.sample(200000)
    
    plt.figure(figsize=(12.5, 4))
    center_trace = mcmc.trace("centers", chain=1)[:]
    prev_center_trace = mcmc.trace("centers", chain=0)[:]
    
    x = np.arange(50000)
    plt.plot(x, prev_center_trace[:, 0], label="previous trace of center 0",
        lw=lw, alpha=0.4, c=colors[1])
    plt.plot(x, prev_center_trace[:, 1], label="previous trace of center 1",
        lw=lw, alpha=0.4, c=colors[0])
    
    x = np.arange(50000, 250000)
    plt.plot(x, center_trace[:, 0], label="new trace of center 0", lw=lw, c="#348ABD")
    plt.plot(x, center_trace[:, 1], label="new trace of center 1", lw=lw, c="#A60628")
    
    plt.title("Traces of unknown center parameters")
    leg = plt.legend(loc="upper right")
    leg.get_frame().set_alpha(0.8)
    plt.xlabel("Steps")


    # Plot the posteriors for each cluster
    plt.figure()
    _i = [1, 2, 3, 0]
    for i in range(2):
        plt.subplot(2, 2, _i[2 * i])
        plt.title("Posterior of center of cluster %d" % i)
        plt.hist(center_trace[:, i], color=colors[i], bins=30,
                histtype="stepfilled")
    
        plt.subplot(2, 2, _i[2 * i + 1])
        plt.title("Posterior of standard deviation of cluster %d" % i)
        plt.hist(std_trace[:, i], color=colors[i], bins=30,
                histtype="stepfilled")
        # plt.autoscale(tight=True)
    
    plt.tight_layout()
    
    
    plt.figure()
    x = np.linspace(borg_pdf.min(), borg_pdf.max(), 500)
    
    # Plot kde
    IQR = np.percentile(borg_pdf, 75) - np.percentile(borg_pdf, 25)
    Npoints = float(len(borg_pdf))
    print Npoints
    width = 2 * IQR * (Npoints)**(-1./3.)
    Nbins = 6./width
    Nbins = int(round(Nbins,-1))
    print Nbins
    
    kde = stats.kde.gaussian_kde(borg_pdf, bw_method=0.1)#/new_pdf.std(ddof=1))
    kde.covariance_factor = lambda : .04
    kde._compute_covariance()
    plt.plot(x, kde(x), color='k', linestyle='dashed', label="Gaussian kde")
    
    # Plot histogram
    plt.hist(borg_pdf, bins=Nbins, histtype="step", normed=True, color="k",
            lw=2, label="histogram of data")
    
    # Plot the 2 distributions
    posterior_center_means = center_trace.mean(axis=0)
    posterior_std_means = std_trace.mean(axis=0)
    posterior_p_mean = mcmc.trace("p")[:].mean()

    # Cluster 1
    y1 = posterior_p_mean * stats.norm.pdf(x, loc=posterior_center_means[0],
                                    scale=posterior_std_means[0])
    plt.plot(x, y1, color=colors[1], label="Cluster 0", lw=3)
    plt.fill_between(x, y1, color=colors[1], alpha=0.3)
   
    # Cluster 2    
    y2 = (1 - posterior_p_mean) * stats.norm.pdf(x, loc=posterior_center_means[1],
                                        scale=posterior_std_means[1])
    plt.plot(x, y2, color=colors[0], label="Cluster 1", lw=3)
    plt.fill_between(x, y2, color=colors[0], alpha=0.3)
    

    
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
    # Joint pdf
    plt.plot(x, y1+y2, label="Joint pdf", lw=3)                
    plt.legend(loc="upper right")

    plt.title(dropout_fields[i]+": Visualizing Clusters using posterior-mean parameters")
    
    plt.xlim(0.0, 3.0)
    plt.show()