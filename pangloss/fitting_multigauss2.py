import pymc as mc
import numpy as np
import pylab as plt
def GaussFunc(x, amplitude, centroid, sigma):
    return amplitude * np.exp(-0.5 * ((x - centroid) / sigma)**2)

# points = np.arange(5000, 5050, 0.02)

# # Profile 1
# centroid_one = 5025.0
# sigma_one = 2.2
# height_one = 0.8
# profile1 = GaussFunc(points, height_one, centroid_one, sigma_one, )

# # Profile 2
# centroid_two = 5027.0
# sigma_two = 1.2
# height_two = 0.5
# profile2 = GaussFunc(points, height_two, centroid_two, sigma_two, )

# # Measured values
# noise = np.random.normal(0.0, 0.02, len(points))
# combined = profile1 + profile2 + noise

# print type(combined), combined.shape

dropout_fields = np.genfromtxt('../BORG/LensingModifications/python/data/borg_dropouts.txt', comments = '#', delimiter="\t", usecols=0, dtype='S30')

pdf_file = "../BORG/LensingModifications/pangloss/figs/borg/"+dropout_fields[0]+"_PofMu.txt"
borg_pdf = np.genfromtxt(pdf_file, comments = '#')
borg_pdf = np.array(borg_pdf)
print type(borg_pdf), borg_pdf.shape
#points = np.linspace(borg_pdf.min(), borg_pdf.max(), 2000)
points = np.linspace(0, 5, len(borg_pdf))
print len(points)

# If you want to plot what this looks like
plt.clf()
# plt.figure()
# plt.plot(points, combined, label="Measured")
# plt.plot(points, profile1, color='red', linestyle='dashed', label="1")
# plt.plot(points, profile2, color='green', linestyle='dashed', label="2")
# plt.title("Feature One and Two")
# plt.legend()

### Suggested one runs the above code first.
### Unknowns we are interested in

est_centroid_one = mc.Uniform("est_centroid_one", 0.6, 1.4 )
est_centroid_two = mc.Uniform("est_centroid_two", 0.6, 1.8 )

est_sigma_one = mc.Uniform( "est_sigma_one", 0, 1 )
est_sigma_two = mc.Uniform( "est_sigma_two", 0, 2 )

est_height_one = mc.Uniform( "est_height_one", 0, 2 ) 
est_height_two = mc.Uniform( "est_height_two", 0, 1 ) 

#std deviation of the noise, converted to precision by tau = 1/sigma**2
precision= 1./mc.Uniform("std", 0, 2)**2

#Set up the model's relationships.

@mc.deterministic( trace = False) 
def est_profile_1(x = points, centroid = est_centroid_one, sigma = est_sigma_one, height= est_height_one):
    return GaussFunc( x, height, centroid, sigma )


@mc.deterministic( trace = False) 
def est_profile_2(x = points, centroid = est_centroid_two, sigma = est_sigma_two, height= est_height_two):
    return GaussFunc( x, height, centroid, sigma )


@mc.deterministic( trace = False )
def mean( profile_1 = est_profile_1, profile_2 = est_profile_2 ):
    return profile_1 + profile_2


observations = mc.Normal("obs", mean, precision, value = borg_pdf, observed = True)


model = mc.Model([est_centroid_one, 
              est_centroid_two, 
                est_height_one,
                est_height_two,
                est_sigma_one,
                est_sigma_two,
                precision])

#always a good idea to MAP it prior to MCMC, so as to start with good initial values
map_ = mc.MAP( model )
map_.fit()

mcmc = mc.MCMC( model )
mcmc.sample(50000) 



# ----------------------------------------------------------------------------
# Plot the traces of the parameters - i.e. do they converge?
plt.figure(figsize=(12.5, 9))
plt.subplot(311)
lw = 1
center1_trace = mcmc.trace("est_centroid_one")[:]
center2_trace = mcmc.trace("est_centroid_two")[:]

std1_trace = mcmc.trace("est_sigma_one")[:]
std2_trace = mcmc.trace("est_sigma_two")[:]

height1_trace = mcmc.trace("est_height_one")[:]
height2_trace = mcmc.trace("est_height_two")[:]

# for pretty colors later in the book.    
if center1_trace[-1] > center2_trace[-1]: colors = ["#348ABD", "#A60628"]
else: colors = ["#A60628", "#348ABD"]

plt.plot(center1_trace[:], label="trace of center 0", c=colors[0], lw=lw)
plt.plot(center2_trace[:], label="trace of center 1", c=colors[1], lw=lw)
plt.title("Traces of unknown parameters")
leg = plt.legend(loc="upper right")
leg.get_frame().set_alpha(0.7)

plt.subplot(312)

plt.plot(std1_trace[:], label="trace of standard deviation of cluster 0",
    c=colors[0], lw=lw)
plt.plot(std2_trace[:], label="trace of standard deviation of cluster 1",
    c=colors[1], lw=lw)
plt.legend(loc="upper left")

plt.subplot(313)
plt.plot(height1_trace, label="trace of height of cluster 0",
    c=colors[0], lw=lw)
plt.plot(height2_trace, label="trace of height of cluster 1",
    c=colors[1], lw=lw)
plt.xlabel("Steps")
plt.ylim(0, 1)
plt.legend()



# sample again
mcmc.sample(100000)

plt.figure(figsize=(12.5, 4))
center1_trace = mcmc.trace("est_centroid_one", chain=1)[:]
center2_trace = mcmc.trace("est_centroid_two", chain=1)[:]

std1_trace = mcmc.trace("est_sigma_one", chain=1)[:]
std2_trace = mcmc.trace("est_sigma_two", chain=1)[:]

height1_trace = mcmc.trace("est_height_one", chain=1)[:]
height2_trace = mcmc.trace("est_height_two", chain=1)[:]

prev_center1_trace = mcmc.trace("est_centroid_one", chain=0)[:]
prev_center2_trace = mcmc.trace("est_centroid_two", chain=0)[:]

x = np.arange(50000)
plt.plot(x, prev_center1_trace[:], label="previous trace of center 0",
    lw=lw, alpha=0.4, c=colors[0])
plt.plot(x, prev_center2_trace[:], label="previous trace of center 1",
    lw=lw, alpha=0.4, c=colors[1])

x = np.arange(50000, 150000)
plt.plot(x, center1_trace[:], label="new trace of center 0", lw=lw, c=colors[0])
plt.plot(x, center2_trace[:], label="new trace of center 1", lw=lw, c=colors[1])

plt.title("Traces of unknown center parameters")
leg = plt.legend(loc="upper right")
leg.get_frame().set_alpha(0.8)
plt.xlabel("Steps")


# Plot the posteriors for each cluster
plt.figure()
_i = [1, 2, 3, 0]
centers = [center1_trace, center2_trace]
stds =[std1_trace, std2_trace]

for i in range(2):
    plt.subplot(2, 2, _i[2 * i])
    plt.title("Posterior of center of cluster %d" % i)
    plt.hist(centers[i], color=colors[i], bins=30,
            histtype="stepfilled")

    plt.subplot(2, 2, _i[2 * i + 1])
    plt.title("Posterior of standard deviation of cluster %d" % i)
    plt.hist(stds[i], color=colors[i], bins=30,
            histtype="stepfilled")
    # plt.autoscale(tight=True)

plt.tight_layout()

plt.figure()

IQR = np.percentile(borg_pdf, 75) - np.percentile(borg_pdf, 25)
Npoints = float(len(borg_pdf))
width = 2 * IQR * (Npoints)**(-1./3.)
Nbins = 6./width
Nbins = int(round(Nbins,-1))
print Nbins
plt.hist(borg_pdf, bins=Nbins, histtype="step", normed=True, color="k",
            lw=2, label="histogram of data")

# Plot the 2 distributions
posterior_center_means = [center1_trace.mean(axis=0), center2_trace.mean(axis=0)]
posterior_std_means = [std1_trace.mean(axis=0), std2_trace.mean(axis=0)]
posterior_height_means = [height1_trace.mean(axis=0), height2_trace.mean(axis=0)]

# Cluster 1
y1 = GaussFunc(points, posterior_height_means[0], posterior_center_means[0], posterior_std_means[0])
plt.plot(points, y1, color=colors[0], label="Cluster 0", lw=3)
plt.fill_between(points, y1, color=colors[0], alpha=0.3)

# Cluster 2    
y2 = GaussFunc(points, posterior_height_means[1], posterior_center_means[1], posterior_std_means[1])
plt.plot(points, y2, color=colors[1], label="Cluster 1", lw=3)
plt.fill_between(points, y2, color=colors[1], alpha=0.3)

plt.show()

