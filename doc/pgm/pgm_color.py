# ============================================================================
# PGM for the SDSI proposal
#
# Phil Marshall, September 2014
# ============================================================================

from matplotlib import rc
rc("font", family="serif", size=10)
rc("text", usetex=True)

import daft

figshape = (8.6, 5.0)
figorigin = (-0.6, -0.4)

# Colors.
red = {"ec": "red"}
orange = {"ec": "orange"}
green = {"ec": "green"}
blue = {"ec": "blue"}
violet = {"ec": "violet"}

# Start the chart:

pgm = daft.PGM(figshape, origin=figorigin)

# Foreground galaxies branch:
pgm.add_node(daft.Node("cosmo1", r"${\bf \Omega}_g$", 1.9, 3.8,plot_params=violet)) 
pgm.add_node(daft.Node("Mhd", r"$M_{{\rm h},i}$", 1.3, 2.2,plot_params=green))
pgm.add_node(daft.Node("zd", r"$z_i$", 2.2, 2.8,plot_params=green))
pgm.add_node(daft.Node("xd", r"${\bf x}_i$", 2.2, 0.7, fixed=True)) 
pgm.add_node(daft.Node("Mstard", r"$M^{*}_i$", 2.2, 1.6,plot_params=orange)) 
pgm.add_node(daft.Node("phid", r"${\bf \phi}_i$", 1.3, 1.1, observed=True))
pgm.add_node(daft.Node("alphad", r"$\alpha$", 0.2, 1.6,plot_params=red))

# Background galaxies:
pgm.add_node(daft.Node("cosmo2", r"${\bf \Omega}_e$", 3.0, 4.0,plot_params=violet)) 
pgm.add_node(daft.Node("sigcrit", r"$\Sigma^{\rm crit}_{ij}$", 3.0, 2.4,plot_params=blue))
pgm.add_node(daft.Node("gamma", r"$\gamma_j$", 4.0, 1.8,plot_params=blue))
pgm.add_node(daft.Node("elens", r"$\epsilon^{\rm lens}_j$", 4.5, 1.2,plot_params=blue))
pgm.add_node(daft.Node("eobs", r"$\epsilon^{\rm obs}_j$", 4.5, 0.4, observed=True))

pgm.add_node(daft.Node("Mh", r"$M_{{\rm h},j}$", 6.0, 3.0,plot_params=green))
pgm.add_node(daft.Node("z", r"$z_j$", 5.0, 3.0,plot_params=green))
pgm.add_node(daft.Node("x", r"${\bf x}_j$", 4.6, 2.4, fixed=True))
pgm.add_node(daft.Node("cosmo3", r"${\bf \Omega}_g$", 5.6, 4.0,plot_params=violet)) 
pgm.add_node(daft.Node("Mstar", r"$M^{*}_j$", 5.6, 2.4,plot_params=orange)) 
pgm.add_node(daft.Node("alpha", r"${\bf \alpha}$", 7.2, 2.4,plot_params=red)) 
pgm.add_node(daft.Node("epsilon", r"$\epsilon_j$", 6.0, 1.4))
pgm.add_node(daft.Node("phi", r"${\bf \phi}_j$", 5.0, 1.8, observed=True))

# Now connect the dots:
pgm.add_edge("cosmo1", "Mhd")
pgm.add_edge("cosmo1", "zd")
pgm.add_edge("Mhd", "Mstard")
pgm.add_edge("zd", "Mstard")
pgm.add_edge("alphad", "Mstard")
pgm.add_edge("zd", "phid")
pgm.add_edge("Mstard", "phid")
pgm.add_edge("zd", "sigcrit")
pgm.add_edge("cosmo2", "sigcrit")
pgm.add_edge("z", "sigcrit")
pgm.add_edge("xd", "gamma")
pgm.add_edge("Mhd", "gamma")
pgm.add_edge("sigcrit", "gamma")
pgm.add_edge("x", "gamma")
pgm.add_edge("gamma", "elens")
pgm.add_edge("elens", "eobs")
pgm.add_edge("cosmo3", "Mh")
pgm.add_edge("cosmo3", "z")
pgm.add_edge("Mh", "Mstar")
pgm.add_edge("z", "Mstar")
pgm.add_edge("alpha", "Mstar")
pgm.add_edge("Mstar", "epsilon")
pgm.add_edge("Mstar", "phi")
pgm.add_edge("z", "phi")
pgm.add_edge("epsilon", "elens")


#Add plate over deflectors
pgm.add_plate(daft.Plate([0.7, 0.2, 2.7, 3.2], label=r"foreground galaxies $i$", label_offset=[5, 5]))

# Add plate over sources
pgm.add_plate(daft.Plate([2.6, 0.0, 4.0, 3.6], label=r"background galaxies $j$", label_offset=[120, 5]))


pgm.render()
pgm.figure.savefig("pgm_color.png", dpi=220)

# ============================================================================

