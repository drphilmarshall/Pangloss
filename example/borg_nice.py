#!/usr/bin/env python
# ===========================================================================
# Outputs BoRG catalogs in better format for Pangloss
# ===========================================================================

__author__ = 'cmason'
# ===========================================================================
# Import useful packages

import os
import mpmath
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy import integrate
from matplotlib import rc_file
rc_file('matplotlibrc')

# ===========================================================================
catalog = np.genfromtxt('borg_0110-0224_multiband.cat')
