"""
NAME
    nocorr

PURPOSE
    Dummy alternative to Mike Jarvis' treecorr module, if it is not yet
    installed. All classes return Null objects.

COMMENTS

INITIALISATION

CLASSES
    Catalog(self)
    GGCorrelation(self)
    NGCorrelation(self)

BUGS

AUTHORS
    This file is part of the Pangloss project, distributed under the
    MIT License, by Tom Collett (IoA) and  Phil Marshall (Oxford).
    Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

HISTORY
    2015-08-05  Started Marshall (SLAC)
"""
# ----------------------------------------------------------------------------

class Catalog(object):
    def __init__(self,**kwargs):
        return None
    def __str__(self):
        return 'Null object, returned in place of a treecorr.Catalog.'

class GGCorrelation(object):
    def __init__(self,**kwargs):
        return None
    def __str__(self):
        return 'Null object, returned in place of a treecorr.GGCorrelation.'

class NGCorrelation(object):
    def __init__(self,**kwargs):
        return None
    def __str__(self):
        return 'Null object, returned in place of a treecorr.NGCorrelation.'

# ----------------------------------------------------------------------------
