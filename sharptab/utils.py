''' Frequently used functions '''
import numpy as np
from numba import njit
from sharptab.constants import MISSING, TOL

@njit
def mag(u, v, missing=MISSING):
    '''
    Compute the magnitude of a vector from its components

    Parameters
    ----------
    u : number, array_like
        U-component of the wind
    v : number, array_like
        V-component of the wind
    missing : number (optional)
        Optional missing parameter. If not given, assume default missing
        value from sharppy.sharptab.constants.MISSING

    Returns
    -------
    mag : number, array_like
        The magnitude of the vector (units are the same as input)

    '''
    return np.sqrt(u**2 + v**2)

@njit
def vec2comp(wdir, wspd, missing=MISSING):
    '''
    Convert direction and magnitude into U, V components

    Parameters
    ----------
    wdir : number, array_like
        Angle in meteorological degrees
    wspd : number, array_like
        Magnitudes of wind vector (input units == output units)
    missing : number (optional)
        Optional missing parameter. If not given, assume default missing
        value from sharppy.sharptab.constants.MISSING

    Returns
    -------
    u : number, array_like (same as input)
        U-component of the wind (units are the same as those of input speed)
    v : number, array_like (same as input)
        V-component of the wind (units are the same as those of input speed)

    '''
    #if not QC(wdir) or not QC(wspd):
    #    return ma.masked, ma.masked

    #wdir = ma.asanyarray(wdir).astype(np.float64)
    #wspd = ma.asanyarray(wspd).astype(np.float64)
    #wdir.set_fill_value(missing)
    #wspd.set_fill_value(missing)
    #assert wdir.shape == wspd.shape, 'wdir and wspd have different shapes'
    #if wdir.shape:
        #wdir[wdir == missing] = ma.masked
        #wspd[wspd == missing] = ma.masked
        #wdir[wspd.mask] = ma.masked
        #wspd[wdir.mask] = ma.masked
    u, v = _vec2comp(wdir, wspd)
    u[np.fabs(u) < TOL] = 0.
    v[np.fabs(v) < TOL] = 0.
    #else:
    #    if wdir == missing:
    #        wdir = ma.masked
    #        wspd = ma.masked
    #    elif wspd == missing:
    #        wdir = ma.masked
    #        wspd = ma.masked
    #    u, v = _vec2comp(wdir, wspd)
    #    if ma.fabs(u) < TOL:
    #        u = 0.
    #    if ma.fabs(v) < TOL:
    #        v = 0.
    return u, v

@njit
def _vec2comp(wdir, wspd):
    '''
    Underlying function that converts a vector to its components

    Parameters
    ----------
    wdir : number, masked_array
        Angle in meteorological degrees
    wspd : number, masked_array
        Magnitudes of wind vector

    Returns
    -------
    u : number, masked_array (same as input)
        U-component of the wind
    v : number, masked_array (same as input)
        V-component of the wind

    '''
    u = wspd * np.sin(np.radians(wdir)) * -1
    v = wspd * np.cos(np.radians(wdir)) * -1
    return u, v
