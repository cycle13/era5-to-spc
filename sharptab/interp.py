from numba import njit
import numpy as np

@njit
def components(prof, p):
    '''
    Interpolates the given data to calculate the U and V components at a
    given pressure

    Parameters
    ----------
    prof : profile object
        Profile object
    p : number, numpy array
        Pressure (hPa) of a level

    Returns
    -------
    U and V components at the given pressure (kts) : number, numpy array
    '''
    # Note: numpy's interpoloation routine expects the interpoloation
    # routine to be in ascending order. Because pressure decreases in the
    # vertical, we must reverse the order of the two arrays to satisfy
    # this requirement.
    U = generic_interp_pres(np.log10(p), prof.logp[::-1], prof.u[::-1])
    V = generic_interp_pres(np.log10(p), prof.logp[::-1], prof.v[::-1])
    return U, V

@njit
def to_agl(prof, h):
    '''
    Convert a height from mean sea-level (MSL) to above ground-level (AGL)

    Parameters
    ----------
    h : number, numpy array
        Height of a level
    prof : profile object
        Profile object

    Returns
    -------
    Converted height (m AGL) : number, numpy array

    '''
    return h - prof.hght[prof.sfc]

@njit
def dwpt(prof, p):
    '''
    Interpolates the given data to calculate a dew point temperature
    at a given pressure

    Parameters
    ----------
    prof : profile object
        Profile object
    p : number, numpy array
        Pressure (hPa) of the level for which dew point temperature is desired

    Returns
    -------
    Dew point tmperature (C) at the given pressure : number, numpy array

    '''
    # Note: numpy's interpoloation routine expects the interpoloation
    # routine to be in ascending order. Because pressure decreases in the
    # vertical, we must reverse the order of the two arrays to satisfy
    # this requirement.
    return generic_interp_pres(np.log10(p), prof.logp[::-1], prof.dwpc[::-1])

@njit
def temp(prof, p):
    '''
    Interpolates the given data to calculate a temperature at a given pressure

    Parameters
    ----------
    prof : profile object
        Profile object
    p : number, numpy array
        Pressure (hPa) of the level for which temperature is desired

    Returns
    -------
    Temperature (C) at the given pressure : number, numpy array

    '''
    # Note: numpy's interpoloation routine expects the interpoloation
    # routine to be in ascending order. Because pressure decreases in the
    # vertical, we must reverse the order of the two arrays to satisfy
    # this requirement.
    return generic_interp_pres(np.log10(p), prof.logp[::-1], prof.tmpc[::-1])

@njit
def hght(prof, p):
    '''
    Interpolates the given data to calculate a height at a given pressure

    Parameters
    ----------
    prof : profile object
        Profile object
    p : number, numpy array
        Pressure (hPa) of the level for which height is desired

    Returns
    -------
    Height (m) at the given pressure : number, numpy array

    '''
    # Note: numpy's interpoloation routine expects the interpoloation
    # routine to be in ascending order. Because pressure decreases in the
    # vertical, we must reverse the order of the two arrays to satisfy
    # this requirement.
    return generic_interp_pres(np.log10(p), prof.logp[::-1], prof.hght[::-1])

@njit
def vtmp(prof, p):
    '''
    Interpolates the given data to calculate a virtual temperature
    at a given pressure

    Parameters
    ----------
    prof : profile object
        Profile object
    p : number, numpy array
        Pressure (hPa) of the level for which virtual temperature is desired

    Returns
    -------
    Virtual tmperature (C) at the given pressure : number, numpy array

    '''
    return generic_interp_pres(np.log10(p), prof.logp[::-1], prof.vtmp[::-1])

@njit
def generic_interp_pres(p, pres, field):
    '''
    Generic interpolation routine

    Parameters
    ----------
    p : number, numpy array
        Pressure (hPa) of the level for which the field variable is desired
    pres : numpy array
        The array of pressure
    field : numpy array
        The variable which is being interpolated
    log : bool
        Flag to determine whether the 'field' variable is in log10 space

    Returns
    -------
    Value of the 'field' variable at the given pressure : number, numpy array

    '''
    not_masked1 = np.ones(pres.shape)
    not_masked1[:] = True

    not_masked2 = np.ones(field.shape)
    not_masked2[:] = True
    not_masked = not_masked1 * not_masked2
    field_intrp = np.interp(p, pres, field)
    """
    if hasattr(p, 'shape') and p.shape == tuple():
        p = p[()]

    if type(p) != type(ma.masked) and np.all(~np.isnan(p)):
        # Bug fix for Numpy v1.10: returns nan on the boundary.
        field_intrp = ma.where(np.isclose(p, pres[not_masked][0]), field[not_masked][0], field_intrp)
        field_intrp = ma.where(np.isclose(p, pres[not_masked][-1]), field[not_masked][-1], field_intrp)

    # Another bug fix: np.interp() returns masked values as nan. We want ma.masked, dangit!
    field_intrp = ma.where(np.isnan(field_intrp), ma.masked, field_intrp)

    # ma.where() returns a 0-d array when the arguments are floats, which confuses subsequent code.
    if hasattr(field_intrp, 'shape') and field_intrp.shape == tuple():
        field_intrp = field_intrp[()]
    """

    return field_intrp