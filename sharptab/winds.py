''' Wind Manipulation Routines '''
import numpy as np
from numba import njit
from sharptab import interp

@njit
def wind_shear(prof, pbot=850, ptop=250):
    '''
    Calculates the shear between the wind at (pbot) and (ptop).

    Parameters
    ----------
    prof: profile object
        Profile object
    pbot : number (optional; default 850 hPa)
        Pressure of the bottom level (hPa)
    ptop : number (optional; default 250 hPa)
        Pressure of the top level (hPa)

    Returns
    -------
    shu : number
        U-component (kts)
    shv : number
        V-component (kts)

    '''
    ubot, vbot = interp.components(prof, pbot)
    utop, vtop = interp.components(prof, ptop)
    shu = utop - ubot
    shv = vtop - vbot
    return shu, shv
