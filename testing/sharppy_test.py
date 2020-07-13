import numpy as np
import xarray as xr
import sharppy.sharptab.profile as SHARPPY_PROFILE
import sharppy.sharptab.params as SHARPPY_PARAMS
from sharptab.constants import G, ZEROCNK

sat_pressure_0c = 6.112
ZEROCNK = 273.15        # Zero C to K
eps = 0.62197           
ROCP = 0.28571426       # Rd over Cp (Kappa)
c1 = 0.0498646455 ; c2 = 2.4082965 ; c3 = 7.07475
c4 = 38.9114 ; c5 = 0.0915 ; c6 = 1.2035
eps = 0.62197

def wind_direction(u, v):
    r"""Compute the wind direction from u and v-components.
    """
    wdir = 90. - np.arctan2(-v, -u)
    origshape = wdir.shape

    wdir[wdir <= 0] += 360. 
    return wdir.reshape(origshape)

def wind_speed(u, v):
    r"""Compute the wind speed from u and v-components.
    """
    speed = np.sqrt(u * u + v * v)
    return speed

def dewpoint_from_rh(temperature, rh):
    r"""Calculate the ambient dewpoint given air temperature and relative
    humidity.
    """
    if np.any(rh > 1.2):
        print('Relative humidity >120%, ensure proper units.')
    return dewpoint(rh * saturation_vapor_pressure(temperature))

def dewpoint(e):
    r"""Calculate the ambient dewpoint given the vapor pressure.
    """
    val = np.log(e / sat_pressure_0c)
    return 243.5 * val / (17.67 - val)

def saturation_vapor_pressure(temperature):
    r"""Calculate the saturation water vapor (partial) pressure.
    """
    # Converted from original in terms of C to use kelvin. Using raw absolute
    # values of C in a formula plays havoc with units support.
    return sat_pressure_0c * np.exp(17.67 * (temperature - 273.15)
                                    / (temperature - 29.65))

def vapor_pressure(pressure, mixing):
    r"""Calculate water vapor (partial) pressure.
    Given total `pressure` and water vapor `mixing` ratio, calculates the
    partial pressure of water vapor.
    """
    return pressure * mixing / (mpconsts.epsilon + mixing)

#import metpy.calc as mpcalc
#from metpy.units import units

# Data Read
fname = '/Users/leecarlaw/scripts/era5-spc/data/20200701_00-20200701_22.nc'
ds = xr.open_dataset(fname)

pres = ds.level.values
hghts = ds.z.values / G
uwnd = ds.u.values 
vwnd = ds.v.values 
wdirs = wind_direction(uwnd, vwnd)
wspds = wind_speed(uwnd, vwnd)

rh = ds.r.values / 100.
tmpc = ds.t.values
dwpc = dewpoint_from_rh(tmpc, rh)
tmpc = tmpc - ZEROCNK

t = 0
for j in range(tmpc.shape[2]):
    for i in range(tmpc.shape[3]):
        prof = SHARPPY_PROFILE.create_profile(pres=pres,tmpc=tmpc[t,:,j,i],
                          hght=hghts[t,:,j,i],dwpc=dwpc[t,:,j,i],
                          wspd=wspds[t,:,j,i],wdir=wdirs[t,:,j,i])
        #pcl = SHARPPY_PARAMS.parcelx( prof, flag=6 )
        eff_inflow = SHARPPY_PARAMS.effective_inflow_layer(prof)
        print(j,i,'$$$$$$$')
#print(pcl.bplus)
#mlpcl = SHARPPY_PARAMS.parcelx( prof, flag=4 )
