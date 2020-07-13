import numpy as np
import xarray as xr
from sharptab.constants import G, ZEROCNK 

import sharptab.thermo as thermo
import sharptab.params as params
import sharptab.profile as profile
import sharptab.interp as interp

# Data Read
fname = '/Users/leecarlaw/scripts/era5-spc/data/20200701_00-20200701_22.nc'
ds = xr.open_dataset(fname)

pres = ds.level.values
hghts = ds.z.values / G
uwnd = ds.u.values
vwnd = ds.v.values
wdirs = thermo.wind_direction(uwnd, vwnd)
wspds = thermo.wind_speed(uwnd, vwnd)

rh = ds.r.values / 100.
tmpc = ds.t.values - ZEROCNK
dwpc = thermo.dewpoint_from_rh(tmpc + ZEROCNK, rh)

# For our jitted sharppy routines, need to be carefuly about units. The float64
# type is standard (type 'double'), but specify here just to be safe. 
tmpc = np.array(tmpc, dtype='float64')
dwpc = np.array(dwpc, dtype='float64')
hghts = np.array(hghts, dtype='float64')
wdirs = np.array(wdirs, dtype='float64')
wspds = np.array(wspds, dtype='float64')
pres = np.array(pres, dtype='int32')

t,j,i = 0,0,70
prof = profile.create_profile(pres=pres,tmpc=tmpc[t,:,j,i],
                              hght=hghts[t,:,j,i],dwpc=dwpc[t,:,j,i],
                              wspd=wspds[t,:,j,i],wdir=wdirs[t,:,j,i])
#print(prof.tmpc, prof.dwpc)
#mupcl = params.parcelx( prof, flag=3 )
#effpcl = params.effective_inflow_layer( prof, mupcl)
#mlpcl = params.parcelx( prof, flag=6 )
eff_inflow = params.effective_inflow_layer(prof)
ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
etop_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))

print(ebot_hght, etop_hght)

#print(mupcl.bplus)
