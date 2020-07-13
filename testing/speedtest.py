import numpy as np
import xarray as xr
from sharptab.constants import G, ZEROCNK 

import sharptab.thermo as thermo
import sharptab.params as params
import sharptab.profile as profile

#import sharppy.sharptab.profile as SHARPPY_PROFILE
#import sharppy.sharptab.params as SHARPPY_PARAMS
#import sharppy.sharptab.thermo as thermo


from datetime import datetime

from tqdm import tqdm

# Data Read
t1 = datetime.now()
fname = '/Users/leecarlaw/scripts/era5-spc/data/20200701_00-20200701_02.nc'
ds = xr.open_dataset(fname)

pres = ds.level.values
hghts = ds.z.values / G
uwnd = ds.u.values
vwnd = ds.v.values
#wdirs = uwnd
#wspds = vwnd
wdirs = thermo.wind_direction(uwnd, vwnd)
wspds = thermo.wind_speed(uwnd, vwnd)

rh = ds.r.values / 100.
tmpc = ds.t.values - ZEROCNK
dwpc = thermo.dewpoint_from_rh(tmpc + ZEROCNK, rh)

t2 = datetime.now()
print("Data read-in completed in: %s seconds" % (t2-t1).seconds)

# For our jitted sharppy routines, need to be carefuly about units. The float64
# type is standard (type 'double'), but specify here just to be safe. 
tmpc = np.array(tmpc, dtype='float64')
dwpc = np.array(dwpc, dtype='float64')
hghts = np.array(hghts, dtype='float64')
wdirs = np.array(wdirs, dtype='float64')
wspds = np.array(wspds, dtype='float64')
pres = np.array(pres, dtype='int32')

npoints = hghts.shape[2] * hghts.shape[3]   
#for t in range(hghts.shape[0]):
for t in range(0,1):
    t3 = datetime.now()
    #for j in range(0,1):
    #    for i in range(0,1):
    for j in range(hghts.shape[2]):
        for i in range(hghts.shape[3]):
            for ii in tqdm(range(npoints)):
                prof = profile.create_profile(pres=pres,tmpc=tmpc[t,:,j,i],
                                      hght=hghts[t,:,j,i],dwpc=dwpc[t,:,j,i],
                                     wspd=wspds[t,:,j,i],wdir=wdirs[t,:,j,i])
                t4 = datetime.now()
                #mupcl = params.parcelx( prof, flag=3 )
                #mlpcl = params.parcelx( prof, flag=4 )
                eff_inflow = params.effective_inflow_layer(prof)
                mupcl = params.parcelx( prof, flag=3 )
                mlpcl = params.parcelx( prof, flag=4 )


    t5 = datetime.now()
    print("Time slice completed in: %s seconds" % (t5-t3).seconds)
'''
for t in range(hghts.shape[0]):
    t3 = datetime.now()
    for j in range(hghts.shape[2]):
        for i in range(hghts.shape[3]):

            for ii in tqdm(range(npoints)):
                prof = SHARPPY_PROFILE.create_profile(pres=pres,tmpc=tmpc[t,:,j,i],
                                          hght=hghts[t,:,j,i],dwpc=dwpc[t,:,j,i],
                                          wspd=wspds[t,:,j,i],wdir=wdirs[t,:,j,i])
                t4 = datetime.now()
                mupcl = SHARPPY_PARAMS.parcelx( prof, flag=3 )
                #print(mupcl.bplus)
    t5 = datetime.now()
    print("Time slice completed in: %s seconds" % (t5-t3).seconds)
''' 







