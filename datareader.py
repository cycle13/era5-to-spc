import requests
import subprocess
import os

from scipy import spatial 
import numpy as np

from airports import *
from mapinfo import domains

from metpy.units import units
import metpy.calc as mpcalc

from datetime import datetime, timedelta

import cdsapi
from era5_config import *
from collections import defaultdict

def download_radar(time, data_path):
    #time += timedelta(hours=1)
    date_str = str(time.year) + str(time.month).zfill(2) + str(time.day).zfill(2)
    hour_str = str(time.hour).zfill(2)
    fname = "SeamlessHSR_00.00_%s-%s0000" % (date_str, hour_str)
    url = "%s/%s/%s/%s/mrms/ncep/SeamlessHSR/%s.grib2.gz"% (_mrms_base_url, str(time.year), 
                                                            str(time.month).zfill(2), 
                                                            str(time.day).zfill(2),
                                                            fname)
    full_name = "%s/%s.nc" % (data_path, fname)
    print(full_name)

    if not os.path.exists(full_name):
        
        # Test for file existence
        r = requests.head(url)
        try:
            resp_length = int(r.headers['content-length'])
        except:
            print("Can't locate MRMS data. Trying earlier time...")
            try:
                time -= timedelta(minutes=2)
                date_str = str(time.year) + str(time.month).zfill(2) + str(time.day).zfill(2)
                hour_str = str(time.hour).zfill(2) + str(time.minute).zfill(2)
                fname2 = "SeamlessHSR_00.00_%s-%s00" % (date_str, hour_str)
                url2 = "%s/%s/%s/%s/mrms/ncep/SeamlessHSR/%s.grib2.gz"% (_mrms_base_url, str(time.year),
                                                                         str(time.month).zfill(2),
                                                                         str(time.day).zfill(2),
                                                                         fname2)
                r = requests.head(url2)
                resp_length = int(r.headers['content-length'])
                url = url2
                fname = fname2
            except:
                #raise ValueError("Can't locate MRMS data")
                pass

        shell = 'wget %s -O %s' % (url, data_path + '/' + fname + '.grib2.gz')

        try:
            subprocess.call(shell, shell=True)
            shell = 'gunzip %s' % (data_path + '/' + fname + '.grib2.gz')
            subprocess.call(shell, shell=True)
            shell = 'wgrib2 %s -netcdf %s' % (data_path + '/' + fname + '.grib2',
                                              full_name)
            subprocess.call(shell, shell=True)
        except:
            full_name = None
            pass
            #raise ValueError("Error Downloading MRMS data from: %s" % (url))

    else:
        print("MRMS data already exists locally.")        
    
    return full_name

def nearest_idx(points, lon, lat, tree=None):
    """
    Search for the nearest grid point using a KDTree
    points = [[LON1, LAT1], [LON2, LAT2]]
    lat: 2-d array of latitudes
    lon: 2-d array of longitudes
    """
    if tree is None:
        lonlat = np.column_stack((lon.ravel(), lat.ravel()))
        tree = spatial.cKDTree(lonlat)
        dist, idx = tree.query(points, k=1)
        ind = np.column_stack(np.unravel_index(idx, lon.shape))
    return [(j,i) for j,i in ind]

def _write_to_sounding(data_cube, idx_locs, ids, date, fmt=None):
    """
    Writes data to sounding files. Added BUFKIT-readable output capabilities. 
    """
    if fmt is None:
        fmt = 'sharppy'

    knt = 0
    for idx in idx_locs:
        t_out = data_cube[0,:,idx[0],idx[1]] - 273.15
        td_out = data_cube[1,:,idx[0],idx[1]]
        u = data_cube[2,:,idx[0],idx[1]] * units('m/s')
        v = data_cube[3,:,idx[0],idx[1]] * units('m/s')
        wdir_out = mpcalc.wind_direction(u, v).magnitude
        wspd_out = mpcalc.wind_speed(u, v).magnitude
        hgt_out = data_cube[4,:,idx[0],idx[1]]

        out_time = "%s%s%s/%s00"%(date[2:4], date[4:6], date[6:8], date[8:10])
        if fmt == 'sharppy':
            out_file = "%s/%s.%s" % (SOUNDING_DIR, date, ids[knt])
            f = open(out_file, 'w')

            f.write("%TITLE%\n")
            f.write(" %s   %s" % (ids[knt], out_time))
            f.write("\n\n")
            f.write('   LEVEL     HGHT     TEMP     DWPT     WDIR     WSPD\n')
            f.write('------------------------------------------------------\n')
            f.write('%RAW%\n')
            
            # This is a weird one. At some point in the past, the GRIB files
            # were ordered differently. Need to check for monotonic increasing
            # or decreasing heights and adjust pressures accordingly
            if strictly_increasing(list(hgt_out)):
                pres_incr = 1
                start_, end_, inc_ = 0, t_out.shape[0], 1
            else:
                pres_incr = -1
                start_, end_, inc_ = t_out.shape[0]-1, -1, -1

            print(levs)
            print(hgt_out)
            print(start_,end_,inc_, pres_incr)
            for row in range(start_, end_, inc_):
                if hgt_out[row] > 0:
                    out_line = "%s,%s,%s,%s,%s,%s" % (levs[::pres_incr][row],
                                                  hgt_out[row],
                                                  t_out[row],
                                                  td_out[row],
                                                  wdir_out[row],
                                                  wspd_out[row])
                    f.write(out_line +'\n')
            f.write('%END%')        
        knt += 1

def get_soundings(data, ids, lats, lons, date):
    """
    ids = [ORD, MDW, etc]
    Gathers the data for writing to sounding files
    """
    date = date.strftime("%Y%m%d%H")
    points = []
    elevs = []
    for id_ in ids:
        entry = metadata[id_]
        points.append([entry[1], entry[0]])

    idx_locs = nearest_idx(points, lons, lats)

    t = data.select(name='Temperature', level=levs)
    rh = data.select(name='Relative humidity', level=levs)
    u = data.select(name='U component of wind', level=levs)
    v = data.select(name='V component of wind', level=levs)
    hgt = data.select(name='Geopotential Height', level=levs)

    n_levs = len(levs)
    data_cube = np.zeros((5,n_levs,t[0].values.shape[0],t[0].values.shape[1]))
    #data_cube[0,0,:,:] = t_2m
    #data_cube[1,0,:,:] = td_2m
    for k in range(0, n_levs):
        t_lev = t[k].values * units.kelvin
        rh_lev = rh[k].values/100.
        rh_lev = rh_lev * units('dimensionless')
        td_lev = mpcalc.dewpoint_rh(t_lev, rh_lev)

        u_lev = u[k].values*units('m/s').to('knots')
        v_lev = v[k].values*units('m/s').to('knots')
        
        data_cube[0,k,:,:] = t_lev
        data_cube[1,k,:,:] = td_lev
        data_cube[2,k,:,:] = u_lev 
        data_cube[3,k,:,:] = v_lev 
        data_cube[4,k,:,:] = hgt[k].values

    _write_to_sounding(data_cube, idx_locs, ids, date, fmt='sharppy')
    #_write_to_sounding(data_cube, idx_locs, ids, date, fmt='bufkit')
