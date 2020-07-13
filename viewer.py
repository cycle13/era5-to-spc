"""
viewer.py
This is the main script controlling the download, plotting, and sounding output
functionalities.

The archived RAP data is accessed via the NCI THREDDS catalog server at 
https://www.ncei.noaa.gov/thredds/catalog/rap130anl/catalog.html. This data 
archive currently goes back to May of 2012. Unforuntately, index files are not
available on this server, so we have to resort to downloading the entire file(s)
on the local disk. Thankfully, the files are relatively small at about 15 MB or
less. A WGET executable must be installed on the system and the pygrib module
must be available for python to read the GRIB2 data. 

Author:         Lee Carlaw
Completed:      September 2019
Modified:
"""

import numpy as np
from collections import defaultdict
import re

from mapinfo import domains
from datareader import download_era5, get_soundings, download_radar
from plot import plot
from plotparms import *

import argparse
import sys, os
from datetime import datetime, timedelta
#from distutils.spawn import find_executable
#if find_executable('wget') is None:
#    raise ValueError("wget executable not found. Program exiting.")
#    sys.exit()
#try:
#    import pygrib
#except:
#    raise ValueError("Python can't find pygrib. Program exiting.")
#    sys.exit()
import xarray as xr

base_data_path = './data/'

def parse_time(time_str):
    plot_time = datetime.strptime(time_str, '%Y-%m-%d/%H')
    return plot_time

def plotter(time, domain, data_path=None, ids=None):
    if data_path is None:
        data_path = base_data_path + '/' + time[:-3]
    if not os.path.exists(data_path):
        os.makedirs(data_path)

    plot_time = parse_time(time)
    #run_time = plot_time - timedelta(hours=1)

    _file = download_era5(plot_time, data_path, domain)
    _radar_file = download_radar(plot_time, data_path)

    # Determine what data needs to be plotted and passed
    #data = pygrib.open(_file)
    #lats, lons = data[1].latlons()
   

    data = xr.open_dataset(_file)
    lats = data.latitude.values
    lons = data.longitude.values
    lons, lats = np.meshgrid(lons, lats)

    if domain is not None:
        plot_bounds = domains[domain]
        plot(data, lons, lats, plot_bounds, domain, plot_time, _radar_file)

    if ids is not None:
        try:
            ids = ids.strip().split(',')
        except:
            raise ValueError("Improperly formatted sounding locations: ORD,MDW,DPA")      
        print(ids)
        get_soundings(data, ids, lats, lons, plot_time)

    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-t', '--time', dest='time', help="YYYY-MM-DD/HH The desired valid time of the archived RAP data to view.")
    ap.add_argument('-d', '--domain', dest='domain', help="[MW|SGP|CGP|NGP|GL|SE|MA|NE|GL|NW|GBSN|SW|CONUS] Plotting domain string. Set to None for no plot.")
    ap.add_argument('-s', '--sounding', dest='ids', help="ORD,MDW,DPA If airport identifiers are specified, outputs SHARPPy-readable sounding files. See https://rucsoundings.noaa.gov/airports.html for airport codes")
    #ap.add_argument('-f', '--forecast', dest='forecast', help="Forecast. Specify run as YYYY-MM-DD/HH")
    ap.add_argument('-p', '--data_path', dest='data_path', help="Path location for downloaded files. If none specified, creates a subdirection YYYYMMDD in the current working directory")
    args = ap.parse_args()
    np.seterr(all='ignore')

    plotter(args.time,
            args.domain,
            data_path = args.data_path,
            ids = args.ids
            )
        
if __name__ == "__main__":
    main()
