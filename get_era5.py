import numpy as np
from collections import defaultdict
import argparse
from datetime import datetime, timedelta
import os
import cdsapi

from utils.era5_config import *
from utils.mapinfo import domains

pwd = os.environ['PWD']
def download_era5(start, end, data_path, domain=None):
    c = cdsapi.Client()

    if domain is not None:
        plot_bounds = domains[domain]
    else:
        plot_bounds = domains['US']

    times = defaultdict(list)
    while start <= end:
        time_string = datetime.strftime(start, '%Y%m%d_%H')
        times['time-str'].append(time_string)
        times['years'].append(str(start.year))
        times['months'].append(str(start.month).zfill(2))
        times['days'].append(str(start.day).zfill(2))
        times['hours'].append(str(start.hour).zfill(2) + ':00')     
        start += timedelta(hours=1)

    years = sorted(set(times['years']))
    months = sorted(set(times['months']))
    days = sorted(set(times['days']))
    hours = sorted(set(times['hours']))
    
    outfile = '%s-%s.nc' % (times['time-str'][0], times['time-str'][-1])
    if not os.path.exists(data_path + '/' + outfile):
        c.retrieve(
            'reanalysis-era5-pressure-levels', 
            {
                'product_type': 'reanalysis', 
                'variable': DEFAULT_DATA[pressure_vars],
                'pressure_level': DEFAULT_DATA[pressure_levs],
                'year': years,
                'month': months,
                'day': days,
                'time': hours,
                'format': 'netcdf',
                'area' : [plot_bounds[-1]+1,plot_bounds[0]-1,plot_bounds[1]-1,plot_bounds[2]+1],
                'grid': [0.25, 0.25],
            }, 
            '%s/%s' % (data_path, outfile))
    '''
    # Surface and other single-level variables
    if not os.path.exists(data_path + '/sfc_' + outfile):
        c.retrieve(
            'reanalysis-era5-single-levels', 
            {
                'product_type': 'reanalysis', 
                'variable': DEFAULT_DATA[surface_vars],
                'year': times['years'],
                'month': times['months'],
                'day': times['days'],
                'time': times['hours'],
                'format': 'netcdf',
                'area' : [plot_bounds[-1]+1,plot_bounds[0]-1,plot_bounds[1]-1,plot_bounds[2]+1],
                'grid': [0.25, 0.25],
            }, 
            '%s/sfc_%s' % (data_path, outfile))
    '''
    return

def parse_time(time_str):
    plot_time = datetime.strptime(time_str, '%Y-%m-%d/%H')
    return plot_time

def download(start_time, end_time, domain):
    data_path = pwd + '/data/'
    start = datetime.strptime(start_time, '%Y-%m-%d/%H')
    end = datetime.strptime(end_time, '%Y-%m-%d/%H')
    _file = download_era5(start, end, data_path, domain=None)
    #_radar_file = download_radar(plot_time, data_path)
    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-s', '--start-time', dest='start_time', help="YYYY-MM-DD/HH The desired valid time of the archived RAP data to view.")
    ap.add_argument('-e', '--end-time', dest='end_time', help="YYYY-MM-DD/HH The desired valid time of the archived RAP data to view.")
    ap.add_argument('-d', '--domain', dest='domain', help="[MW|SGP|CGP|NGP|GL|SE|MA|NE|GL|NW|GBSN|SW|CONUS] Plotting domain string. Set to None for no plot.")
    args = ap.parse_args()
    np.seterr(all='ignore')

    download(args.start_time,
             args.end_time,
             args.domain,
            )
    
if __name__ == "__main__":
    main()
