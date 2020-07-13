import argparse
from plotparms import *
from utils.mapinfo import domains
from plot_funcs import make_plots

def plotter(filename, domain=None, ids=None):
    '''
    '''
    if not domain: domain='US'
    plot_bounds = domains[domain]
    make_plots(filename, plot_bounds)
    

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(dest='filename', help="Full path to netCDF4 file")
    ap.add_argument('-d', '--domain', dest='domain', help="[MW|SGP|CGP|NGP|GL|SE|MA|NE|GL|NW|GBSN|SW|CONUS] Plotting domain string. Set to None for no plot.")
    ap.add_argument('-s', '--ids', dest='ids', help="Plot soundings: ORD,VPZ,DFW")
    args = ap.parse_args()

    plotter(args.filename, args.domain, ids=args.ids)

if __name__ == "__main__": main()
