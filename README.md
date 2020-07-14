# era5-to-spc
Scripts in this repo attempt to re-create SPC's [Mesoanalysis Archive](https://www.spc.noaa.gov/exper/mesoanalysis/) using the ECMWF's [ERA5 Dataset](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5).

The main overhaul here was to massively accelerate the CPU-intensive thermodynamic calculations performed by SHARPpy (mainly due to constant lifting of parcels) using the [Python Numba](http://numba.pydata.org/) module. This is a non-trivial task, as several standard Python modules and code are not supported by Numba. In addition, the nature of the "Just-in-time" compilation requires explicit `type` declarations within Python `Classes`. As a result, the original SHARPpy code had to be parsed out line-by-line to allow it to work with Numba and the jit decorator, and some flexibility has certainly been lost here. The biggest issues were the lack of `**kwarg` support and numpy masked arrays. In this current iteration of code, it's assumed that the input meteorological arrays are full without any missing/masked data. 

## Code Execution Time Improvements
The initial "Just-in-time compiling" of the associated SHARPpy code (creating the Profile and Parcel classes) takes about 30 seconds, and this will likely be consistent across many cases. Each subsequent time slice [for a full CONUS domain] takes about 10-20 seconds, and about 4 seconds for a smaller "Midwest" domain test. Compare this with the straight-out-of-the-box SHARPpy code (running in serial) where each time slice takes upwards of 25 to 30 minutes!

* Simple `tqdm` output for a time slice (~10 seconds):
![](https://github.com/lcarlaw/era5-to-spc/raw/master/readme_images/numba.png?raw=true)

* Same for straight serial out of the box run (~30 minutes):
![](https://github.com/lcarlaw/era5-to-spc/raw/master/readme_images/serial.png?raw=true)

## Basic Setup Notes
You must first create an account on the [CDS Site](https://cds.climate.copernicus.eu/#!/home). Once you register, you'll be directed to copy off an API key and a url string that you'll save into a file: `$HOME/.cdsapirc`.

### Creating the base environment
The setup here proceeds using Anaconda, as well as assuming a completely vanilla Python3 install.  I've also edited my `~/.condarc` file to add conda-forge to the default channels. The order is important here for versioning control as newer versions of matplotlib and basemap cause issues with some of the map features.

```
conda create --name era5 python=3.7
conda activate era5

conda install matplotlib=2.2.4
conda install basemap=1.2.0 basemap-data-hires=1.2.0
conda install xarray netcdf4
conda install cdsapi
conda install numba
conda install metpy
```

Export the conda `.yml` file with:

```
conda env export --from-history -f environment.yml
```

## Useage
General usage information.

### Downloading ERA5 Data
```
python get_era5.py -s START_TIME -e END_TIME [-d DOMAIN]
```
`START_TIME`: Initial (earliest) slice of data request. Form is `YYYY-MM-DD/HH`
`END_TIME`: Latest (last) slice of data request. Form is `YYYY-MM-DD/HH`
`DOMAIN`: Optional argument to download a specific domain subset. If left off, this will default to the CONUS region. The acceptable values are located in the `./utils/mapinfo.py` file. There is a helper script called `./utils/view_map.py` to visualize the domains that will be downloaded and plotted.

NetCDF4 files will be saved into a folder called `data` in your present working directory.

Note that depending on the size of your file request (length of time, domain size), download/queuing times can become quite long, potentially on the order of days. You can see the live [queue here](https://cds.climate.copernicus.eu/live/queue), and also actively monitor [your requests here](https://cds.climate.copernicus.eu/cdsapp#!/yourrequests).

### Image Plotting
Plotting is controlled at the top level by `main.py` which accepts several command-line arguments and passes them to the various plotting functions. The general use is:
```
python main.py FILENAME [-d DOMAIN]
```
`FILENAME`: Required argument specifying the full path to the ERA5 netCDF file.
`DOMAIN`: Optional argument specifying the domain to be plotted. These are defined in the `./utils/mapinfo.py` file and can be edited.

The `plotparms.py` file controls most of the plotting configuration, from linestyles and colors, to wind barbs and text labels. It's a goofy file since it has to control so much, and there are undoubtedly much better ways to handle this...
