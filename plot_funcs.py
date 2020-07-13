import xarray as xr
import pandas as pd
import numpy as np
import numpy.ma as ma
from plotparms import *

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pylab
from scipy.ndimage import gaussian_filter

from sharptab.constants import *
import sharptab.thermo as thermo
import sharptab.profile as profile
import sharptab.params as params
import sharptab.interp as interp
import sharptab.winds as winds
import sharptab.utils as utils

import metpy.calc as mpcalc

from datetime import datetime

def _make_basemap(bounds):
    m_ = Basemap(projection='stere', llcrnrlon=bounds[0], llcrnrlat=bounds[1],
                 urcrnrlon=bounds[2], urcrnrlat=bounds[3], lat_ts=50,
                 lat_0=50, lon_0=-97., resolution='i')
    m_.drawcoastlines(color='#a09ea0', linewidth=2)
    m_.drawstates(color='#a09ea0', linewidth=2)
    m_.drawcountries(color='#a09ea0', linewidth=2)
    m_.drawcounties(color="#e6e6e6", linewidth=1)
    return m_

def make_plots(filename, bounds):

    # Data read in
    data = xr.open_dataset(filename)
    # We don't need the full dataset. Figure out what data lies in our specified
    # domain and subset it.
    ds = data.sel(latitude=slice(bounds[3]+1, bounds[1]-1),
                  longitude=slice(bounds[0]-1, bounds[2]+1))
    pres = ds.level.values
    hght = ds.z.values / G
    uwnd = ds.u.values
    vwnd = ds.v.values
    wdir = thermo.wind_direction(uwnd, vwnd)
    wspd = thermo.wind_speed(uwnd, vwnd) * MS2KTS

    rh = ds.r.values / 100.
    tmpc = ds.t.values - ZEROCNK
    dwpc = thermo.dewpoint_from_rh(tmpc + ZEROCNK, rh)

    lons = ds.longitude.values
    lats = ds.latitude.values
    lons, lats = np.meshgrid(lons, lats)

    fig_aspect = 1.5
    fig_wid = 15
    fig_hght = fig_wid / fig_aspect

    fig = pylab.figure(figsize=(fig_wid, fig_hght), dpi=100)
    ax_left = 0.01
    ax_bot = 0.0475
    ax_hght = 0.935
    ax_wid = ax_hght / fig_aspect
    ax = pylab.axes([ax_left, ax_bot, 0.975, ax_hght], xticks=[], yticks=[])

    # Create the re-usable basemap plotting object
    m = _make_basemap(bounds)
    x, y = m(lons, lats)
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

    for t in range(tmpc.shape[0]):
        dt = pd.to_datetime(ds.time.values[t])
        valid_date = datetime.strftime(dt, '%Y%m%d/%H00')
        save_date = datetime.strftime(dt, '%Y%m%d_%H')

        # SHARPpy plots and routines. Includes thermodynamics and any wind shear
        # calculations
        prof_data = {'pres':pres, 'tmpc':tmpc[t], 'dwpc':dwpc[t],
                     'hght':hght[t], 'wdir':wdir[t], 'wspd':wspd[t]}
        print("Entering SHARPpy for calculations")
        arrs = sharppy_calcs(**prof_data)

        for parm in sharppy_dict.keys():
            cf_data = None; cf = None; cb = None
            c1_data = None; c1 = None; clab1 = None
            c2_data = None; c2 = None; clab2 = None
            barb = None
            meta = sharppy_dict[parm]

            if meta['cf_data']:
                cf_data = gaussian_filter(arrs[meta['cf_data']], sigma=sigma)
                cf = m.contourf(x, y, cf_data, meta['plot_levs'][0],
                                **plot_kwargs[parm][0])

            if meta['c1_data']:
                if parm not in ['eshr']:    # Careful with nan-possible parms
                    c1_data = gaussian_filter(arrs[meta['c1_data']], sigma=sigma)
                else:
                    c1_data = arrs[meta['c1_data']]
                c1 = m.contour(x, y, c1_data, meta['plot_levs'][1],
                               **plot_kwargs[parm][1])
                clab1 = pylab.clabel(c1, fmt="%d", inline_spacing=8)

            if meta['c2_data']:
                c2_data = gaussian_filter(arrs[meta['c2_data']], sigma=sigma)
                c2 = m.contour(x, y, c2_data, meta['plot_levs'][2],
                               **plot_kwargs[parm][2])
                clab2 = pylab.clabel(c2, fmt="%d", inline_spacing=8)

            # Wind barbs
            if meta['plot_barbs'][0]:
                u = arrs[meta['plot_barbs'][1][0]]
                v = arrs[meta['plot_barbs'][1][1]]
                if meta in ['500mb', '700mb', '300mb', '250mb']:
                    barb_color='k'
                else:
                    barb_color = '#d5ae71'
                barb = m.barbs(x[::skip,::skip],y[::skip,::skip],u[::skip,::skip],
                               v[::skip,::skip], length=6.25, linewidth=0.45,
                               zorder=99, sizes=dict(width=0.15,emptybarb=0.15,
                               spacing=0.1375, height=0.3), flagcolor=None,
                               color=barb_color)

            # Colorbar
            if meta['cf_data']:
                cax = inset_axes(ax, width="25%", height="2.5%", loc="lower left",
                                 bbox_to_anchor=(-0.002,-0.055,1,1),
                                 bbox_transform=ax.transAxes)
                cb = pylab.colorbar(cf, cax=cax, orientation='horizontal',
                                    extend='max')
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_label_position('top')
                cax.xaxis.set_tick_params(pad=0.001)
                cax.tick_params(labelsize=10)

            # Image text and information headers
            img_info = "%s %s" % (valid_date, meta['plot_info'])
            t1 = ax.annotate(img_info, xy=(0.5, -0.031), va='top', ha='center',
                             xycoords='axes fraction', fontsize=11)
            t2 = ax.annotate('ERA5 0.25 x 0.25 deg Global Reanalysis', xy=(0,1),
                             va='bottom', xycoords='axes fraction', fontsize=12,
                             color='b')
            save_name = "%s/%s_%s_%s.png" % (PLOT_DIR, 'MW', parm, save_date)
            #pylab.savefig(save_name, dpi=pylab.gcf().dpi, transparent=True)
            pylab.savefig(save_name, dpi=pylab.gcf().dpi)
            _clean_objects(cf, cb, c1, c2, clab1, clab2, t2, t1, barb)
    pylab.close()

def _clean_objects(cf, cb, c1, c2, clab1, clab2, t2, t1, barb):
    """Remove previous plotting objects"""
    if cf:
        for member in cf.collections:
            member.remove()
    if cb: cb.remove()
    if c1:
        for member in c1.collections:
            member.remove()
    if c2:
        for member in c2.collections:
            member.remove()
    #if c3 is not None:
    #    for member in c3.collections:
    #        member.remove()
    #if c4 is not None:
    #    for member in c4.collections:
    #        member.remove()
    if barb:
        barb[0].remove()
        barb[1].remove()
    if clab1:
        for label in clab1:
            label.remove()
    if clab2:
        for label in clab2:
            label.remove()

    #if clab3 is not None:
    #    for label in clab3:
    #        label.remove()
    t1.remove()
    t2.remove()


def sharppy_calcs(**kwargs):
    tmpc = kwargs.get('tmpc')
    dwpc = kwargs.get('dwpc')
    hght = kwargs.get('hght')
    wdir = kwargs.get('wdir')
    wspd = kwargs.get('wspd')
    pres = kwargs.get('pres')

    # For our jitted sharppy routines, need to be REALLY careful with units or
    # things will break. Probably overkill here, but worth it to be safe!
    tmpc = np.array(tmpc, dtype='float64')
    dwpc = np.array(dwpc, dtype='float64')
    hght = np.array(hght, dtype='float64')
    wdir = np.array(wdir, dtype='float64')
    wspd = np.array(wspd, dtype='float64')
    pres = np.array(pres, dtype='int32')

    mucape = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    mlcape = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    mlcin = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    mulpl = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    ebot = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    etop = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    ebwd_u = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    ebwd_v = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    eshr = np.zeros((tmpc.shape[1], tmpc.shape[2]))

    t1 = datetime.now()
    for j in range(tmpc.shape[1]):
        for i in range(tmpc.shape[2]):
            prof = profile.create_profile(pres=pres, tmpc=tmpc[:,j,i],
                                          hght=hght[:,j,i], dwpc=dwpc[:,j,i],
                                          wspd=wspd[:,j,i], wdir=wdir[:,j,i])

            # Effective inflow and shear calculations
            eff_inflow = params.effective_inflow_layer(prof)
            ebot[j,i] = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
            etop[j,i] = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))

            # This isn't quite right...need to find midpoint between eff Inflow
            # bottom and EL
            ebwd_u[j,i], ebwd_v[j,i] = winds.wind_shear(prof, pbot=eff_inflow[0],
                                                        ptop=500)
            eshr[j,i] = utils.mag(ebwd_u[j,i], ebwd_v[j,i])

            # Parcel buoyancy calculations
            mupcl = params.parcelx(prof, flag=3)
            mlpcl = params.parcelx(prof, flag=4)
            mucape[j,i] = mupcl.bplus
            mulpl[j,i] = mupcl.lplhght
            mlcape[j,i] = mlpcl.bplus
            mlcin[j,i] = mlpcl.bminus

    eshr = ma.where(eshr < 25., np.nan, eshr)
    ebwd_u = ma.where(eshr < 25., np.nan, ebwd_u)
    ebwd_v = ma.where(eshr < 25., np.nan, ebwd_v)

    t2 = datetime.now()
    print("Calculations took: ", str((t2-t1).seconds), " seconds")
    ret = {'mucape' : mucape,
           'mlcape' : mlcape,
           'mlcin' : mlcin,
           'mulpl' : mulpl,
           'ebot' : ebot,
           'etop': etop,
           'ebwd_u': ebwd_u,
           'ebwd_v': ebwd_v,
           'eshr': eshr
    }
    return ret
