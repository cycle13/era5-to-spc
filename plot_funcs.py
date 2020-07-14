import xarray as xr
import pandas as pd
import numpy as np
from plotparms import *
from sharptab.constants import G, ZEROCNK, MS2KTS
import sharptab.sharppy_funcs as sharppy_funcs
import sharptab.thermo as thermo

import pylab
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import gaussian_filter

import metpy.calc as mpcalc
from datetime import datetime

"""
    plot_funcs.py

    Primary function controlling the actual plot creation. Plot parameters are
    specified in the plotparms.py file within this top-level directory.

"""

def make_basemap(bounds):
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
    # domain and subset it. Give a litte buffer for the map projection.
    ds = data.sel(latitude=slice(bounds[3]+2, bounds[1]-2),
                  longitude=slice(bounds[0]-2, bounds[2]+2))
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
    m = make_basemap(bounds)
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
        arrs = sharppy_funcs.calcs(**prof_data)
        for parm in PLOT_DICT.keys():
            cf_data = None; cf = None; cb = None
            c1_data = None; c1 = None; clab1 = None
            c2_data = None; c2 = None; clab2 = None
            barb = None
            meta = PLOT_DICT[parm]

            if meta['cf_data']:
                cf_data = gaussian_filter(arrs[meta['cf_data']], sigma=sigma)
                cf = m.contourf(x, y, cf_data, meta['plot_levs'][0],
                                **PLOT_KWARGS[parm][0])

            if meta['c1_data']:
                # Carefule with nan-possible parms such as effective inflow and
                # other shear calculations.
                if parm not in ['eshr', '']:
                    c1_data = gaussian_filter(arrs[meta['c1_data']], sigma=sigma)
                else:
                    c1_data = arrs[meta['c1_data']]

                c1 = m.contour(x, y, c1_data, meta['plot_levs'][1],
                               **PLOT_KWARGS[parm][1])
                clab1 = pylab.clabel(c1, fmt="%d", inline_spacing=8)

            if meta['c2_data']:
                c2_data = gaussian_filter(arrs[meta['c2_data']], sigma=sigma)
                c2 = m.contour(x, y, c2_data, meta['plot_levs'][2],
                               **PLOT_KWARGS[parm][2])
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
                cax = inset_axes(ax, width="18%", height="2.5%", loc="lower left",
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
            t2 = ax.annotate('ERA5 0.25 x 0.25 deg Global Reanalysis-SPC Meld Project',
                             xy=(0,1), va='bottom', xycoords='axes fraction',
                             fontsize=12, color='b')
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
