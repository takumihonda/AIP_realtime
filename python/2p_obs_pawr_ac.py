from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import sys
import os

from tools_AIP import prep_proj_multi_cartopy, get_cfeature, setup_grids_cartopy, read_obs_nc, read_mask_full
import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

quick = True
#quick = False


def main( INFO={}, vtime=datetime( 2019, 8, 24, 15, 10, 0), 
          height=3000, exp_l=[],
          texp_l=[], odz=300 ):

    mask, mlon2d, mlat2d = read_mask_full()
    mzmax = mask.shape[0]
    mz1d = np.arange( 0, 500*(mzmax+1), 500 )
    mzlev = np.argmin( np.abs( mz1d - height ) )
   
    mask2d = mask[mzlev,:,:]
    mask2d = np.where( mask2d < 0.5, np.nan, 1 ) 


    odat_l = []
    olon_l = []
    olat_l = []

    for exp_ in exp_l:
        ofn = '{0:}/{1:}/pawr_grads/pawr_obs_{2:}.nc'.format( INFO["TOP"], exp_,
                     vtime.strftime('%Y%m%d-%H%M%S') )
        odat, olon, olat, olev, rlon, rlat = read_obs_nc( fn=ofn, dz=odz )

        odat_l.append( odat )
        olon_l.append( olon )
        olat_l.append( olat )

    lons = 138.90161
    lone = 140.31638
    lats = 35.285637
    late = 36.43222

    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )

    res = '10m'
    if quick:
       res = '50m'

    land = get_cfeature( typ='land', res=res )
    coast = get_cfeature( typ='coastline', res=res )

    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    # radar location
    lon_r = 139.609
    lat_r = 35.861

    fig = plt.figure( figsize=(10, 5) )
    fig.subplots_adjust( left=0.06, bottom=0.03, right=0.95, top=0.95,
                         wspace=0.2, hspace=0.01)

    ax_l = prep_proj_multi_cartopy( fig, xfig=2, yfig=1, proj='merc', 
                         latitude_true_scale=lat_r )
 
    levels= np.array( [ 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 ] )
    cmap = mcolors.ListedColormap(['cyan', 'b', 'dodgerblue',
                                   'lime','yellow',
                                   'orange', 'red', 'firebrick', 'magenta',
                                    'purple'])
    cmap.set_over('gray', alpha=1.0)
    cmap.set_under('w', alpha=0.0)
    cmap.set_bad( color='gray', alpha=0.5 )

    norm = BoundaryNorm( levels, ncolors=cmap.N, clip=False)


    for i, ax in enumerate( ax_l ):

        ax.set_extent([ lons, lone, lats, late ] )
        ax.add_feature( coast, zorder=0 )
 
        setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                             fs=9, lw=0.0 )
 
        SHADE = ax.scatter( olon_l[i], olat_l[i], s=5.0, c=odat_l[i], 
                    norm=norm, 
                    cmap=cmap,
                    transform=data_crs ) 
 
        SHADE2 = ax.contourf( mlon2d, mlat2d, mask2d,
                             transform=data_crs, 
                             levels=[0.5, 1.5],
                             colors=['k', 'gray', 'b'],
                             alpha=0.3,
                             extend='both', )

        ax.plot( rlon, rlat, marker='o', color='k', markersize=10,
                 transform=data_crs )
 
        ax.text( 0.5, 1.01, texp_l[i],
                va='bottom', 
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=12, )
 
        bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
                 'edgecolor':'w' }

        ax.text( 0.99, 0.01, r'Z={0:.1f}$\pm${1:.1f} km'.format( height*0.001, odz*0.001 ), 
                va='bottom', 
                ha='right',
                bbox=bbox,
                transform=ax.transAxes,
                color='k', fontsize=10, )


        if i == 1:

           cvtime = vtime.strftime( '%H:%M:%S' )
           ax.text( 1.0, 1.01, 'Valid at {0:}'.format( cvtime ),
                   va='bottom', 
                   ha='right',
                   transform=ax.transAxes,
                   color='k', fontsize=10, )

           pos = ax.get_position()
           cb_width = 0.006
           cb_height = pos.height*0.9
           ax_cb = fig.add_axes( [ pos.x1+0.003, pos.y0, 
                                   cb_width, cb_height] )
           cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                              ticks=levels[::1], extend='both' )
    
           ax.text( 1.01, 0.95, '(dBZ)', 
                   va='top', 
                   ha='left',
                   transform=ax.transAxes,
                   color='k', fontsize=10, )

           


    tit = 'MP-PAWR observations'
    fig.suptitle( tit, fontsize=14 )

    ofig = "2p_obs_ac_h{0:}km_v{1:}.png".format( height*0.001, vtime.strftime('%Y%m%d%H%M%S') )

    print( ofig )
    if not quick:
       opath = "png/attenuation"
       os.makedirs( opath, exist_ok=True )
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
    else:
       plt.show()




INFO = { "TOP": "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test" }

vtime = datetime( 2019, 8, 24, 15, 30, 0 )
#vtime = datetime( 2019, 8, 24, 15, 40, 0 )
#vtime = datetime( 2019, 8, 24, 15, 50, 0 )
#vtime = datetime( 2019, 8, 24, 16,  0, 0 )

exp_l = [ "d4", "d4_attenuation_corrected"]
texp_l = [ "CTRL", "TEST" ]

height = 3000.0
odz = 300

stime = datetime( 2019, 8, 24, 15, 30, 0 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

etime = stime

time = stime
while time <= etime:
   main( INFO=INFO, vtime=time, exp_l=exp_l,
         texp_l=texp_l, height=height, odz=odz )

   time += timedelta( minutes=10 )
