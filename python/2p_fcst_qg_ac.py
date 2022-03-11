from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import sys
import os

from tools_AIP import prep_proj_multi_cartopy, get_cfeature, setup_grids_cartopy, read_fcst_nc
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

def main( INFO={}, vtime=datetime( 2019, 8, 24, 15, 10, 0), height=3000,
          fstime=datetime( 2019, 8, 24, 15, 10, 0 ), exp_l=[],
          texp_l=[], nvar='qg' ):


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

        ffn = '{0:}/{1:}/dafcst_nc/{2:}.nc'.format( INFO["TOP"], exp_l[i],
                     fstime.strftime('%Y%m%d-%H%M%S') )
        print( ffn )
        fdat2d, flon1d, flat1d = read_fcst_nc( fn=ffn, fstime=fstime, vtime=vtime )
     
        flon2d, flat2d = np.meshgrid( flon1d, flat1d )
    
        lons = np.min( flon2d )
        lone = np.max( flon2d )
    
        lats = np.min( flat2d )
        late = np.max( flat2d )


        ax.set_extent([ lons, lone, lats, late ] )
        ax.add_feature( coast, zorder=0 )
 
        setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                             fs=9, lw=0.0 )
 
#        if i == 0:
        SHADE = ax.contourf( flon2d, flat2d, fdat2d[:,:], 
                     norm=norm,
                     levels=levels,
                     cmap=cmap,
                     extend='both',
                     transform=data_crs )
#        elif i == 1:
#           SHADE = ax.scatter( olon, olat, s=5.0, c=odat, 
#                       norm=norm, 
#                       cmap=cmap,
#                       #extend='both',
#                       #vmin=np.min( levels ), 
#                       #vmax=np.max( levels ), 
#                       transform=data_crs ) 
#           print( np.max( odat ) )

        ax.plot( lon_r, lat_r, marker='o', color='k', markersize=10,
                 transform=data_crs )

        ax.text( 0.5, 1.01, texp_l[i],
                va='bottom', 
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=12, )

        bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
                 'edgecolor':'w' }
        ax.text( 0.99, 0.01, 'Z={0:.1f} km'.format( height*0.001 ), 
                va='bottom', 
                ha='right',
                bbox=bbox,
                transform=ax.transAxes,
                color='k', fontsize=10, )

           

        if i == 1:

           cfstime = fstime.strftime( '%H:%M:%S' )
           cvtime = vtime.strftime( '%H:%M:%S' )
           ax.text( 1.0, 1.01, 'Init at {0:}\nValid at {1:}'.format( cfstime, cvtime ),
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




    if fstime == vtime:
       tit_ = 'analysis'
    else:
       tit_ = 'forecast'

    tit = 'SCALE-LETKF {0:}'.format( tit_ )
    fig.suptitle( tit, fontsize=14 )


    ofig = "2p_fcst_ac_h{0:}km_s{1:}_v{2:}_{3:}.png".format( height*0.001, 
                    fstime.strftime('%Y%m%d%H%M%S'),
                    vtime.strftime('%Y%m%d%H%M%S'),
                    nvar,
                     )

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

exp_l = [ "d4", "d4_attenuation_corrected"]
texp_l = [ "CTRL", "TEST" ]

vtime = datetime( 2019, 8, 24, 15, 30, 0 )
#vtime = datetime( 2019, 8, 24, 15, 40, 0 )
#vtime = datetime( 2019, 8, 24, 15, 50, 0 )
#vtime = datetime( 2019, 8, 24, 16, 0, 0 )
fstime = datetime( 2019, 8, 24, 15, 30, 0 )

height = 3000
nvar = 'qg'


main( INFO=INFO, vtime=vtime, fstime=fstime, exp_l=exp_l, 
      texp_l=texp_l, height=height, nvar=nvar )
