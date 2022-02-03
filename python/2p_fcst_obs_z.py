from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import sys

from tools_AIP import prep_proj_multi_cartopy, get_cfeature, setup_grids_cartopy
import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

quick = True

def read_obs_nc( fn='', height=3000.0, dz=100.0, otyp=4001, DBZ=False ):
    nc = Dataset( fn, "r", format="NETCDF4" )
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    lev = nc.variables['lev'][:]
    dat = nc.variables['dat'][:]
    elm = nc.variables['elm'][:]
    nc.close()

    print( dat.shape )
    #dat = dat[ ( elm == 4003 ) * ( lev == height ) ]
    dat = dat[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]
    lon = lon[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]
    lat = lat[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]

    lev = lev[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]
    #print( lev )

    if not DBZ:
       dat = 10 * np.log10( np.where( dat>0.0001, dat, 0.0001 ) )

    return( dat, lon, lat )

def read_nc( fn='', fstime=datetime( 2019, 8, 24, 15, 10, 0 ),
             vtime=datetime( 2019, 8, 24, 15, 10, 0 ), 
             fdtsec=60, height=3000.0 ):

    ftsec_ = ( vtime - fstime ).total_seconds()

    tlev = int( ftsec_ / fdtsec )
    print( tlev )

    nc = Dataset( fn, "r", format="NETCDF4" )
    var4d = nc.variables['Reflectivity'][:]
    lon1d = nc.variables['Longitude'][:]
    lat1d = nc.variables['Latitude'][:]
    z1d = nc.variables['Height'][:]
    ftimes = nc.variables['time'][:]
    nc.close()

    zlev = np.argmin( np.abs( z1d - height ) )

    return( var4d[tlev,:,:,zlev], lon1d, lat1d )

def main( INFO={}, vtime=datetime( 2019, 8, 24, 15, 10, 0), 
          fstime=datetime( 2019, 8, 24, 15, 10, 0 ) ):

    dz = 300

    ofn = '{0:}/pawr_grads/pawr_obs_{1:}.nc'.format( INFO["TOP"],
                 vtime.strftime('%Y%m%d-%H%M%S') )
    odat, olon, olat = read_obs_nc( fn=ofn, dz=dz )
    print( odat.shape )

    ffn = '{0:}/dafcst_nc/{1:}.nc'.format( INFO["TOP"],
                 fstime.strftime('%Y%m%d-%H%M%S') )
    print( ffn )
    fdat2d, flon1d, flat1d = read_nc( fn=ffn, fstime=fstime, vtime=vtime )
 
    flon2d, flat2d = np.meshgrid( flon1d, flat1d )

    lons = np.min( flon2d )
    lone = np.max( flon2d )

    lats = np.min( flat2d )
    late = np.max( flat2d )

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
                            fs=8, lw=0.0 )

       if i == 0:
          SHADE = ax.contourf( flon2d, flat2d, fdat2d[:,:], 
                       norm=norm,
                       levels=levels,
                       cmap=cmap,
                       extend='both',
                       transform=data_crs )
       elif i == 1:
          SHADE = ax.scatter( olon, olat, s=5.0, c=odat, 
                      norm=norm, 
                      cmap=cmap,
                      #extend='both',
                      #vmin=np.min( levels ), 
                      #vmax=np.max( levels ), 
                      transform=data_crs ) 
          print( np.max( odat ) )

       pos = ax.get_position()
       cb_width = 0.006
       cb_height = pos.height*0.9
       ax_cb = fig.add_axes( [ pos.x1+0.002, pos.y0, 
                               cb_width, cb_height] )
       cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                          ticks=levels[::1], extend='both' )

#    levels = np.arange( 0, 6500, 500 )
    plt.show()

#    sys.exit()
#    ref4d = read_nc( fn=fn )
#
#    print( ref4d.shape )
#    import matplotlib.pyplot as plt
#    plt.contourf( ref4d[10,:,:,18], levels=[15, 20, 30 ] )
#    plt.show()

INFO = { "TOP": "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test/d4" }

vtime = datetime( 2019, 8, 24, 15, 30, 0 )
fstime = datetime( 2019, 8, 24, 15, 10, 0 )


main( INFO=INFO, vtime=vtime, fstime=fstime )
