from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import sys
import os

from tools_AIP import prep_proj_multi_cartopy, get_cfeature, setup_grids_cartopy, read_fcst_nc
import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm
#from matplotlib import cm

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh


quick = True
#quick = False


def main( INFO={}, vtime=datetime( 2019, 8, 24, 15, 10, 0), height=3000,
          fstime=datetime( 2019, 8, 24, 15, 10, 0 ), exp_l=[],
          clat=36.1, clon=136,
          texp_l=[], nvar='Reflectivity' ):


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

    if clon < 0.0:
       clon = lon_r

#    fig = plt.figure( figsize=(10, 5) )
#    fig.subplots_adjust( left=0.06, bottom=0.03, right=0.95, top=0.95,
#                         wspace=0.2, hspace=0.01)
#
#    ax_l = prep_proj_multi_cartopy( fig, xfig=2, yfig=1, proj='merc', 
#                         latitude_true_scale=lat_r )
 
    fig = plt.figure( figsize=(10, 9) )
    fig.subplots_adjust( left=0.06, bottom=0.03, right=0.95, top=0.95,
                         wspace=0.1, hspace=0.2)
    projection = ccrs.Mercator( latitude_true_scale=lat_r, )
    ax_l = []
    xfig = 2
    yfig = 2
    for i in range( 1, 5 ):
       if i <= 2:
          ax_l.append( fig.add_subplot( yfig,xfig,i, projection=projection ) )
       else:
          ax_l.append( fig.add_subplot( yfig,xfig,i, projection=None ) )

    if nvar == 'Reflectivity':
       levels = np.array( [ 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 ] )
       cmap = mcolors.ListedColormap(['cyan', 'b', 'dodgerblue',
                                      'lime','yellow',
                                      'orange', 'red', 'firebrick', 'magenta',
                                       'purple'])
       cmap.set_over('gray', alpha=1.0)
       cmap.set_under('w', alpha=0.0)
       cmap.set_under('mistyrose', alpha=1.0 )
       cmap.set_bad( color='gray', alpha=0.5 )
       unit = 'Z (dBZ)'

    elif nvar == 'Vr':
       #levels = np.array( [  -15, -10, -5, -1, 1, 5, 10, 15] ) 
       levels = np.arange( -6, 6.6, 0.5 )
       levels = np.array( [ -15,  -12, -9, -6, -3, -2, -1, -0.5, 0.5, 1,  2, 3, 6, 9, 12, 15 ] ) 
       #levels = np.arange( -5, 5.5, 0.5 )
       cmap = plt.cm.get_cmap("RdBu_r")
       #cmap = plt.cm.get_cmap("Spectral_r")
       cmap.set_over('gray', alpha=1.0)
       cmap.set_under('aqua', alpha=1.0)
       cmap.set_bad( color='gray', alpha=0.5 )
       #cmap = plt.cm.get_cmap("Spectral_r")
       unit = r'Vr (m s$^{-1}$)'

#    norm = None
    norm = BoundaryNorm( levels, ncolors=cmap.N, clip=False)
#    norm = cm.colors.Normalize( vmax=np.max( levels ), vmin=np.min( levels ) )

    l_l = []
    b_l = []
    w_l = []
    h_l = []

    for i, ax in enumerate( ax_l ):

        if i <= 1:
#           #print( "TEST\n", ax.get_position[0] )
#           l, b, w, h = ax.get_position().bounds
#           l_l.append( l )
#           b_l.append( b )
#           w_l.append( w )
#           h_l.append( h )

           ffn = '{0:}/{1:}/dafcst_nc/{2:}.nc'.format( INFO["TOP"], exp_l[i],
                        fstime.strftime('%Y%m%d-%H%M%S') )
           print( ffn )
           fdat2d, flon1d, flat1d = read_fcst_nc( fn=ffn, fstime=fstime, vtime=vtime, nvar=nvar, height=height )
         
           flon2d, flat2d = np.meshgrid( flon1d, flat1d )
 
           rxlev = np.argmin( np.abs( flon1d - lon_r ) )
           rylev = np.argmin( np.abs( flat1d - lat_r ) )
           i1d = ( np.arange( 0, len( flon1d ), 1 ) - rxlev ) * 0.5
           j1d = ( np.arange( 0, len( flon1d ), 1 ) - rylev ) * 0.5
       
           i2d, j2d = np.meshgrid( i1d, j1d )
           dist2d = np.sqrt( np.square( i2d ) + np.square( j2d ) )

           lons = np.min( flon2d )
           lone = np.max( flon2d )
        
           lats = np.min( flat2d )
           late = np.max( flat2d )
    
           ax.set_extent([ lons, lone, lats, late ] )
           ax.add_feature( coast, zorder=0 )
     
           setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                                fs=9, lw=0.0 )
     
           SHADE = ax.contourf( flon2d, flat2d, fdat2d[:,:], 
                        norm=norm,
                        levels=levels,
                        cmap=cmap,
                        extend='both',
                        transform=data_crs )

           ax.plot( lon_r, lat_r, marker='o', color='k', markersize=10,
                    transform=data_crs )

           CONT = ax.contour( flon2d, flat2d, dist2d, 
                        levels=[10, 20, 30, 40,  50, 60 ],
                        colors='gray',
                        linestyles='dashed',
                        linewidths=1.0, 
                        transform=data_crs,
                            )
           if clat > 0.0:
              ax.plot( [ lons, lone ], [ clat, clat ], transform=data_crs,
                       color='k', ls='dashed', lw=1.0 )
           elif clon > 0.0:
              ax.plot( [ clon, clon ], [ lats, late ], transform=data_crs,
                       color='k', ls='dashed', lw=1.0 )

        else:

#           l, b, w, h = ax.get_position().bounds
#           ax.set_position( [ l_l[i-2], b, w_l[i-2], h_l[i-2] ] )

           ffn = '{0:}/{1:}/dafcst_nc/{2:}.nc'.format( INFO["TOP"], exp_l[i-2],
                        fstime.strftime('%Y%m%d-%H%M%S') )
           fdat2d, fx1d, fy1d = read_fcst_nc( fn=ffn, fstime=fstime, vtime=vtime , VERT=True, clon=clon, clat=clat, nvar=nvar )
           fy2d, fx2d = np.meshgrid( fy1d, fx1d )
           print( fx2d.shape, fy2d.shape, fdat2d.shape)
           SHADE = ax.contourf( fx2d, fy2d, fdat2d[:,:], 
                        norm=norm,
                        levels=levels,
                        cmap=cmap,
                        extend='both', )


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
           ax.text( 1.1, 1.01, 'Init at {0:}\nValid at {1:}'.format( cfstime, cvtime ),
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

           ax.text( 1.01, 0.95, unit, 
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


    if clon > 0.0:
       ref_ = clon
    else:
       ref_ = clat

    ofig = "4p_fcst_ac_{0:}_{1:}_h{2:}km_s{3:}_v{4:}_crs{5:.2f}.png".format( 
                    exp_l[0],
                    exp_l[1],
                    height*0.001, 
                    fstime.strftime('%Y%m%d%H%M%S'),
                    vtime.strftime('%Y%m%d%H%M%S'),
                    ref_, 
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
texp_l = [ "CTRL", "TEST", "CTRL", "TEST" ]

exp_l = [ "d4", "d4_500m_ac"]
texp_l = [ "CTRL", "TEST", "CTRL", "TEST" ]

#exp_l = [ "d4_500m_ac", "d4_500m_ac_vr"]
#texp_l = [ "TEST", "TEST (VR only)", "TEST", "TEST (VR only)" ]

#exp_l = [ "d4_500m_ac", "d4_500m_ac_z"]
#texp_l = [ "TEST", "TEST (Z only)", "TEST", "TEST (Z only)" ]

#exp_l = [ "d4_500m_ac_z", "d4_500m_ac_vr"]
#texp_l = [ "TEST (Z only)", "TEST (VR only)", "TEST (Z only)", "TEST (VR only)" ]

vtime = datetime( 2019, 8, 24, 15, 30, 0 )
#vtime = datetime( 2019, 8, 24, 15, 40, 0 )
#vtime = datetime( 2019, 8, 24, 15, 50, 0 )
#vtime = datetime( 2019, 8, 24, 16, 0, 0 )
fstime = datetime( 2019, 8, 24, 15, 30, 0 )
#fstime = datetime( 2019, 8, 24, 15, 20, 0 )

fstime = datetime( 2019, 8, 24, 15, 0, 30 )
fstime = datetime( 2019, 8, 24, 15, 1, 0 )
fstime = datetime( 2019, 8, 24, 15, 3, 0 )
fstime = datetime( 2019, 8, 24, 15, 5, 0 )
fstime = datetime( 2019, 8, 24, 15, 10, 0 )
fstime = datetime( 2019, 8, 24, 15, 7, 0 )
fstime = datetime( 2019, 8, 24, 15, 6, 0 )
fstime = datetime( 2019, 8, 24, 15, 5, 30 )
fstime = datetime( 2019, 8, 24, 15, 20, 0 )
fstime = datetime( 2019, 8, 24, 15, 30, 0 )

fstime = datetime( 2019, 8, 19, 13, 30, 0 )

#fstime = datetime( 2019, 8, 24, 15, 1, 30 )
vtime = fstime + timedelta( minutes=10 )
#vtime = fstime

height = 3000
#height = 1000
#height = 500

clat =  36.15
clat =  36.2
clat =  -1
#clat =  35.96
#clon =  36.2
clon = 139.37
clon = 139.5
#clon = 139.4

nvar = "Reflectivity"
#nvar = "Vr"

main( INFO=INFO, vtime=vtime, fstime=fstime, exp_l=exp_l, 
      texp_l=texp_l, height=height, clat=clat, clon=clon, nvar=nvar )
