from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import sys
import os

from tools_AIP import prep_proj_multi_cartopy, get_cfeature, setup_grids_cartopy, read_fcst_nc, read_fcst_grads_all, get_cz
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

def set_cmap_fcst_grads_all( nvar='qg' ):

    if nvar == "qr" or nvar == "qg":
       levels = np.arange( 0.2, 3.0, 0.2 )
 
       cmap = plt.cm.get_cmap( 'plasma_r' )

       fac = 1.e3
       unit = r'g kg$^{-1}$'
       extend = 'max'
       tvar = nvar.capitalize()

    elif nvar == "w" or nvar == "u" or nvar == "v":
       #levels = np.arange( 0.2, 3.0, 0.2 )
       levels = np.array( [ -8, -6, -4, -2, -1.0, -0.5, 0.5, 1.0, 2, 4, 6, 8 ] )
 
       cmap = plt.cm.get_cmap( 'RdBu_r' )

       fac = 1.0
       unit = r'm s$^{-1}$'
       extend = 'both'
       tvar = nvar.upper()
       cmap.set_over('gray', alpha=1.0)
       cmap.set_under('k', alpha=1.0)

    norm = BoundaryNorm( levels, ncolors=cmap.N, clip=False)
#    cmap.set_bad( color='gray', alpha=0.5 )
    
    return( cmap, levels, norm, fac, unit, extend, tvar )   

def main( INFO={}, vtime=datetime( 2019, 8, 24, 15, 10, 0), height=3000,
          fstime=datetime( 2019, 8, 24, 15, 10, 0 ), exp_l=[],
          clat=36.1, clon=136, dtsec=30,
          texp_l=[], nvar='qg', ZSKIP=2 ):

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }


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

    cmap, levels, norm, fac, unit, extend, tvar = set_cmap_fcst_grads_all( nvar=nvar )

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


 
    l_l = []
    b_l = []
    w_l = []
    h_l = []

    # get grid data
    ffn = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test/d4_500m_ac/dafcst_nc/20190824-153000.nc'
    _, flon1d, flat1d = read_fcst_nc( fn=ffn, fstime=fstime, vtime=fstime )
    z1d_ = get_cz()
    z1d = []
    for i in range( 0, len( z1d_ ), ZSKIP ):
        z1d.append( z1d_[i] )
    z1d = np.array( z1d )

    zlev = np.argmin( np.abs( z1d - height ) )
    xlev = np.argmin( np.abs( flon1d - clon ) )
    ylev = np.argmin( np.abs( flat1d - clat ) )


    flon2d, flat2d = np.meshgrid( flon1d, flat1d )
        
    lons = np.min( flon2d )
    lone = np.max( flon2d )
        
    lats = np.min( flat2d )
    late = np.max( flat2d )
    
    for i, ax in enumerate( ax_l ):

        if i <= 1:
#           #print( "TEST\n", ax.get_position[0] )
#           l, b, w, h = ax.get_position().bounds
#           l_l.append( l )
#           b_l.append( b )
#           w_l.append( w )
#           h_l.append( h )

           fn = '{0:}/{1:}/dafcst/fcst_all_{2:}.grd'.format( INFO["TOP"], exp_l[i],
                        fstime.strftime('%Y%m%d-%H%M%S') )
           #print( fn )
           tlev = int( ( vtime - fstime ).total_seconds() / dtsec )
           fdat3d = read_fcst_grads_all( INFO, fn=fn, itime=fstime, tlev=tlev, nvar=nvar, gz=len( z1d )  )
           fdat2d = fdat3d[zlev,:,:]*fac
         
           ax.set_extent([ lons, lone, lats, late ] )
           ax.add_feature( coast, zorder=0 )
     
           setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                                fs=9, lw=0.0 )
     
           SHADE = ax.contourf( flon2d, flat2d, fdat2d[:,:], 
                        norm=norm,
                        levels=levels,
                        cmap=cmap,
                        extend=extend,
                        transform=data_crs )

           ax.plot( lon_r, lat_r, marker='o', color='k', markersize=10,
                    transform=data_crs )

           if clat > 0.0:
              ax.plot( [ lons, lone ], [ clat, clat ], transform=data_crs,
                       color='k', ls='dashed', lw=1.0 )
           elif clon > 0.0:
              ax.plot( [ clon, clon ], [ lats, late ], transform=data_crs,
                       color='k', ls='dashed', lw=1.0 )

           ax.text( 0.99, 0.01, 'Z={0:.1f} km'.format( height*0.001 ), 
                   va='bottom', 
                   ha='right',
                   bbox=bbox,
                   transform=ax.transAxes,
                   color='k', fontsize=10, )

        else:

#           l, b, w, h = ax.get_position().bounds
#           ax.set_position( [ l_l[i-2], b, w_l[i-2], h_l[i-2] ] )

#           ffn = '{0:}/{1:}/dafcst_nc/{2:}.nc'.format( INFO["TOP"], exp_l[i-2],
#                        fstime.strftime('%Y%m%d-%H%M%S') )
#           fdat2d, fx1d, fy1d = read_fcst_nc( fn=ffn, fstime=fstime, vtime=vtime , VERT=True, clon=clon, clat=clat )

           fn = '{0:}/{1:}/dafcst/fcst_all_{2:}.grd'.format( INFO["TOP"], exp_l[i-2],
                        fstime.strftime('%Y%m%d-%H%M%S') )
           #print( fn )
           tlev = int( ( vtime - fstime ).total_seconds() / dtsec )
           fdat3d = read_fcst_grads_all( INFO, fn=fn, itime=fstime, tlev=tlev, nvar=nvar, gz=len( z1d )  )


           fx1d = z1d*1.e-3
           if clat > 0.0:
              fdat2d = fdat3d[:,ylev,:]*fac
              fy1d = flon1d

           elif clon > 0.0:

              fdat2d = fdat3d[:,:,xlev]*fac
              fy1d = flat1d
         

           fx2d, fy2d = np.meshgrid( fy1d, fx1d )
           print( fx2d.shape, fy2d.shape, fdat2d.shape)
           SHADE = ax.contourf( fx2d, fy2d, fdat2d[:,:], 
                        norm=norm,
                        levels=levels,
                        cmap=cmap,
                        extend=extend, )


        ax.text( 0.5, 1.01, texp_l[i],
                va='bottom', 
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=12, )

           

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

    tit = '{0:} SCALE-LETKF {1:}'.format( tvar, tit_ )
    fig.suptitle( tit, fontsize=14 )


    if clon > 0.0:
       ref_ = clon
    else:
       ref_ = clat

    ofig = "4p_fcst_ac_{0:}_h{1:}km_s{2:}_v{3:}_crs{4:.2f}.png".format( 
                    nvar,
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


INFO = { "TOP": "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test",
         "gx": 256,
         "gy": 256,
         "gz": 45,
        }

exp_l = [ "d4", "d4_attenuation_corrected"]
texp_l = [ "CTRL", "TEST", "CTRL", "TEST" ]

exp_l = [ "d4", "d4_500m_ac"]
texp_l = [ "CTRL", "TEST", "CTRL", "TEST" ]


#exp_l = [ "d4", "d4_500m_ac"]
#texp_l = [ "CTRL", "TEST", "CTRL", "TEST" ]

exp_l = [ "d4_500m_ac", "d4_500m_ac_z"]
texp_l = [ "TEST", "TEST Z", "TEST", "TEST Z" ]

exp_l = [ "d4_500m_ac_z", "d4_500m_ac_vr"]
texp_l = [ "TEST (Z only)", "TEST (VR only)", "TEST (Z only)", "TEST (VR only)" ]

vtime = datetime( 2019, 8, 24, 15, 30, 0 )
#vtime = datetime( 2019, 8, 24, 15, 40, 0 )
#vtime = datetime( 2019, 8, 24, 15, 50, 0 )
#vtime = datetime( 2019, 8, 24, 16, 0, 0 )
fstime = datetime( 2019, 8, 24, 15, 30, 0 )
#fstime = datetime( 2019, 8, 24, 15, 20, 0 )

#vtime = fstime

height = 3000
height = 3000
height = 500
#height = 1000

clat =  36.15
clat =  36.2
clat =  36.06
clat =  -1
#clon =  36.2
clon = 139.37

nvar = "qg"
nvar = "qr"
nvar = "v"
#nvar = "w"

main( INFO=INFO, vtime=vtime, fstime=fstime, exp_l=exp_l, 
      texp_l=texp_l, height=height, clat=clat, clon=clon,
      nvar=nvar )
