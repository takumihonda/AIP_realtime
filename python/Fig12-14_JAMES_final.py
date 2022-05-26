import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat, dist, get_cfeature, setup_grids_cartopy, prep_proj_multi_cartopy, draw_rec_4p

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm

from scipy.interpolate import griddata

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

quick = True
quick = False

USE_ARCH_DAT = True
#USE_ARCH_DAT = False

def main( INFO, time_l=[], hgt=3000.0, tlev_l=[] ):

    # If directly produce a pdf, the file size becomes very large (~22 MB)
    # As a remedy, a temporary png file is converted to pdf on Mac
    if time_l[0] == datetime(2019, 8, 24, 15, 40, 0 ):
       data_path = "../../dat4figs_JAMES/Fig12"
       ofig = "Fig12.png"
    elif time_l[0] == datetime(2019, 8, 19, 13, 40, 0 ):
       data_path = "../../dat4figs_JAMES/Fig14"
       ofig = "Fig14.png"
    else:
       sys.exit()
    os.makedirs( data_path, exist_ok=True )

    # radar location
    lon_r = 139.609
    lat_r = 35.861

#    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    if not USE_ARCH_DAT:
       mz1d = INFO["obsz"]
       mzidx = np.argmin( np.abs( mz1d - hgt ) )
   
       mask, mlon2d, mlat2d = read_mask_full()
       mask2d = mask[mzidx,:,:]

#    fig.subplots_adjust( left=0.0, bottom=0.0, right=1.0, top=1.0,
#                         wspace=0.0, hspace=0.0 )
    fig = plt.figure( figsize=(13, 8.5) )
    fig.subplots_adjust( left=0.07, bottom=0.03, right=0.96, top=0.97,
                         wspace=0.02, hspace=0.02 )
 
    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    ax_l = prep_proj_multi_cartopy( fig, xfig=3, yfig=2, proj='merc', 
                         latitude_true_scale=lat_r )
 
    res = '10m'
#    if quick:
#       res = '50m'

    land = get_cfeature( typ='land', res=res )
    coast = get_cfeature( typ='coastline', res=res )

    if not USE_ARCH_DAT:
       time = datetime( 2019, 8, 24, 15, 0, 30 )
       obs3d, olon2d, olat2d, oz1d = read_obs_grads( INFO, itime=time )
       ozidx = np.argmin( np.abs( oz1d - hgt ) )
       mzidx = np.argmin( np.abs( mz1d - hgt ) )

       imask2d = griddata( ( mlon2d.ravel(), mlat2d.ravel() ), mask[mzidx,:,:].ravel(),
                          (flon2d, flat2d),
                          #method='cubic',
                          method='nearest',
                         )

#    # for pcolormesh
#    olon2d -= np.abs( olon2d[1,0] - olon2d[0,0] )
#    olat2d -= np.abs( olat2d[0,1] - olat2d[0,0] )
#    mlon2d -= np.abs( mlon2d[1,0] - mlon2d[0,0] )
#    mlat2d -= np.abs( mlat2d[0,1] - mlat2d[0,0] )

    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    fz1d = INFO["cz"]
    fzidx = np.argmin( np.abs( fz1d - hgt ) )

    i1d = np.arange( flon2d.shape[0] ) + 1.0
    j1d = np.arange( flon2d.shape[1] ) + 1.0

    i1d -= np.mean( i1d )
    j1d -= np.mean( j1d )

    # 0.5km mesh
    j2d, i2d = np.meshgrid( i1d*0.5, j1d*0.5 )

    dist2d = np.sqrt( np.square(i2d) + np.square(j2d) )
#    dist2d_ = dist( lon_r, lat_r, lon2d_4, lat2d_4 ) * 0.001


#    ox2d, oy2d = m_l[0]( olon2d, olat2d )
#    mx2d, my2d = m_l[0]( mlon2d, mlat2d )
#    fx2d, fy2d = m_l[0]( flon2d, flat2d )

#    x2d_, y2d_ = m_l[0]( lon2d_4, lat2d_4 )

    levs_dbz= np.array( [ 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 ] )
    cmap_dbz = mcolors.ListedColormap(['cyan', 'b', 'dodgerblue',
                                       'lime','yellow',
                                       'orange', 'red', 'firebrick', 'magenta',
                                       'purple'])
    cmap_dbz.set_over('gray', alpha=1.0)
    cmap_dbz.set_under('w', alpha=0.0)
    cmap_dbz.set_bad( color='gray', alpha=0.5 )

    norm = BoundaryNorm( levs_dbz, ncolors=cmap_dbz.N, clip=False)

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)" ]

    if not USE_ARCH_DAT:
       # for pcolor mesh
       xlen2 = flon2d.shape[0] // 2
       ylen2 = flon2d.shape[1] // 2
       xlen = flon2d.shape[0] 
       ylen = flon2d.shape[1] 
       x2d = flon2d - ( flon2d[xlen2+1,ylen2] - flon2d[xlen2,ylen2] )
       y2d = flat2d - ( flat2d[xlen2,ylen2+1] - flat2d[xlen2,ylen2] )


    lons = flon2d[0,0]
    lone = flon2d[-2,-2]

    lats = flat2d[0,0]
    late = flat2d[-2,-2]
 
#    halo = 8
#    lons = flon2d[halo,halo]
#    lone = flon2d[-halo,-halo]

#    lats = flat2d[halo,halo]
#    late = flat2d[-halo,-halo]

    print( lons, lone, lats, late )

    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )

    for i, ax in enumerate( ax_l ):
       #ax.set_adjustable('datalim')
       ax.set_aspect('auto')
       itime = time_l[i]
       tlev = tlev_l[i]

       fn = '{0:}/data{1:}.npz'.format( data_path, i )

       ax.set_extent([ lons, lone, lats, late ] )
#       ax.add_feature( land, zorder=0 )
       ax.add_feature( coast, zorder=0 )

       xfs = 12
       yfs = 12
       if i % 3 != 0:
          yfs = 0
       if i <= 2:
          xfs = 0
       setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                            xfs=xfs, yfs=yfs, lw=0.0 )

#       ax.add_feature(cfeature.LAND, color='g') 
#       ax.add_feature(cfeature.COASTLINE, linewidth=10.8)
#       ax.coastlines( color='k', linestyle='solid', linewidth=10.5, zorder=1 )
      
       if not USE_ARCH_DAT:
          if i<= 2: 
             obs3d, _, _, _ = read_obs_grads( INFO, itime=itime )
             obs2d_ = griddata( ( olon2d.ravel(), olat2d.ravel() ), 
                                obs3d[ozidx,:,:].ravel(),
                                (flon2d, flat2d),
                                #method='cubic',
                                method='nearest',
                               )
   
             var2d = np.where( ( imask2d < 1.0 ) , obs2d_, np.nan )
          else:
             print( "fcst", itime, tlev )
             fcst3d, _ = read_fcst_grads( INFO, itime=itime, tlev=tlev , FT0=True, )
             var2d = fcst3d[fzidx,:,: ]

          data = var2d[:xlen-1,:ylen-1]
          np.savez( fn, x2d=x2d, y2d=y2d, data=data )
       else:

          data = np.load( fn )['data']
          x2d = np.load( fn )['x2d']
          y2d = np.load( fn )['y2d']

       #SHADE = ax.pcolormesh( x2d, y2d, var2d[:xlen-1,:ylen-1], 
       SHADE = ax.pcolormesh( x2d, y2d, data, 
                       cmap=cmap_dbz, vmin=np.min(levs_dbz),
                       vmax=np.max(levs_dbz),
#                       norm=norm, 
                       transform=data_crs, 
                       )

       ax.plot( lon_r, lat_r, ms=8.0, marker='o', color='r',
                 markeredgecolor='w', transform=data_crs, )

       CONT = ax.contour( flon2d, flat2d, dist2d, 
                          levels=[20, 40, 60], zorder=1,
                          colors='k', linewidths=0.5,
                          linestyles='dashed',
                          transform=data_crs,
                          )

       ax.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
                   fontsize=12, fmt='%.0f km', colors="k" )

       if i == 5:
          pos = ax.get_position()
          cb_width = 0.006
          cb_height = pos.height*1.0
          ax_cb = fig.add_axes( [ pos.x1+0.002, pos.y1-cb_height*0.5, 
                                  cb_width, cb_height] )
          cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                             ticks=levs_dbz[::1], extend='both' )
          cb.ax.tick_params( labelsize=12 )

          ax.text( 1.01, 1.51, "(dBZ)",
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=12, )


       ax.text( 0.01, 0.99, pnum_l[i],
               va='top', 
               ha='left',
               transform=ax.transAxes,
               color='k', fontsize=14, 
               bbox=bbox )

       if i <= 2:
          ax.text( 0.5, 1.01, itime.strftime('%H:%M:%S UTC %m/%d') ,
                  va='bottom', 
                  ha='center',
                  transform=ax.transAxes,
                  color='k', fontsize=14, )

          if i == 2:
             ax.text( 0.9, 1.01, "Z={0:.0f} km".format( hgt/1000 ),
                     va='bottom', 
                     ha='left',
                     transform=ax.transAxes,
                     color='k', fontsize=14, )

       if i <= 2:
          tit = "MP-PAWR obs"
       else:
          #tit = "Forecast"
          tit = "Forecast\n(FT={0:.0f} min)".format( tlev*30/60 )
   
       ax.text( 0.5, 0.99, tit,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=16, 
               bbox=bbox )

       lon_l = [ 139.35, 139.9 ]
       lat_l = [ 35.9, 36.22 ]

       lon_l = [ 139.2, 139.7 ]
       lat_l = [ 35.7, 36.1 ]
#       draw_rec_4p( ax, lon_l=lon_l, lat_l=lat_l, lc='magenta', lw=2.0, 
#                    transform=data_crs )

#    ofig = "6p_obs_fcst_" + itime.strftime('%m%d') + ".png"
    print(ofig)

#    ofig = "test.png"

    if not quick:
       opath = "pdf/png"
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       plt.show()





############3

data_path = "../../dat4figs_JAMES/info"
os.makedirs( data_path, exist_ok=True )
fn_info = '{0:}/data.npz'.format( data_path, )


TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"
EXP = "20201117/D4_500m_CTRL"

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/dafcst".format( EXP )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

fcst_zmax = 43

if not USE_ARCH_DAT:
   obsz, olon2d, olat2d = read_obs_grads_latlon()
   lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )
   np.savez( fn_info, obsz=obsz, olon2d=olon2d, olat2d=olat2d,
                      lon2d=lon2d, lat2d=lat2d, hgt3d=hgt3d,
                      cz=cz, ohgt3d=ohgt3d,
            )
else:
   obsz = np.load( fn_info )['obsz']
   olon2d = np.load( fn_info )['olon2d']
   olat2d = np.load( fn_info )['olat2d']
   lon2d = np.load( fn_info )['lon2d']
   lat2d = np.load( fn_info )['lat2d']
   hgt3d = np.load( fn_info )['hgt3d']
   cz = np.load( fn_info )['cz']
   ohgt3d = np.load( fn_info )['ohgt3d']

INFO = { "TOP": TOP,
         "EXP": EXP,
         "time0": time0,
         "FCST_DIR": FCST_DIR,
         "gz": fcst_zmax,
         "gy": lon2d.shape[0],
         "gx": lon2d.shape[1],
         "lon2d": lon2d,
         "lat2d": lat2d,
         "cz": cz,
         "olon2d": olon2d,
         "olat2d": olat2d,
         "obsz": obsz,
       }


itime = datetime( 2019, 8, 19, 13, 30 )
#itime = datetime( 2019, 8, 24, 15, 30 )

tlev1 = 20
tlev2 = 40
tlev3 = 60

#tlev1 = 0
#tlev2 = 10
#tlev3 = 20

time_l = [
          itime + timedelta( seconds=tlev1*30 ),
          itime + timedelta( seconds=tlev2*30 ),
          itime + timedelta( seconds=tlev3*30 ),
          itime, # scale
          itime, # scale
          itime, # scale
         ]

hgt = 3000.0

tlev_l = [ tlev1, tlev2, tlev3,
           tlev1, tlev2, tlev3, ]

main( INFO, time_l=time_l, hgt=hgt, tlev_l=tlev_l )

