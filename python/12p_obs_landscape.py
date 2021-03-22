import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat, dist, get_cfeature, setup_grids_cartopy, prep_proj_multi_cartopy, read_nowcast_hires, dbz2rain

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl

from scipy.interpolate import griddata

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

quick = True
#quick = False

def main( INFO, time_l=[], hgt=3000.0, tlev_l=[], lab_l=[], CRS=False, 
          clon=139.75, clat=36.080 ):

    # cross section
    clons = 139.5 + 0.001
    clone = 140.1 - 0.001

    clats = 35.85 + 0.001 + 0.1
    clate = 36.25 - 0.001

    # radar location
    lon_r = 139.609
    lat_r = 35.861

    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    mz1d, _, _ = read_obs_grads_latlon()
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    mask, mlon2d, mlat2d = read_mask_full()
    mask2d = mask[mzidx,:,:]

    fig = plt.figure( figsize=(13, 9) )
    fig.subplots_adjust( left=0.04, bottom=0.03, right=0.96, top=0.97,
                         wspace=0.2, hspace=0.06)
 
    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    ax_l = prep_proj_multi_cartopy( fig, xfig=4, yfig=3, proj='merc', 
                         latitude_true_scale=lat_r )
 
    res = '10m'
    if quick:
       res = '50m'

    land = get_cfeature( typ='land', res=res )
    coast = get_cfeature( typ='coastline', res=res )


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


    cmap = plt.cm.get_cmap("jet")
    levs = np.array( [ 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40 ] )

#    cmap = mpl.colors.ListedColormap(['cyan','dodgerblue',
#                                      'blue', 'yellow',
#                                      'orange', 'red', 'magenta'])
#    levs = np.array( [ 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0 ] )

    cmap.set_under( 'w', alpha=0.0 )
    cmap.set_over( 'k', alpha=1.0 )
    cmap.set_bad( color='gray', alpha=0.3 )
    extend = 'max'


    norm = BoundaryNorm( levs, ncolors=cmap.N, clip=False)

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    pnum_l = [ "(a)", "(b)", "(c)", 
               "(d)", "(e)", "(f)", 
               "(g)", "(h)", "(i)",
               "(j)", "(k)", "(l)",
               ]


    # for pcolor mesh
    fxlen2 = flon2d.shape[0] // 2
    fylen2 = flon2d.shape[1] // 2
    fxlen = flon2d.shape[0] 
    fylen = flon2d.shape[1] 
    fx2d = flon2d - ( flon2d[fxlen2+1,fylen2] - flon2d[fxlen2,fylen2] )
    fy2d = flat2d - ( flat2d[fxlen2,fylen2+1] - flat2d[fxlen2,fylen2] )

    # for pcolor mesh
    oxlen2 = olon2d.shape[0] // 2
    oylen2 = olon2d.shape[1] // 2
    oxlen = olon2d.shape[0] 
    oylen = olon2d.shape[1] 
    ox2d = olon2d - ( olon2d[oxlen2+1,oylen2] - olon2d[oxlen2,oylen2] )
    oy2d = olat2d - ( olat2d[oxlen2,oylen2+1] - olat2d[oxlen2,oylen2] )



    lons = flon2d[0,0] + 0.5
    lone = flon2d[-2,-2] - 0.15 - 0.1

    lats = flat2d[0,0] + 0.5
    late = flat2d[-2,-2] -0.05 - 0.1
 
    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )


    for i, ax in enumerate( ax_l ):
       itime = time_l[i]
       tlev = tlev_l[i]
       lab_ = lab_l[i]

       ax.set_extent([ lons, lone, lats, late ] )
       ax.add_feature( land, zorder=0 )
       ax.add_feature( coast, zorder=0 )

       setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                            fs=9, lw=0.0 )

#       ax.add_feature(cfeature.LAND, color='g') 
#       ax.add_feature(cfeature.COASTLINE, linewidth=10.8)
#       ax.coastlines( color='k', linestyle='solid', linewidth=10.5, zorder=1 )
      
       if lab_ == "obs":
          obs3d, _, _, _ = read_obs_grads( INFO, itime=itime )
          obs3d[ obs3d == -9.99e33 ] = np.nan 
#          obs2d_ = griddata( ( olon2d.ravel(), olat2d.ravel() ), 
#                             obs3d[ozidx,:,:].ravel(),
#                             (flon2d, flat2d),
#                             method='nearest',
#                            )

          #var2d = np.where( ( imask2d < 1.0 ) , obs2d_, np.nan )
          var2d = np.where( ( mask[mzidx,:,:] < 1.0 ) , obs3d[ozidx,:,:], np.nan )
          x2d = ox2d
          y2d = oy2d
          xlen = oxlen
          ylen = oylen

       elif lab_ == "scale":
          fcst3d, _ = read_fcst_grads( INFO, itime=itime, tlev=tlev , FT0=True, )
          var2d = fcst3d[fzidx,:,: ]
          x2d = fx2d
          y2d = fy2d
          xlen = fxlen
          ylen = fylen

       elif lab_ == "nowcast":
          jma2d, _, _ = read_nowcast_hires( stime=itime, ft=timedelta(seconds=tlev*30) )

          x2d = jx2d
          y2d = jy2d
          xlen = jxlen
          ylen = jylen
          var2d = jma2d

#          var2d = griddata( ( jlon2d.ravel(), jlat2d.ravel() ), 
#                            jma2d.ravel(),
#                            (flon2d, flat2d),
#                            method='nearest',
#                            )

    
       if lab_ != "nowcast":
          var2d = dbz2rain( var2d )


       SHADE = ax.pcolormesh( x2d, y2d, var2d[:xlen-1,:ylen-1], 
                       cmap=cmap, vmin=np.min(levs),
                       vmax=np.max(levs),
                       norm=norm, 
                       transform=data_crs, #extend=extend,
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
                   fontsize=8, fmt='%.0f km', colors="k" )

       if CRS:

#          ax.plot( [ clons, clone ], [ clat, clat ], 
#                   color='r', linewidth=1.0, linestyle='dotted',
#                   transform=data_crs )

          ax.plot( [ clon, clon ], [ clats, clate ], 
                   color='r', linewidth=1.0, linestyle='dotted',
                   transform=data_crs )

#          ctit_l = [ "A", "B" ]
#          clon_l = [ clons, clone ]
#          dlon_l = [ -0.0, 0.0 ]
#          ha_l = [ "left", "right" ]
#          for j in range( 2 ):
#              ax.text( clon_l[j] + dlon_l[j], clat, ctit_l[j],
#                       fontsize=10,
#                       ha=ha_l[j],
#                       va='bottom',
#                       color='r',
#                       bbox=bbox,     
#                       transform=data_crs )

#          ctit_l = [ "C", "D" ]
          ctit_l = [ "A", "B" ]
          clat_l = [ clats, clate ]
          dlat_l = [ -0.0, 0.0 ]
          va_l = [ "bottom", "top" ]
          for j in range( 2 ):
              print( clon )
              ax.text( clon, clat_l[j] + dlat_l[j], ctit_l[j], 
                       fontsize=10,
                       ha='left',
                       va=va_l[j], 
                       color='r',
                       bbox=bbox,     
                       transform=data_crs )


       if i == 7:
          pos = ax.get_position()
          cb_width = 0.008
          cb_height = pos.height*2.0
          ax_cb = fig.add_axes( [ pos.x1+0.005, pos.y1-cb_height*0.5, 
                                  cb_width, cb_height] )
          cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                             ticks=levs[::1], extend=extend )
          cb.ax.tick_params( labelsize=9 )

          ax.text( 1.01, 2.01, "(mm/h)",
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=9, )


       ax.text( 0.01, 0.99, pnum_l[i],
               va='top', 
               ha='left',
               transform=ax.transAxes,
               color='k', fontsize=10, 
               bbox=bbox )

       vtime = itime + timedelta( seconds=tlev*30 )
       ax.text( 0.5, 0.01, vtime.strftime('%H:%M:%S UTC') ,
               va='bottom', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=11,
               bbox=bbox, 
               zorder=4 )

       if i == 3:
          ax.text( 0.9, 1.01, "Z={0:.0f} km".format( hgt/1000 ),
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, )

       if lab_ == "obs":
          tit = "MP-PAWR obs"
       elif lab_ == "scale":
          tit = "Forecast (FT={0:.1f} min)".format( tlev*30/60 )
       elif lab_ == "nowcast":
          tit = "JMA nowcast (FT={0:.0f} min)".format( tlev*30/60 )
   
       ax.text( 0.5, 0.99, tit,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=11, 
               bbox=bbox )


    ofig = "12p_obs_{0:}_clon{1:.2f}_clat{2:.2f}_landscape_stime{3:}.png".format( itime.strftime('%m%d'), clon, clat, time_l[4].strftime('%m%d%H%M%S')  )
    print(ofig)

    if not quick:
       opath = "png"
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       plt.show()





############3

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"
EXP = "20201117/D4_500m_CTRL"

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/dafcst".format( EXP )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

fcst_zmax = 43

obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )

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
       }


itime = datetime( 2019, 8, 19, 13, 30 )
sitime = datetime( 2019, 8, 24, 15, 30 )

itime = datetime( 2019, 8, 24, 15, 25 )

stlev1 = 0
stlev2 = 10
stlev3 = 20
stlev4 = 30


time_l = [
          itime,
          itime,
          itime,
          itime,
          itime,
          itime,
          itime,
          itime,
          itime,
          itime,
          itime,
          itime,
         ]

di = 2
for i, time in enumerate( time_l ):
    time_l[i] = itime + timedelta( seconds=i*di*30 )


hgt = 2000.0

tlev_l = [ 0, 0, 0, 0, 
           0, 0, 0, 0,
           0, 0, 0, 0,
         ]


lab_p = "MP-PAWR obs"
lab_n = "JMA nowcast"
lab_s = "SCALE-LETKF Forecast"
lab_p = "obs"
lab_n = "nowcast"
lab_s = "scale"

lab_l = [ lab_p, lab_p, lab_p, lab_p,
          lab_p, lab_p, lab_p, lab_p,
          lab_p, lab_p, lab_p, lab_p,
        ]

clon = 139.75
clat = 36.080
clat = 36.09

CRS = True
main( INFO, time_l=time_l, hgt=hgt, tlev_l=tlev_l, lab_l=lab_l, CRS=CRS )

