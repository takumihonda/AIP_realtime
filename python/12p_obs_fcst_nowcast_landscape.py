import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat, dist, get_cfeature, setup_grids_cartopy, prep_proj_multi_cartopy, read_nowcast_hires, dbz2rain, draw_rec_4p, read_obs,read_mask

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

data_path = "../../dat4figs_GRL/Fig01"
os.makedirs( data_path, exist_ok=True )

USE_ARCH_DAT = True
#USE_ARCH_DAT = False

def main( INFO, time_l=[], hgt=3000.0, tlev_l=[], lab_l=[], CRS=False, 
          clon=139.75, clat=36.080 ):


    rec_lats = 36.05
    rec_late = 36.1
   
    rec_lons = 139.7
    rec_lone = 139.8

    # cross section
#    clons = 139.5 + 0.001
#    clone = 140.1 - 0.001
#
#    clats = 35.85 + 0.001 + 0.1
#    clate = 36.25 - 0.001
#  
#    clats = 35.951 
#    clate = 36.2

    clats = 35.951 + 0.05
    clate = 36.2   - 0.02

    # radar location
    lon_r = 139.609
    lat_r = 35.861

    fn = '{0:}/topo.npz'.format( data_path, )

    if USE_ARCH_DAT:
       topo2d_4 = np.load( fn )['topo']
       lat2d_4 = np.load( fn )['lat']
       lon2d_4 = np.load( fn )['lon']
       mask = np.load( fn )['mask']
    else:
       lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
       mask = read_mask()
       np.savez( fn, lon=lon2d_4, lat=lat2d_4, topo=topo2d_4, mask=mask )

    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    mz1d = INFO["obsz"]
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    #mask, mlon2d, mlat2d = read_mask_full()
    mlon2d = INFO["olon2d"]
    mlat2d = INFO["olat2d"]
    mask2d = mask[mzidx,:,:]

    fig = plt.figure( figsize=(14, 7) )
    fig.subplots_adjust( left=0.3, bottom=0.03, right=0.96, top=0.97,
                         wspace=0.05, hspace=0.01)
 
    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    xfig = 4
    yfig = 3
    ax_l = prep_proj_multi_cartopy( fig, xfig=xfig, yfig=yfig, proj='merc', 
                         latitude_true_scale=lat_r )
 
    res = '10m'
    if quick:
       res = '50m'

    land = get_cfeature( typ='land', res=res )
    coast = get_cfeature( typ='coastline', res=res )


    time = datetime( 2019, 8, 24, 15, 0, 30 )
    #obs3d, olon2d, olat2d, oz1d = read_obs_grads( INFO, itime=time )
    #obs3d = read_obs( utime=time, mask=mask )
    olon2d = INFO["olon2d"]
    olat2d = INFO["olat2d"]
    oz1d = INFO["obsz"]
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


    # for pcolor mesh
    _, jlon2d, jlat2d = read_nowcast_hires( stime=time_l[0], ft=timedelta(seconds=0) )
    jxlen2 = jlon2d.shape[0] // 2
    jylen2 = jlon2d.shape[1] // 2
    jxlen = jlon2d.shape[0] 
    jylen = jlon2d.shape[1] 
    jx2d = jlon2d - ( jlon2d[jxlen2+1,jylen2] - jlon2d[jxlen2,jylen2] )
    jy2d = jlat2d - ( jlat2d[jxlen2,jylen2+1] - jlat2d[jxlen2,jylen2] )


    lons = flon2d[0,0] + 0.5
    lone = flon2d[-2,-2] - 0.15 - 0.2

    lats = flat2d[0,0] + 0.5 + 0.05
    late = flat2d[-2,-2] -0.05 - 0.08
 
    print( lons, lone )
    print( lats, late )

    lons = 139.40159606933594 
    lone = 139.9608367919922

    lats = 35.84
    late = 36.22

    lons = 139.4
    lone = 139.9

    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )


    for i, ax in enumerate( ax_l ):

       fn = '{0:}/data{1:}.npz'.format( data_path, i )

       itime = time_l[i]
       tlev = tlev_l[i]
       lab_ = lab_l[i]

       ax.set_extent([ lons, lone, lats, late ] )
       ax.add_feature( land, zorder=0 )
       ax.add_feature( coast, zorder=0 )

       yfs = 0
       xfs = 0
       if i % xfig == 0:
          yfs = 10
       if i >= (xfig*2):
          xfs = 10
#       yfs = 10
#       xfs = 10
       setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                            xfs=xfs, yfs=yfs, lw=0.0 )

#       ax.add_feature(cfeature.LAND, color='g') 
#       ax.add_feature(cfeature.COASTLINE, linewidth=10.8)
#       ax.coastlines( color='k', linestyle='solid', linewidth=10.5, zorder=1 )
      
       if lab_ == "obs":
          if not USE_ARCH_DAT:
             obs3d, _ = read_obs( utime=itime, mask=mask )
             np.savez( fn, obs3d=obs3d )
          else:
             obs3d = np.load( fn )['obs3d']
#          if not USE_ARCH_DAT:
#          #obs3d, _, _, _ = read_obs_grads( INFO, itime=itime )
#             obs3d, _ = read_obs( utime=itime, mask=mask )
#             np.savez( fn, obs3d=obs3d )
#          else:
#             obs3d = np.load( fn )['obs3d']

#          obs3d[ obs3d == -9.99e33 ] = np.nan 
#          obs2d_ = griddata( ( olon2d.ravel(), olat2d.ravel() ), 
#                             obs3d[ozidx,:,:].ravel(),
#                             (flon2d, flat2d),
#                             method='nearest',
#                            )

          #var2d = np.where( ( imask2d < 1.0 ) , obs2d_, np.nan )
          print( obs3d.shape )
          print( mask.shape,  )
          var2d = np.where( ( mask[mzidx,:,:] < 1.0 ) , obs3d[ozidx,:,:], np.nan )
          x2d = ox2d
          y2d = oy2d
          xlen = oxlen
          ylen = oylen

       elif lab_ == "scale":
          if not USE_ARCH_DAT:
             fcst3d, _ = read_fcst_grads( INFO, itime=itime, tlev=tlev , FT0=True, )
             np.savez( fn, fcst3d=fcst3d )
          else:
             fcst3d = np.load( fn )['fcst3d']
          var2d = fcst3d[fzidx,:,: ]
          x2d = fx2d
          y2d = fy2d
          xlen = fxlen
          ylen = fylen

       elif lab_ == "nowcast":
          if not USE_ARCH_DAT:
             jma2d, _, _ = read_nowcast_hires( stime=itime, ft=timedelta(seconds=tlev*30) )
             np.savez( fn, jma2d=jma2d )
          else:
             jma2d = np.load( fn )['jma2d']

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
                   fontsize=10, fmt='%.0f km', colors="k" )

       if CRS and ( i <= 2 or ( i >= 4 and i <= 6 ) ):

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
              ax.text( clon+0.01, clat_l[j] + dlat_l[j], ctit_l[j], 
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
               color='k', fontsize=12, 
               bbox=bbox )

       vtime = itime + timedelta( seconds=tlev*30 )
       if i < xfig:
          ax.text( 0.5, 1.01, vtime.strftime('%H:%M:%S UTC') ,
                  va='bottom', 
                  ha='center',
                  transform=ax.transAxes,
                  color='k', fontsize=13,
                  #bbox=bbox, 
                  zorder=4 )


       if i == 3:
          ax.text( 0.9, 1.01, "Z={0:.0f} km".format( hgt/1000 ),
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, )


       if i % xfig == 0:
          if i <= 9:
             lon_l = [ rec_lons, rec_lone ]
             lat_l = [ rec_lats, rec_late ]
             draw_rec_4p( ax, lon_l=lon_l, lat_l=lat_l, 
                    lc='magenta', lw=2.0, transform=data_crs )
   
   

          if lab_ == "obs":
             tit = "MP-PAWR\n obs"
          elif lab_ == "scale":
             tit = "SCALE-LETKF\nforecast"
          elif lab_ == "nowcast":
             tit = "JMA\nnowcast"

          ax.text( -0.7, 0.5, tit,
                  va='center', 
                  ha='left',
                  transform=ax.transAxes,
                  weight='bold',
#                  rotation=90,
                  color='k', fontsize=14, )

       if lab_ == "scale" or lab_ == "nowcast":
          if lab_ == "scale":
             tit = "FT={0:.1f} min".format( tlev*30/60 )
          elif lab_ == "nowcast":
             tit = "FT={0:.1f} min".format( tlev*30/60 )
      
          ax.text( 0.5, 0.98, tit,
                  va='top', 
                  ha='center',
                  transform=ax.transAxes,
                  color='k', fontsize=10, 
                  bbox=bbox )



#    ofig = "12p_obs_fcst_nowcast_{0:}_clon{1:.2f}_clat{2:.2f}_landscape_stime{3:}.png".format( itime.strftime('%m%d'), clon, clat, time_l[4].strftime('%m%d%H%M%S')  )
#    print(ofig)
    ofig = "Fig01_GRL.png"

    if not quick:
       opath = "pdf_GRL"
       os.makedirs( opath, exist_ok=True )
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


data_path_info = "../../dat4figs_GRL/info"
os.makedirs( data_path_info, exist_ok=True )
fn_info = '{0:}/data.npz'.format( data_path, )
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



#obsz, olon2d, olat2d = read_obs_grads_latlon()
#lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )

INFO = { "TOP": TOP,
         "EXP": EXP,
         "time0": time0,
         "FCST_DIR": FCST_DIR,
         "gz": fcst_zmax,
         "gy": lon2d.shape[0],
         "gx": lon2d.shape[1],
         "lon2d": lon2d,
         "lat2d": lat2d,
         "olon2d": olon2d,
         "olat2d": olat2d,
         "obsz": obsz,
         "cz": cz,
       }


itime = datetime( 2019, 8, 19, 13, 30 )
itime = datetime( 2019, 8, 24, 15, 30 )
sitime = datetime( 2019, 8, 24, 15, 30 )

#itime = datetime( 2019, 8, 24, 15, 25 )
#sitime = datetime( 2019, 8, 24, 15, 25 )

tlev1 = 0
tlev2 = 10
tlev3 = 20
tlev4 = 30

stlev1 = 0
stlev2 = 10
stlev3 = 20
stlev4 = 30

#sitime = datetime( 2019, 8, 24, 15, 25, 30 )
#stlev1 = 9
#stlev2 = 19
#stlev3 = 29
#stlev4 = 39

#sitime = datetime( 2019, 8, 24, 15, 26, 30 )
#stlev1 = 7
#stlev2 = 17
#stlev3 = 27
#stlev4 = 37

#sitime = datetime( 2019, 8, 24, 15, 27, 30 )
#stlev1 = 5
#stlev2 = 15
#stlev3 = 25
#stlev4 = 35

time_l = [
          itime + timedelta( seconds=tlev1*30 ),
          itime + timedelta( seconds=tlev2*30 ),
          itime + timedelta( seconds=tlev3*30 ),
          itime + timedelta( seconds=tlev4*30 ),
          sitime,  # SCALE
          sitime,  # SCALE
          sitime,  # SCALE
          sitime,  # SCALE
          itime,  # nowcast
          itime,  # nowcast
          itime,  # nowcast
          itime,  # nowcast
         ]



hgt = 2000.0

tlev_l = [ 0, 0, 0, 0, 
           stlev1, stlev2, stlev3, stlev4, 
           tlev1, tlev2, tlev3, tlev4, 
         ]


lab_p = "MP-PAWR obs"
lab_n = "JMA nowcast"
lab_s = "SCALE-LETKF Forecast"
lab_p = "obs"
lab_n = "nowcast"
lab_s = "scale"

lab_l = [ lab_p, lab_p, lab_p, lab_p,
          lab_s, lab_s, lab_s, lab_s,
          lab_n, lab_n, lab_n, lab_n,
        ]

clon = 139.75
clat = 36.080
clat = 36.09

CRS = True
main( INFO, time_l=time_l, hgt=hgt, tlev_l=tlev_l, lab_l=lab_l, CRS=CRS )

