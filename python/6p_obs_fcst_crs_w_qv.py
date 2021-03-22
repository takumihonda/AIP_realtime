import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat, dist, get_cfeature, setup_grids_cartopy, prep_proj_multi_cartopy, read_fcst_grads_all

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import FormatStrFormatter

from scipy.interpolate import griddata

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

quick = True
#quick = False

def main( INFO, time_l=[], hgt=3000.0, tlev_l=[], clat=40.0, clon=139.75,
          CRS="ZONAL", lab_l=[], ores='500m', nvar="w" ):

    if CRS == "ZONAL":
       clons = 139.5 + 0.001
       clone = 140.1 - 0.001
    elif CRS == "MERID":
       clon = 139.75
       clats = 35.85 + 0.001 + 0.1
       clate = 36.25 - 0.001

    # radar location


    # radar location
    lon_r = 139.609
    lat_r = 35.861

    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    mz1d, _, _ = read_obs_grads_latlon( ores=ores )
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    mask_, mlon2d, mlat2d = read_mask_full()
    mask = mask_[:22,:,:]
    mask2d = mask[mzidx,:,:]

    fig = plt.figure( figsize=(13, 8.5) )
    fig.subplots_adjust( left=0.04, bottom=0.03, right=0.96, top=0.97,
                         wspace=0.15, hspace=0.25 )
 
    # original data is lon/lat coordinate

    yfig = 2
    xfig = 3
    ax_l = []
    for i in range( 1, xfig*yfig+1 ):
       ax_l.append( fig.add_subplot( yfig,xfig,i, ) ) #projection=projection ) )

#    ax_l = prep_proj_multi_cartopy( fig, xfig=1, yfig=1, proj='merc', 
#                         latitude_true_scale=lat_r )
 
    res = '10m'
    if quick:
       res = '50m'

    land = get_cfeature( typ='land', res=res )
    coast = get_cfeature( typ='coastline', res=res )


    time = datetime( 2019, 8, 24, 15, 0, 30 )
    obs3d, olon2d, olat2d, oz1d = read_obs_grads( INFO, itime=time, ores=ores )
    ozidx = np.argmin( np.abs( oz1d - hgt ) )
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    olen2 = olat2d.shape[1] // 2
    olen2_ = olon2d.shape[0] // 2
    oyidx = np.argmin( np.abs( olat2d[:,olen2] - clat ) )
    oxidx = np.argmin( np.abs( olon2d[olen2_,:] - clon ) )

    flen2 = flat2d.shape[1] // 2
    flen2_ = flon2d.shape[0] // 2
    fyidx = np.argmin( np.abs( flat2d[:,flen2] - clat ) )
    fxidx = np.argmin( np.abs( flon2d[flen2_,:] - clon ) )



    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    fz1d = INFO["cz"][:INFO["gz"]]
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
    cmap_dbz.set_over('k', alpha=1.0)
    cmap_dbz.set_under('w', alpha=0.0)
    cmap_dbz.set_bad( color='gray', alpha=0.3 )

    cmap_qv = plt.cm.get_cmap("RdPu")
    cmap_qv.set_over('k', alpha=1.0)
    cmap_qv.set_under('w', alpha=0.0)
    #levs_qv = np.arange( 0, 19, 1)
    levs_qv = np.arange( 5, 19, 1)
    norm_qv = BoundaryNorm( levs_qv, ncolors=cmap_qv.N, clip=False)

    norm_dbz = BoundaryNorm( levs_dbz, ncolors=cmap_dbz.N, clip=False)

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

#    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)" ]
    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(e)" ]


#    lons = flon2d[0,0]
#    lone = flon2d[-2,-2]

    lats = flat2d[0,0]
    late = flat2d[-2,-2]
 
    lons = flon2d[0,0] + 0.5 
    lone = flon2d[-2,-2] - 0.15

    print( "Cross at lon: {0:.3f}-{1:.3f}".format( lons, lone ) )
    print( "Cross at lat: {0:.3f}-{1:.3f}".format( clat, clat ) )

    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )

    ymin = 0.0
    ymax = 9.0

#    xmin = lons + 0.05
#    xmax = lone - 0.05
    if CRS == "ZONAL":
       xmin = clons
       xmax = clone
    elif CRS == "MERID":
       xmin = clats
       xmax = clate
  

    levs_wu = np.arange( 1, 21, 1 )
    levs_wu = np.arange( 2, 22, 2 )
    levs_wd = np.arange( -21, 0, 1 )

    dfz1d = np.diff( fz1d ) * 0.5
    dfz1d = np.append( dfz1d, dfz1d[-1] )

    ylab = 'Height (km)'

    cmap = cmap_dbz
    levs = levs_dbz
    norm = norm_dbz

    for i, ax in enumerate( ax_l ):

       if i == 1 or i == 4:
          vtime = itime + timedelta( seconds=tlev*30 )
          ax.text( 0.5, 1.01, vtime.strftime('Valid: %H:%M:%S UTC %m/%d') ,
                  va='bottom', 
                  ha='center',
                  transform=ax.transAxes,
                  color='k', fontsize=13, )

       if i == 3:
          ax.axis( 'off' )
          continue

       itime = time_l[i]
       tlev = tlev_l[i]
       lab_ = lab_l[i]

#       ax.set_extent([ lons, lone, lats, late ] )
#       ax.add_feature( land, zorder=0 )
#       ax.add_feature( coast, zorder=0 )
#
#       setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
#                            fs=8, lw=0.0 )

#       ax.add_feature(cfeature.LAND, color='g') 
#       ax.add_feature(cfeature.COASTLINE, linewidth=10.8)
#       ax.coastlines( color='k', linestyle='solid', linewidth=10.5, zorder=1 )
      
       if lab_ == "obs": 
          obs3d, _, _, _ = read_obs_grads( INFO, itime=itime, ores=ores )
#          obs2d_ = griddata( ( olon2d.ravel(), olat2d.ravel() ), 
#                             obs3d[ozidx,:,:].ravel(),
#                             (flon2d, flat2d),
#                             #method='cubic',
#                             method='nearest',
#                            )
          obs3d[ obs3d == -9.99e33 ] = np.nan


          if CRS == "ZONAL":
             var2d = obs3d[:,oyidx,:]
             x2d, y2d = np.meshgrid( olon2d[oyidx,:] - ( olon2d[oyidx,1] - olon2d[oyidx,0] ), # pcolormesh
                                     oz1d - ( oz1d[1] - oz1d[0] ) * 0.5 )                     # pcolormesh
             xlen = x2d.shape[0] 
             ylen = x2d.shape[1] 
             var2d = np.where( ( mask[:,oyidx,:] < 1.0 ) , var2d, np.nan )

          elif CRS == "MERID":
             var2d = obs3d[:,:,oxidx]
             x2d, y2d = np.meshgrid( olat2d[:,oxidx] - ( olat2d[1,oxidx] - olat2d[0,oxidx] ), # pcolormesh
                                     oz1d - ( oz1d[1] - oz1d[0] ) * 0.5 )                     # pcolormesh

       else:
          print( "fcst", itime, tlev )
          fcst3d, _ = read_fcst_grads( INFO, itime=itime, tlev=tlev , FT0=True, )
          if i >= 4:
             nvar = "qv"
          fcst3d_nvar = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar=nvar ) 
          if CRS == "ZONAL":
             var2d = fcst3d[:,fyidx,: ]
             x2d, y2d = np.meshgrid( flon2d[fyidx,:], 
                                     fz1d - dfz1d )   # pcolormesh
          elif CRS == "MERID":
             var2d = fcst3d[:,:,fxidx ]
             var2d_c = fcst3d_nvar[:INFO["fcst_zmax"],:,fxidx ]
             x2d, y2d = np.meshgrid( flat2d[:,fxidx], 
                                     fz1d - dfz1d )   # pcolormesh
          
       if i >= 4:
          var2d = var2d_c * 1.e3
          cmap = cmap_qv
          levs = levs_qv
          norm = norm_qv

       xlen = x2d.shape[0] 
       ylen = x2d.shape[1] 


       SHADE = ax.pcolormesh( x2d, y2d*0.001, var2d[:xlen-1,:ylen-1], 
                       cmap=cmap, vmin=np.min(levs),
                       vmax=np.max(levs),
                       norm=norm, 
                       )
 
       if lab_ != "obs" and i <= 2: 
          CONT = ax.contour( x2d, y2d*0.001, var2d_c[:,:], 
                          levels=levs_wu, 
                          colors='k',
                          linewidths=1.0,     
                          linestyles='solid',     
                          )
#          ax.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
#                      fontsize=8, fmt='%.0f', colors="k" )

#          CONT2 = ax.contour( x2d, y2d*0.001, var2d_c[:,:], 
#                           levels=levs_wd, 
#                           colors='b',
#                           linewidths=1.0,     
#                           linestyles='solid',     
#                          )


       ax.set_ylim( ymin, ymax )
       ax.set_xlim( xmin, xmax )

       ctitx_l = [ 0.03, 0.97 ]

       if CRS == "ZONAL":
          ax.xaxis.set_major_formatter( FormatStrFormatter( '%.2fE' ) )
          ctit_l = [ "A", "B"]
       elif CRS == "MERID":
          ax.xaxis.set_major_formatter( FormatStrFormatter( '%.2fN' ) )
          #ctit_l = [ "C", "D"]
          ctit_l = [ "A", "B"]


       if i == 0 or i == 3:
          ax.set_ylabel( ylab, fontsize=10 )

#       ax.plot( lon_r, lat_r, ms=8.0, marker='o', color='r',
#                 markeredgecolor='w', transform=data_crs, )
#
#       CONT = ax.contour( flon2d, flat2d, dist2d, 
#                          levels=[20, 40, 60], zorder=1,
#                          colors='k', linewidths=0.5,
#                          linestyles='dashed',
#                          transform=data_crs,
#                          )
#
#       ax.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
#                   fontsize=8, fmt='%.0f km', colors="k" )

#       ax.plot( [ lons, lone ], [ clat, clat ], transform=data_crs, 
#                linewidth=5.0, color='r',
#                  )
  
       if i== 2 or i == 5:
          pos = ax.get_position()
          cb_width = 0.006
          cb_height = pos.height*0.9
          ax_cb = fig.add_axes( [ pos.x1+0.002, pos.y0, #-cb_height*0.5, 
                                  cb_width, cb_height] )
          cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                             ticks=levs[::1], extend='both' )
          cb.ax.tick_params( labelsize=8 )

          ax.text( 1.01, 0.99, "(dBZ)",
                  va='top', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=9, )


       ax.text( 0.01, 0.99, pnum_l[i],
               va='top', 
               ha='left',
               transform=ax.transAxes,
               color='k', fontsize=10, 
               bbox=bbox )

       for j in range( 2 ):
           ax.text( ctitx_l[j], 0.01, ctit_l[j],
                   va='bottom', 
                   ha='center',
                   transform=ax.transAxes,
                   color='r', 
                   fontsize=11, 
                   bbox=bbox )
    


       if i == 2:

          if CRS == "ZONAL":
             ctit_ = '{0:.2f}'.format( clat ) + 'N'
          elif CRS == "MERID":
             ctit_ = '{0:.2f}'.format( clon ) + 'E'
          ax.text( 0.9, 1.01, ctit_,
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, )

       if lab_ == "obs":
          tit = "MP-PAWR obs"
       else:
          tit = "Forecast (FT={0:.1f} min)".format( tlev*30/60 )
          vtime = itime + timedelta( seconds=tlev*30 )
          print( "vtime for SCALE", vtime )
   
       ax.text( 0.5, 0.99, tit,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=12, 
               bbox=bbox )

    cll = clat
    if CRS == "MERID":
       cll = clon

    ofig = "6p_obs_fcst_w_qv_crs_{0:}_{1:}_cll{2:.3f}_{3:}.png".format(  itime.strftime('%m%d'), CRS, cll, time_l[3].strftime('%m%d%H%M%S') )
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

obsz, olon2d, olat2d = read_obs_grads_latlon( )
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
         "fcst_zmax": fcst_zmax,
       }


itime = datetime( 2019, 8, 19, 13, 30 )
itime = datetime( 2019, 8, 24, 15, 30 )


tlev1 = 0
tlev2 = 10
tlev3 = 20



sitime = datetime( 2019, 8, 24, 15, 25, 30 )
stlev1 = 9
stlev2 = 19
stlev3 = 29

sitime = datetime( 2019, 8, 24, 15, 25, 0 )
stlev1 = 10
stlev2 = 20
stlev3 = 30

sitime = datetime( 2019, 8, 24, 15, 30 )
stlev1 = 0
stlev2 = 10
stlev3 = 20


time_l = [
          itime + timedelta( seconds=tlev1*30 ),
          itime + timedelta( seconds=tlev2*30 ),
          itime + timedelta( seconds=tlev3*30 ),
          sitime, # scale
          sitime, # scale
          sitime, # scale
         ]

hgt = 3000.0

tlev_l = [ 0, 0, 0, 
           stlev1, stlev2, stlev3, ]

lab_p = "obs"
lab_s = "scale"

lab_l = [ lab_p, lab_s, lab_s, 
          lab_p, lab_s, lab_s, 
        ]

tlev_l = [ 0, 0, 10, 
           0, 0, 0, ]

vtime1 = datetime( 2019, 8, 24, 15, 30, 0 )
vtime2 = datetime( 2019, 8, 24, 15, 25, 0 )
sitime1 = datetime( 2019, 8, 24, 15, 30, 0 )
sitime2 = datetime( 2019, 8, 24, 15, 25, 0 )

tlev_l = [ 0, 0, 10, 
           0, 10, 20, ]

vtime1 = datetime( 2019, 8, 24, 15, 30, 0 )
vtime2 = datetime( 2019, 8, 24, 15, 35, 0 )
sitime1 = datetime( 2019, 8, 24, 15, 30, 0 )
sitime2 = datetime( 2019, 8, 24, 15, 25, 0 )

time_l = [
          vtime1, 
          sitime1, # scale
          sitime2, # scale
          vtime2,
          sitime1, # scale
          sitime2, # scale
         ]


#lab_l = [ lab_s, lab_s, lab_s, 
#          lab_s, lab_s, lab_s, 
#        ]
#
#
#stime = datetime( 2019, 8, 24, 15, 30, 0 )
#vtime = datetime( 2019, 8, 24, 15, 30, 0 )
#time_l = [ stime, 
#           stime - timedelta( seconds=60*1 ), 
#           stime - timedelta( seconds=60*2 ), 
#           stime - timedelta( seconds=60*3 ), 
#           stime - timedelta( seconds=60*4 ), 
#           stime - timedelta( seconds=60*5 ), 
#         ]
#for i, stime_ in enumerate( time_l ):
#    tlev_l[i] = int( ( vtime - stime_ ).total_seconds() / 30 )
##    tlev_l[i] = 0

clon = 139.75
clat = 36.09

CRS = "ZONAL"
CRS = "MERID"

ores = "500m"
#ores = "100m"

nvar = "w"

main( INFO, time_l=time_l, hgt=hgt, tlev_l=tlev_l, clat=clat, CRS=CRS, 
      clon=clon, lab_l=lab_l, ores=ores, nvar=nvar )

