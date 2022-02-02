import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat, dist, get_cfeature, setup_grids_cartopy, prep_proj_multi_cartopy, read_ga_grads_all, draw_rec_4p, read_obs

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
quick = False

data_path = "../../dat4figs_GRL/Fig04"
os.makedirs( data_path, exist_ok=True )

USE_ARCH_DAT = False
USE_ARCH_DAT = True

def read_qh_grads_all( INFO, itime=datetime( 2019,8,24,15,30,0 ), typ='g'):

    qh = read_ga_grads_all( INFO, itime=itime, nvar="qc", typ=typ ) 
    for nvar in [ "qr", "qi", "qs", "qg" ]:
        qh += read_ga_grads_all( INFO, itime=itime, nvar=nvar, typ=typ ) 

    return( qh )

def main( INFO, time_l=[], hgt_l=[3000.0], clat=40.0, nvar_l=["w"], 
      CRS="ZONAL", clon=139.8, otime=datetime( 2019, 8, 24, 15, 30, 0 ) ):

    rec_lats = 36.05
    rec_late = 36.1
   
    rec_lons = 139.7
    rec_lone = 139.8

    data_crs = ccrs.PlateCarree()

    if CRS == "ZONAL":
       clons = 139.5 + 0.001
       clone = 140.1 - 0.001
    elif CRS == "MERID":
#       clats = 35.8 + 0.001
#       clate = 36.3 - 0.001
       clats = 35.85 + 0.001 + 0.1
       clate = 36.25 - 0.001

       clats = 35.951
       clate = 36.2

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
       mask_ = np.load( fn )['mask_']
    else:
       lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
       mask_, mlon2d, mlat2d = read_mask_full()
       np.savez( fn, lon=lon2d_4, lat=lat2d_4, topo=topo2d_4, mask_=mask_ )

#    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
#    mz1d, _, _ = read_obs_grads_latlon()
#    mzidx = np.argmin( np.abs( mz1d - hgt_l[i] ) )

    mask = mask_[:22,:,:]
#    mask2d = mask[mzidx,:,:]

    fig, ( (ax1,ax2,ax3),(ax4,ax5,ax6) ) = plt.subplots( 2, 3, 
                              gridspec_kw={
                              'width_ratios': [1, 1, 1],
                              'height_ratios': [1, 1.0]},
                              figsize=(12,9),
                             )
    ax_l = [ ax1, ax2, ax3, ax4, ax5, ax6 ]
    fig.subplots_adjust( left=0.06, bottom=0.1, right=0.99, top=0.96,
                         wspace=0.05, hspace=0.1 )


    time = datetime( 2019, 8, 24, 15, 0, 30 )
    #obs3d, olon2d, olat2d, oz1d = read_obs_grads( INFO, itime=otime )
    olon2d = INFO["olon2d"]
    olat2d = INFO["olat2d"]
    oz1d = INFO["obsz"]
    ozidx = np.argmin( np.abs( oz1d - hgt ) )

    obs3d, _ = read_obs( utime=otime, mask=mask )
   
#    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    olen2 = olat2d.shape[1] // 2
    oyidx = np.argmin( np.abs( olat2d[:,olen2] - clat ) )
    oxidx = np.argmin( np.abs( olon2d[olen2,:] - clon ) )

    flen2 = flat2d.shape[1] // 2
    fyidx = np.argmin( np.abs( flat2d[:,flen2] - clat ) )
    fxidx = np.argmin( np.abs( flon2d[flen2,:] - clon ) )



    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    fz1d = INFO["cz"][:INFO["gz"]]
#    fzidx = np.argmin( np.abs( fz1d - hgt ) )

    i1d = np.arange( flon2d.shape[0] ) + 1.0
    j1d = np.arange( flon2d.shape[1] ) + 1.0

    i1d -= np.mean( i1d )
    j1d -= np.mean( j1d )

    # 0.5km mesh
    j2d, i2d = np.meshgrid( i1d*0.5, j1d*0.5 )

    dist2d = np.sqrt( np.square(i2d) + np.square(j2d) )

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
    ox2d_ = olon2d - ( olon2d[oxlen2+1,oylen2] - olon2d[oxlen2,oylen2] )
    oy2d_ = olat2d - ( olat2d[oxlen2,oylen2+1] - olat2d[oxlen2,oylen2] )

    #levs_w = np.arange( -1, 1.2, 0.2)
    levs_w = [ -1., -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0] 
    levs_hdiv = np.arange( -5, 5.5, 0.5)

    cmap_w = plt.cm.get_cmap("RdBu_r")

    levs_u = np.arange( -1, 1.2, 0.2)
    cmap_u = plt.cm.get_cmap("RdBu_r")


    unit_ms = "(m/s)"

    levs_t = np.arange( -1, 1.2, 0.2 )

    cmap_t = plt.cm.get_cmap("RdBu_r")

    unit_t = "(K)"
    unit_hdiv = r'(10$^{-4}$s)'
    unit_qv = "(g/kg)"
    unit_qg = "(g/kg)"
    unit_qr = "(g/kg)"
    unit_qs = "(g/kg)"


    cmap_qv = plt.cm.get_cmap("BrBG")
    #levs_qv = np.arange( -1, 1.2, 0.2 )
    levs_qv = np.arange( -0.6, 0.7, 0.1 )
    levs_qv = np.arange( -0.3, 0.35, 0.05 )
    levs_qv = [ -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.05, 0.1, 0.15,  0.2, 0.25, 0.3 ]


    cmap_qg = cmap_qv
#    levs_qg = np.arange( -1, 1.2, 0.2 )
    levs_qg = levs_qv

    cmap_qr = cmap_qv
#    levs_qr = np.arange( -1, 1.2, 0.2 )
    levs_qr = levs_qv

    cmap_qs = cmap_qv
    levs_qs = np.arange( -1, 1.2, 0.2 )

    cmap_p = cmap_t
    levs_p = np.arange( -10, 12, 2 )

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)" ]


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


#    xmin = lons + 0.05
#    xmax = lone - 0.05

    dfz1d = np.diff( fz1d ) * 0.5
    dfz1d = np.append( dfz1d, dfz1d[-1] )

    ylab1 = 'Latitude'
    ylab2 = 'Height (km)'


    for i, ax in enumerate( ax_l ):

       fn = '{0:}/data{1:}.npz'.format( data_path, i )

       fzidx = np.argmin( np.abs( INFO["cz"][:] - hgt_l[i] ) )
       print( "chk", i )
     
       if CRS == "ZONAL":
          xmin = clons
          xmax = clone
       elif CRS == "MERID":
          xmin = clats
          xmax = clate

       if i <= 2:
#          xmin = 139.65
#          xmax = 139.8
#          ymin = 36
#          ymax = 36.15
          ymin = 35.99
          ymax = 36.21
          xmin = 139.62 #- 0.1
          xmax = 139.89 #+ 0.1
       else:
          ymin = 0.0 + 0.01
          ymax = 9.0 - 0.01

 
       nvar = nvar_l[i]

       if nvar == "w":
          fac = 1.0
          cmap = cmap_w
          levs = levs_w
          unit = unit_ms
          tvar = "W"

       elif nvar == "qv":
          fac = 1.e3
          tvar = "QV"
          cmap = cmap_qv
          levs = levs_qv
          unit = unit_qv


       elif nvar == "qh":
          fac = 1.e3
          tvar = "QHYD"
          cmap = cmap_qg
          levs = levs_qg
          unit = unit_qs


       cmap.set_over('k', alpha=1.0 )
       cmap.set_under('gray', alpha=1.0 )

       norm = BoundaryNorm( levs, ncolors=cmap.N, clip=False )

 
       if not USE_ARCH_DAT:
          for j, vtime in enumerate( time_l ):
             if j == 0:
                if nvar == "qh":
                   g3d = read_qh_grads_all( INFO, itime=vtime, typ='g')
                   a3d = read_qh_grads_all( INFO, itime=vtime, typ='a')
                elif nvar == "hdiv":
                   g3d = ( np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='u', typ='g' ), axis=2)  \
                         + np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='v', typ='g' ), axis=1) ) / ( 500.0*2 )
   
                   a3d = ( np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='u', typ='a' ), axis=2) \
                         + np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='v', typ='a' ), axis=1) ) / ( 500.0*2 ) 
   
                else:
                   g3d = read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='g' ) 
                   a3d = read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='a' ) 
             else:
                if nvar == "qh":
                   g3d += read_qh_grads_all( INFO, itime=vtime, typ='g')
                   a3d += read_qh_grads_all( INFO, itime=vtime, typ='a')
                elif nvar == "hdiv":
                   g3d += ( np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='u', typ='g' ), axis=2) \
                          + np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='v', typ='g' ), axis=1) ) / ( 500.0*2 )
   
                   a3d += ( np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='u', typ='a' ), axis=2) \
                          + np.gradient( read_ga_grads_all( INFO, itime=vtime, nvar='v', typ='a' ), axis=1) ) / ( 500.0*2 )
   
                else:
                   g3d += read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='g' ) 
                   a3d += read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='a' ) 
   
          var3d = ( a3d - g3d ) / len( time_l )
          np.savez( fn, var3d=var3d )

       else:
          var3d = np.load( fn )['var3d']
      

       if i <= 2:
          yfac = 1.0
          var2d = var3d[fzidx,:,:]
          x2d, y2d = fx2d, fy2d 

          ozidx = np.argmin( np.abs( oz1d - hgt_l[i] ) )
          print( "check\n", hgt_l[i], hgt, "\n" )
          ovar2d = obs3d[ozidx,:,:]
          ox2d, oy2d = ox2d_, oy2d_ 
       else:
          yfac = 1.e-3
          if CRS == "ZONAL":
             var2d = var3d[:,fyidx,: ]
             x2d, y2d = np.meshgrid( flon2d[fyidx,:], 
                                     fz1d - dfz1d )   # pcolormesh
   
          elif CRS == "MERID":
             var2d = var3d[:,:,fxidx ]
             x2d, y2d = np.meshgrid( fy2d[:,fxidx], 
                                     fz1d - dfz1d )   # pcolormesh
   
             ovar2d = obs3d[:,:,oxidx ]
             ox2d, oy2d = np.meshgrid( olat2d[:,oxidx] - ( olat2d[1,oxidx] - olat2d[0,oxidx] ), # pcolormesh
                                     oz1d - ( oz1d[1] - oz1d[0] ) * 0.5 )                     # pcolormesh
#             ox2d, oy2d = np.meshgrid( olat2d[:,oxidx], 
#                                       oz1d  )

       xlen = x2d.shape[0] 
       ylen = x2d.shape[1] 
       print( np.max( var2d ), np.min( var2d) )

       print( "debug", np.max( var2d*fac ), np.min( var2d*fac ), var2d.shape )
       SHADE = ax.pcolormesh( x2d, y2d*yfac, var2d[:xlen-1,:ylen-1]*fac, 
                       cmap=cmap, vmin=np.min(levs),
                       vmax=np.max(levs),
                       norm=norm, 
                       )
 
       # draw obs
#       CONT = ax.contour( ox2d, oy2d*yfac, ovar2d,
#                          levels=[ 20, ],
#                          linewidths=0.5, colors='k' )
       ovar2d[ ovar2d<20.0] = np.nan
       ocmap = mcolors.ListedColormap( ['w', 'gray', ] )
       OPLOT = ax.pcolormesh( ox2d, oy2d*yfac, ovar2d,
                             vmin=15, vmax=20, cmap=ocmap, alpha=0.2,
                             edgecolors='face', linewidths=0.0 )

       ctitx_l = [ 0.03, 0.97 ]

       if i >= 3:
          xticks = np.arange( 30, 40+0.05, 0.05 )
          if CRS == "ZONAL":
             ax.xaxis.set_major_formatter( FormatStrFormatter( '%.2fE' ) )
             ctit_l = [ "A", "B"]
          elif CRS == "MERID":
             ax.xaxis.set_major_formatter( FormatStrFormatter( '%.2fN' ) )
             #ctit_l = [ "C", "D"]
             ctit_l = [ "A", "B"]

          for j in range( 2 ):
              ax.text( ctitx_l[j], 0.1, ctit_l[j],
                      va='bottom', 
                      ha='center',
                      transform=ax.transAxes,
                      color='r', fontsize=11, 
                      bbox=bbox )
       else:
          xticks = np.arange( 130, 140+0.05, 0.05 )
          yticks = np.arange( 30, 40+0.05, 0.05 )
          ax.set_yticks( yticks )

          ax.plot( [ clon, clon ], [ clats, clate ],
                  color='r', linewidth=1.0, linestyle='dotted', )
                  #transform=data_crs )
          ax.xaxis.set_major_formatter( FormatStrFormatter( '%.2fE' ) )
          ax.yaxis.set_major_formatter( FormatStrFormatter( '%.2fN' ) )

          ctit_l = [ "A", "B" ]
          clat_l = [ clats, clate ]
          dlat_l = [ -0.0, 0.0 ]
          va_l = [ "bottom", "top" ]
          for j in range( 2 ):
              ax.text( clon+0.01, clat_l[j] + dlat_l[j], ctit_l[j], 
                       fontsize=10,
                       ha='left',
                       va=va_l[j], 
                       color='r',
                       bbox=bbox, )

          lon_l = [ rec_lons, rec_lone ]
          lat_l = [ rec_lats, rec_late ]
          #draw_rec_4p( ax, lon_l=lon_l, lat_l=lat_l, lc='magenta', lw=2.0, 
          lc = 'magenta'
          lw = 2.0
          ax.plot( [ lon_l[0], lon_l[0] ], [ lat_l[0], lat_l[1] ],  color=lc, lw=lw, )
          ax.plot( [ lon_l[1], lon_l[1] ], [ lat_l[0], lat_l[1] ],  color=lc, lw=lw, )
          ax.plot( [ lon_l[0], lon_l[1] ], [ lat_l[0], lat_l[0] ],  color=lc, lw=lw, )
          ax.plot( [ lon_l[0], lon_l[1] ], [ lat_l[1], lat_l[1] ],  color=lc, lw=lw, )

       ax.set_xticks( xticks )

       ax.set_ylim( ymin, ymax )
       ax.set_xlim( xmin, xmax )
       print( "loc", xmin, xmax, ymin, ymax )


       if i == 3: 
          ax.set_ylabel( ylab2, fontsize=10 )

 
       if i >= 3:
          pos = ax.get_position()
          cb_width = pos.width
          cb_height = 0.01
          ax_cb = fig.add_axes( [ pos.x0, pos.y0-cb_height-0.04,
                                  cb_width, cb_height] )
          cb = plt.colorbar( SHADE, cax=ax_cb, orientation='horizontal',  
                             ticks=levs[::1], extend='both' )
          cb.ax.tick_params( labelsize=7 )
 
          ax.hlines( y=hgt_l[i-3]/1000, xmin=xmin, xmax=xmax, ls='dashed',
                     lw=1.5, color='gray' ) 
        
#          ax.text( xmax-0.01, hgt_l[i-3]/1000+0.1, pnum_l[i-3],
#                  va='bottom', 
#                  ha='right',
#                  transform=ax.transData, 
#                  color='k', fontsize=10, 
#                  bbox=bbox )
       else:
          ax.text( 0.99, 0.99, 'Z={0:.1f} km'.format( hgt_l[i]/1000 ),
                  va='top', 
                  ha='right',
                  transform=ax.transAxes,
                  color='k', fontsize=12, 
                   )

       ax.text( 0.01, 0.99, pnum_l[i],
               va='top', 
               ha='left',
               transform=ax.transAxes,
               color='k', fontsize=12, 
               bbox=bbox )

       ax.text( 0.5, 0.98, tvar,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=13, 
               bbox=bbox )

       ax.text( 0.6, 0.98, unit,
               va='top', 
               ha='left',
               transform=ax.transAxes,
               color='k', fontsize=10, 
               bbox=bbox )

       if i != 0 and i != 3:
          ax.tick_params( axis='y', labelsize=0 )


       if i == 2:
#          if CRS == "ZONAL":
#             ctit_ = '{0:.2f}N'.format( clat )
#          elif CRS == "MERID":
#             ctit_ = '{0:.2f}E'.format( clon )
#          ax.text( 1.0, 1.08, ctit_,
#                  va='bottom', 
#                  ha='right',
#                  transform=ax.transAxes,
#                  color='k', fontsize=10, )
#
          trange = 'Period: {0:}-{1:}'.format( 
                      time_l[0].strftime('%H%M:%S'), 
                      time_l[-1].strftime('%H%M:%S UTC %m/%d/%Y'), 
                    )
          ax.text( 0.0, 1.01, trange,
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=12, 
                   )



#       tit = "Forecast (FT={0:.1f} min)".format( tlev*30/60 )
#   
#    tit = 'Time-averaged analysis increment'
#    fig.suptitle( tit, fontsize=14 ) 

    cll = clat
    if CRS == "MERID":
       cll = clon

    ofig = "6p_ainc_crs_{0:}_cll{1:.3f}.png".format( CRS, cll, )
    ofig = "Fig04_GRL.png"
    print(ofig)

    if not quick:
       opath = "png"
       opath = "pdf_GRL"
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
fn_info = '{0:}/data.npz'.format( data_path_info, )
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
         "olon2d": olon2d,
         "olat2d": olat2d,
         "obsz": obsz,
         "cz": cz,
       }



hgt = 3000.0
hgt = 4000.0
hgt = 500.0
#hgt = 200.0


clon = 139.8
clat = 36.09

clon = 139.75
clat = 36.09


CRS = 'ZONAL'
CRS = 'MERID'


nvar_l = [
          "qh",
          "w",
          "qv",
          "qh",
          "w",
          "qv",
         ]

stime = datetime( 2019, 8, 24, 15, 25, 30 )
etime = datetime( 2019, 8, 24, 15, 30, 0 )

time_l = []
time = stime
while time <= etime:
      time_l.append( time )
      time += timedelta( seconds=30 ) 

hgt_l = [ 4000, 4000, 500, 4000, 4000, 4000  ]
main( INFO, time_l=time_l, hgt_l=hgt_l, clat=clat, nvar_l=nvar_l, CRS=CRS, clon=clon )

