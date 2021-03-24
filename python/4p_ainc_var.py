import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat, dist, get_cfeature, setup_grids_cartopy, prep_proj_multi_cartopy, read_ga_grads_all, read_ga_grads_dbz

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

def read_qh_grads_all( INFO, itime=datetime( 2019,8,24,15,30,0 ), typ='g'):

    qh = read_ga_grads_all( INFO, itime=itime, nvar="qc", typ=typ ) 
    for nvar in [ "qr", "qi", "qs", "qg" ]:
        qh += read_ga_grads_all( INFO, itime=itime, nvar=nvar, typ=typ ) 

    return( qh )

def main( INFO, time_l=[], hgt=3000.0, clat=40.0, nvar_l=["w"], 
      CRS="ZONAL", clon=139.8 ):

    if CRS == "ZONAL":
       clons = 139.5 + 0.001
       clone = 140.1 - 0.001
    elif CRS == "MERID":
#       clats = 35.8 + 0.001
#       clate = 36.3 - 0.001
       clats = 35.85 + 0.001 + 0.1
       clate = 36.25 - 0.001

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

    mask_, mlon2d, mlat2d = read_mask_full()
    mask = mask_[:22,:,:]
    mask2d = mask[mzidx,:,:]

    fig = plt.figure( figsize=( 9.5, 8.5) )
    fig.subplots_adjust( left=0.05, bottom=0.03, right=0.93, top=0.93,
                         wspace=0.2, hspace=0.15 )
 
    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    ax_l = prep_proj_multi_cartopy( fig, xfig=2, yfig=2, proj='merc', 
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

    olen2 = olat2d.shape[1] // 2
    oyidx = np.argmin( np.abs( olat2d[:,olen2] - clat ) )
    oxidx = np.argmin( np.abs( olon2d[olen2,:] - clon ) )

    flen2 = flat2d.shape[1] // 2
    fyidx = np.argmin( np.abs( flat2d[:,flen2] - clat ) )
    fxidx = np.argmin( np.abs( flon2d[flen2,:] - clon ) )



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


    levs_w = np.arange( -1, 1.2, 0.2)
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


    lons = flon2d[0,0]
    lone = flon2d[-2,-2]

    lats = flat2d[0,0]
    late = flat2d[-2,-2]
 
#    lons = flon2d[0,0] + 0.5 
#    lone = flon2d[-2,-2] - 0.15

    print( "Cross at lon: {0:.3f}-{1:.3f}".format( lons, lone ) )
    print( "Cross at lat: {0:.3f}-{1:.3f}".format( clat, clat ) )

    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )

#    ymin = 0.0 + 0.01
#    ymax = 9.0 - 0.01

#    xmin = lons + 0.05
#    xmax = lone - 0.05
    if CRS == "ZONAL":
       xmin = clons
       xmax = clone
    elif CRS == "MERID":
       xmin = clats
       xmax = clate

 

    dfz1d = np.diff( fz1d ) * 0.5
    dfz1d = np.append( dfz1d, dfz1d[-1] )

    ylab = 'Height (km)'

    # for pcolor mesh
    xlen2 = flon2d.shape[0] // 2
    ylen2 = flon2d.shape[1] // 2
    xlen = flon2d.shape[0] 
    ylen = flon2d.shape[1] 
    x2d = flon2d - ( flon2d[xlen2+1,ylen2] - flon2d[xlen2,ylen2] )
    y2d = flat2d - ( flat2d[xlen2,ylen2+1] - flat2d[xlen2,ylen2] )

    for i, ax in enumerate( ax_l ):
 
       ax.set_extent([ lons, lone, lats, late ] )
       ax.add_feature( land, zorder=0 )
       ax.add_feature( coast, zorder=0 )

       setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                            fs=8, lw=0.0 )
     
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

       elif nvar == "qg":
          fac = 1.e3
          tvar = "QG"
          cmap = cmap_qg
          levs = levs_qg
          unit = unit_qg

       elif nvar == "qr":
          fac = 1.e3
          tvar = "QR"
          cmap = cmap_qr
          levs = levs_qr
          unit = unit_qr

       elif nvar == "qs":
          fac = 1.e3
          tvar = "QS"
          cmap = cmap_qs
          levs = levs_qs
          unit = unit_qs

       elif nvar == "qh":
          fac = 1.e3
          tvar = "QHYD"
          cmap = cmap_qg
          levs = levs_qg

       elif nvar == "t":
          fac = 1.0
          cmap = cmap_t
          levs = levs_t
          tvar = "T"
          unit = unit_t

       elif nvar == "p":
          fac = 1.0
          cmap = cmap_p
          levs = levs_p
          tvar = "P"
          unit = "Pa"

       elif nvar == "hdiv":
          fac = 1.e4
          cmap = cmap_w
          levs = levs_hdiv
          unit = unit_hdiv
          tvar = "HDIV"

       elif nvar == "u" or nvar == "v":
          fac = 1.0
          cmap = cmap_u
          levs = levs_u
          unit = unit_ms
          if nvar == "u":
             tvar = "U"
          elif nvar == "v":
             tvar = "V"

       elif nvar == "dbz":
          fac = 1.0
          cmap = cmap_w
          levs = levs_hdiv
          unit = unit_hdiv
          tvar = "DBZ"

       cmap.set_over('k', alpha=1.0 )
       cmap.set_under('gray', alpha=1.0 )

       norm = BoundaryNorm( levs, ncolors=cmap.N, clip=False )

 
       for j, vtime in enumerate( time_l ):
          if j == 0:
             if nvar == "qh":
                g3d = read_qh_grads_all( INFO, itime=vtime, typ='g')
                a3d = read_qh_grads_all( INFO, itime=vtime, typ='a')
             elif nvar == "dbz":
                g3d = read_ga_grads_dbz( INFO, itime=vtime, typ='g')
                a3d = read_ga_grads_dbz( INFO, itime=vtime, typ='a')
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
             elif nvar == "dbz":
                g3d += read_ga_grads_dbz( INFO, itime=vtime, typ='g' ) 
                a3d += read_ga_grads_dbz( INFO, itime=vtime, typ='a' ) 
             else:
                g3d += read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='g' ) 
                a3d += read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='a' ) 

       var3d = ( a3d - g3d ) / len( time_l )
      
       var2d = var3d[fzidx,:,: ]

       xlen = x2d.shape[0] 
       ylen = x2d.shape[1] 
       print( np.max( var2d ), np.min( var2d) )


       SHADE = ax.pcolormesh( x2d, y2d, var2d[:xlen-1,:ylen-1]*fac, 
                       cmap=cmap, vmin=np.min(levs),
                       vmax=np.max(levs),
                       norm=norm, 
                       transform=data_crs,
                       )
 

       ctitx_l = [ 0.03, 0.97 ]

       if CRS == "ZONAL":
          ax.xaxis.set_major_formatter( FormatStrFormatter( '%.1fE' ) )
          ctit_l = [ "A", "B"]
       elif CRS == "MERID":
          ax.xaxis.set_major_formatter( FormatStrFormatter( '%.1fN' ) )
          #ctit_l = [ "C", "D"]
          ctit_l = [ "A", "B"]

#       if i == 0: 
#          ax.text( -0.06, -0.05, ylab,
#                  va='center', 
#                  ha='right',
#                  rotation=90,
#                  transform=ax.transAxes,
#                  color='k', fontsize=12, )

#          ax.set_ylabel( ylab, fontsize=10 )

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

#       ax.plot( [ lons, lone ], [ clat, clat ], transform=data_crs, 
#                linewidth=5.0, color='r',
#                  )
  
       if i == 1 or i == 3:
          pos = ax.get_position()
          cb_width = 0.006
          cb_height = pos.height*0.9
          ax_cb = fig.add_axes( [ pos.x1+0.003, pos.y0, #-cb_height*0.5, 
                                  cb_width, cb_height] )
          cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                             ticks=levs[::1], extend='both' )
          cb.ax.tick_params( labelsize=10 )
   
          ax.text( 1.01, 0.99, unit,
                  va='top', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, )


       ax.text( 0.01, 0.99, pnum_l[i],
               va='top', 
               ha='left',
               transform=ax.transAxes,
               color='k', fontsize=10, 
               bbox=bbox )

       ax.text( 0.5, 0.9, tvar,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=14, 
               bbox=bbox )

       ax.plot( [ clon, clon ], [ clats, clate ], 
                 color='r', linewidth=1.0, linestyle='dotted',
                 transform=data_crs )


       for j in range( 2 ):
           ax.text( ctitx_l[j], 0.01, ctit_l[j],
                   va='bottom', 
                   ha='center',
                   transform=ax.transAxes,
                   color='r', fontsize=11, 
                   bbox=bbox )
    


       if i == 1: 
          trange = '{0:}-{1:}'.format( 
                      time_l[0].strftime('%H%M:%S'), 
                      time_l[-1].strftime('%H%M:%S UTC %m/%d/%Y'), 
                    )
          ax.text( 0.0, 1.01, trange,
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=12, 
                   )

       if i == 1:
          ax.text( 0.9, 1.01, "Z={0:.1f} km".format( hgt/1000 ),
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, )


#       tit = "Forecast (FT={0:.1f} min)".format( tlev*30/60 )
#   
    tit = 'Time-averaged analysis increment'
    fig.suptitle( tit, fontsize=14 ) 


    ofig = "4p_ainc_{0:}_z{1:.1f}km.png".format( CRS, hgt/1000, )
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



hgt = 500.0
hgt = 4000.0
hgt = 8000.0


clon = 139.8
clat = 36.09

clon = 139.75
clat = 36.09


CRS = 'ZONAL'
CRS = 'MERID'


nvar_l = [
          #"dbz",
          "qv",
          "qh",
          "w",
          "hdiv",
#          "t",
         ]

stime = datetime( 2019, 8, 24, 15, 25, 0 )
etime = datetime( 2019, 8, 24, 15, 30, 0 )

time_l = []
time = stime
while time <= etime:
      time_l.append( time )
      time += timedelta( seconds=30 ) 

main( INFO, time_l=time_l, hgt=hgt, clat=clat, nvar_l=nvar_l, CRS=CRS, clon=clon )

