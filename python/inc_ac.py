import os
import sys
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import get_cfeature, setup_grids_cartopy, prep_proj_multi_cartopy, read_ga_grads_all, read_ga_grads_dbz, get_cz, read_mask_full

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

def read_nc( fn='', fstime=datetime( 2019, 8, 24, 15, 10, 0 ),
             vtime=datetime( 2019, 8, 24, 15, 10, 0 ), 
             fdtsec=30, height=3000.0, VERT=False, clat=36.1, clon=136.1 ):

    ftsec_ = ( vtime - fstime ).total_seconds()

#    tlev = int( ftsec_ / fdtsec )
#    print( tlev )

    nc = Dataset( fn, "r", format="NETCDF4" )
#    var4d = nc.variables['Reflectivity'][:]
    lon1d = nc.variables['Longitude'][:]
    lat1d = nc.variables['Latitude'][:]
    z1d = nc.variables['Height'][:]
    ftimes = nc.variables['time'][:]
    nc.close()

    tlev = np.argmin( np.abs( ftimes - ftsec_ ) )

    return( lon1d, lat1d, z1d )

def read_qh_grads_all( INFO, itime=datetime( 2019,8,24,15,30,0 ), typ='g'):

    qh = read_ga_grads_all( INFO, itime=itime, nvar="qc", typ=typ ) 
    for nvar in [ "qr", "qi", "qs", "qg" ]:
        qh += read_ga_grads_all( INFO, itime=itime, nvar=nvar, typ=typ ) 

    return( qh )

def set_cmap( nvar='qg' ):
    fac = 1.e3
    cmap = plt.cm.get_cmap("BrBG")
    levs = np.array( [ -2, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 2 ] )
    unit = r'g jg$^{{-1}}$'
    tvar = "TEST"
    if nvar == "dbz":
       fac = 1.e0
       unit = 'dBZ'
       #levs = np.array( [ -8, -6, -4, -2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 4, 6, 8, ] )
       levs = np.array( [ -4, -2, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 2, 4,] )
       tvar = "Radar reflectivity"
    elif nvar == "vr":
       fac = 1.e0
       unit = r'm s$^{{-1}}$'
       levs = np.array( [ -4, -2, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 2, 4,] )
       tvar = "Doppler velocity"
       cmap = plt.cm.get_cmap("RdBu_r")

    return( fac, cmap, levs, unit, tvar )

def main( INFO, time_l=[], hgt=3000.0, clat=40.0, nvar="w", exp_l=[],
      CRS="ZONAL", clon=139.8, gz=45, texp_l=[] ):

    ftime = datetime(2019, 8, 24, 15, 30 )
    ffn = '{0:}/{1:}/dafcst_nc/{2:}.nc'.format( INFO["TOP"], "d4",
              ftime.strftime('%Y%m%d-%H%M%S') )
    flon1d, flat1d, _ = read_nc( fn=ffn, fstime=ftime, vtime=ftime )
    flon2d, flat2d = np.meshgrid( flon1d, flat1d )

    fz1d = get_cz()

     
    zlev = np.argmin( np.abs( fz1d - hgt ) )

    mask, mlon2d, mlat2d = read_mask_full()
    mzmax = mask.shape[0]
    mz1d = np.arange( 0, 500*(mzmax+1), 500 )
    mzlev = np.argmin( np.abs( mz1d - hgt ) )
   
    mask2d = mask[mzlev,:,:]
    mask2d = np.where( mask2d < 0.5, np.nan, 1 ) 

    if CRS == "ZONAL":
       clons = 139.5 + 0.001
       clone = 140.1 - 0.001
    elif CRS == "MERID":
#       clats = 35.8 + 0.001
#       clate = 36.3 - 0.001
       clats = 35.85 + 0.001 + 0.1
       clate = 36.25 - 0.001

    # radar location
    lon_r = 139.609
    lat_r = 35.861


    fig = plt.figure( figsize=( 10.5, 5.5) )
    fig.subplots_adjust( left=0.05, bottom=0.03, right=0.93, top=0.93,
                         wspace=0.2, hspace=0.15 )
 
    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    ax_l = prep_proj_multi_cartopy( fig, xfig=2, yfig=1, proj='merc', 
                         latitude_true_scale=lat_r )


 
    res = '10m'
    if quick:
       res = '50m'

    land = get_cfeature( typ='land', res=res )
    coast = get_cfeature( typ='coastline', res=res )


    time = datetime( 2019, 8, 24, 15, 0, 30 )



    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)" ]



    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )

    if CRS == "ZONAL":
       xmin = clons
       xmax = clone
    elif CRS == "MERID":
       xmin = clats
       xmax = clate


    ylab = 'Height (km)'


    if nvar == "vr":
       VR = True
    else:
       VR = False


    for i, ax in enumerate( ax_l ):
 
       fac, cmap, levs, unit, tvar = set_cmap( nvar=nvar )
   
       cmap.set_over('k', alpha=1.0 )
       cmap.set_under('gray', alpha=1.0 )
       cmap.set_bad('r', alpha=1.0 )

#       ax.set_extent([ lons, lone, lats, late ] )
       ax.add_feature( land, zorder=0 )
       ax.add_feature( coast, zorder=0 )

       setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                            fs=8, lw=0.0 )
     

       norm = BoundaryNorm( levs, ncolors=cmap.N, clip=False )

       INFO["EXP"] = exp_l[i]
 

       for j, vtime in enumerate( time_l ):
          if nvar == "dbz" or nvar == "vr":
             if j == 0:
                g3d = read_ga_grads_dbz( INFO, itime=vtime, typ='g', VR=VR, gz=gz )
                a3d = read_ga_grads_dbz( INFO, itime=vtime, typ='a', VR=VR, gz=gz )
             else:
                g3d += read_ga_grads_dbz( INFO, itime=vtime, typ='g', VR=VR, gz=gz )
                a3d += read_ga_grads_dbz( INFO, itime=vtime, typ='a', VR=VR, gz=gz )
          else:
             if j == 0:
                g3d = read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='g', gz=gz )
                a3d = read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='a', gz=gz )
             else:
                g3d += read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='g', gz=gz )
                a3d += read_ga_grads_all( INFO, itime=vtime, nvar=nvar, typ='a', gz=gz )

       print( a3d.shape, g3d.shape )
       var3d = ( a3d - g3d ) / len( time_l )
      
       var2d = var3d[zlev,:,: ]*fac
       print( np.max( var2d ), np.min( var2d ) )

       SHADE = ax.contourf( flon2d, flat2d, var2d,
                            transform=data_crs, 
                            levels=levs,
                            norm=norm,
                            extend='both',
                            cmap=cmap )

       SHADE2 = ax.contourf( mlon2d, mlat2d, mask2d,
                            transform=data_crs, 
                            levels=[0.5, 1.5],
                            colors=['k', 'gray', 'b'],
                            alpha=0.3,
                            extend='both', )

       ax.plot( lon_r, lat_r, marker='o', color='k', markersize=10,
                transform=data_crs )

       ax.text( 0.5, 1.01, texp_l[i],
               va='bottom', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=12, )


       if i == 1:
          pos = ax.get_position()
          cb_width = 0.012
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

          ax.text( 0.98, 1.01, 'Z={0:.1f} km'.format( hgt*0.001 ),
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, )

          ax.text( 1.04, 1.04, '{0:}-{1:}'.format( 
                   time_l[0].strftime('%Y/%m/%d %H:%M'), 
                   time_l[-1].strftime('%H:%M') ),
                  va='bottom', 
                  ha='right',
                  transform=ax.transAxes,
                  color='k', fontsize=10, )



#   
    tit = 'Time-averaged analysis increment\n{0:}'.format( tvar )
    fig.suptitle( tit, fontsize=14 ) 


    ofig = "2p_ainc_{0:}_z{1:.1f}km_{2:}_{3:}.png".format( nvar, hgt/1000, 
           exp_l[0], exp_l[1] )
    print(ofig)

    if not quick:
       opath = "png/attenuation"
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       plt.show()





############3

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test"
EXP = "d4"
EXP ="d4_attenuation_corrected"

exp_l = [ "d4", "d4_attenuation_corrected" ]
texp_l = [ "CTRL", "TEST" ]
exp_l = [ "d4_500m_ac", "d4_500m_ac_z" ]
texp_l = [ "TEST", "TEST (Z only)" ]

exp_l = [ "d4_500m_ac_vr", "d4_500m_ac_vr" ]
texp_l = [ "TEST", "TEST (VR only)" ]

FCST_DIR = "{0:}/{1:}/dafcst".format( TOP, EXP )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

fcst_zmax = 43

#obsz, olon2d, olat2d = read_obs_grads_latlon()
#lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )

INFO = { "TOP": TOP,
         "time0": time0,
         "FCST_DIR": FCST_DIR,
#         "gz": fcst_zmax,
         "gy": 256, #lon2d.shape[0],
         "gx": 256, #lon2d.shape[1],
#         "lon2d": lon2d,
#         "lat2d": lat2d,
#         "cz": cz,
       }



hgt = 500.0
hgt = 2500.0


clon = 139.8
clat = 36.09

clon = 139.75
clat = 36.09


CRS = 'ZONAL'
CRS = 'MERID'


nvar = "dbz"
nvar = "vr"

stime = datetime( 2019, 8, 24, 15, 30, 0 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

stime = datetime( 2019, 8, 24, 15, 0, 30 )
etime = datetime( 2019, 8, 24, 15, 30, 0 )
stime = etime

gz = 45

time_l = []
time = stime
while time <= etime:
      time_l.append( time )
      time += timedelta( seconds=30 ) 

hgt_l = np.arange( 1, 11, 1 ) * 1000
hgt_l = [ 4000 ]
hgt_l = [ 3000 ]
#hgt_l = [ 2000 ]
hgt_l = [ 1000, 2000, 3000, 4000 ]
hgt_l = [ 1000 ]
for hgt in hgt_l:

    main( INFO, time_l=time_l, hgt=hgt, clat=clat, nvar=nvar, CRS=CRS, clon=clon, gz=gz, 
          exp_l=exp_l, texp_l=texp_l )

