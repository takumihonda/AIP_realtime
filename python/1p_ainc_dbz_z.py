import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat, read_ga_grads_all, read_ga_grads_dbz

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

def main( INFO, time_l=[], nvar='dbz', 
          lons=120, lone=150, lats=0, late=90, ):


    fig, (ax1) = plt.subplots(1, 1, figsize=(5,7) )
    fig.subplots_adjust( left=0.15, bottom=0.05, right=0.93, top=0.93,
                         wspace=0.2, hspace=0.15 )
 
    ax_l = [ ax1 ]


 

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)" ]




    ylab = 'Height (km)'

    ymin = 0.0
    ymax = 12.0
    dy = 1.0
    yticks = np.arange( ymin, ymax+dy, dy )
   
    xmin = -0.2
    xmax = 0.2


    g3d = read_ga_grads_dbz( INFO, itime=time_l[0], typ='g' )
    lon3d = np.resize( INFO["lon2d"], g3d.shape )
    lat3d = np.resize( INFO["lat2d"], g3d.shape )
    afac = np.where( ( lon3d >= lons ) & ( lon3d <= lone ) &
              ( lat3d >= lats ) & ( lat3d <= late ), 1.0, np.nan )

    fac = 1.0
    if nvar == "qh":
       fac = 1.e3

    for i, ax in enumerate( ax_l ):
 
     
       unit = "dBZ"
       tvar = "DBZ"

       halo = 0

 
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
       xlen = var3d.shape[2]
       ylen = var3d.shape[1]

       var1d = np.nanmean( var3d[:,halo:ylen-halo,halo:xlen-halo]*
                           afac[:,halo:ylen-halo,halo:xlen-halo], axis=(1,2) )
       ax.plot( var1d*fac, INFO["cz"]*1.e-3, )

       ax.vlines( x=0.0, ymin=ymin, ymax=ymax, linewidth=1.0, color='gray', 
                  linestyle='dashed' )

       ax.set_ylim( ymin, ymax )
       ax.set_xlim( xmin, xmax )

       ax.set_ylabel( ylab, fontsize=12 )
 
       ax.set_yticks( yticks )

       ax.text( 0.01, 0.99, pnum_l[i],
               va='top', 
               ha='left',
               transform=ax.transAxes,
               color='k', fontsize=10, 
               bbox=bbox )

       tvar = nvar
       if nvar == "dbz":
          tvar = "Radar reflectivity (dBZ)"

       ax.text( 0.5, 0.9, tvar,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=13, 
               bbox=bbox )



       if i == 0: 
          trange = '{0:}-\n{1:}'.format( 
                      time_l[0].strftime('%H:%M:%S UTC'), 
                      time_l[-1].strftime('%H:%M:%S UTC\n%m/%d/%Y'), 
                    )
          ax.text( 0.0, 0.5, trange,
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, 
                   )

          arange = 'lon: {0:.2f}-{1:.2f}\nlat: {2:.2f}-{3:.2f}'.format( lons, lone, lats, late, 
                    )
          ax.text( 0.0, 0.4, arange,
                  va='bottom', 
                  ha='left',
                  transform=ax.transAxes,
                  color='k', fontsize=10, 
                   )



    tit = 'Area-averaged analysis increment'
    fig.suptitle( tit, fontsize=14 ) 


    ofig = "1p_ainc_zprof_{0:}_{1:}.png".format( nvar, time_l[0].strftime('%Y%m%d%H') )
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








stime = datetime( 2019, 8, 24, 15, 30, 0 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )
stime = datetime( 2019, 8, 24, 15, 20, 0 )
etime = datetime( 2019, 8, 24, 15, 30, 0 )
stime = datetime( 2019, 8, 24, 15, 30, 0 )
etime = datetime( 2019, 8, 24, 15, 40, 0 )
#etime = stime

#stime = datetime( 2019, 8, 19, 13, 30, 30 )
#etime = datetime( 2019, 8, 19, 14,  0,  0 )
#INFO["time0"] = datetime( 2019, 8, 19, 13, 0 )

time_l = []
time = stime
while time <= etime:
      time_l.append( time )
      time += timedelta( seconds=30 ) 


lons = lon2d[0,0] + 0.1
lone = lon2d[-1,-1] - 0.1

lats = lat2d[0,0] + 0.1
late = lat2d[-1,-1] - 0.1

lons = 139.35
lone = 139.9
lats = 35.9
late = 36.22

nvar = "w"
nvar = "dbz"
#nvar = "u"
#nvar = "v"
#nvar = "qv"
#nvar = "qh"
main( INFO, time_l=time_l, nvar=nvar, lons=lons, lone=lone, lats=lats, late=late )

