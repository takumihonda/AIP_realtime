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

def get_theight( tk1d, hgt1d, tk=273.0 ):
    idx0 = np.argmin( np.abs( tk1d - tk ) )
    idx1 = idx0 - 1

    rat = ( tk1d[idx0] - tk ) / ( tk1d[idx0] - tk1d[idx1] )

    hgt = hgt1d[idx0] * rat + hgt1d[idx1] * ( 1.0 - rat )

    return( hgt )

def main( INFO, itime=datetime( 2019, 8, 24, 15, 0, 0), tlev_l=[], 
       lons=120, lone=150, lats=0, late=90,
       ):

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    cz = INFO["cz"]


    qc3d = read_fcst_grads_all( INFO, itime=itime, tlev=0 , FT0=True, nvar='qc' ) 
    lon3d = np.resize( flon2d, qc3d.shape )
    lat3d = np.resize( flat2d, qc3d.shape )

    afac = np.where( ( lon3d >= lons ) & ( lon3d <= lone ) &
              ( lat3d >= lats ) & ( lat3d <= late ), 1.0, np.nan )

    fig = plt.figure( figsize=(13, 8.5) )
    fig.subplots_adjust( left=0.04, bottom=0.1, right=0.96, top=0.97,
                         wspace=0.15, hspace=0.15 )

    yfig = 1
    xfig = 3
    ax_l = []
    for i in range( 1, xfig*yfig+1 ):
       ax_l.append( fig.add_subplot( yfig,xfig,i, ) ) #projection=projection ) )

    cu = 'k'

    labu = 'Wmax'

    lw = 2.0

    ymin = 0
    ymax = 13
    xmin = 0
    xmax = 20
    dx = 4    

    xticks = np.arange( xmin+dx, xmax+dx, dx )

    xlab = r'Maximum updraft (m s$^{-1}$)'
    ylab = 'Height (km)'

    tk0 = 273
    tk_l = [ tk0, ] #tk0-10, tk0-20 ]

    fac = 1.0
    zfac = 1.e-3
    for i, ax in enumerate( ax_l ):
       tlev = tlev_l[i]

       w3d = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar='w' ) 

       tk3d = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar='t' ) 

       w1d = np.nanmax( w3d*afac, axis=(1,2) )
       #w1d = -np.nanmin( w3d*afac, axis=(1,2) )

       tk1d = np.nanmean( tk3d*afac, axis=(1,2) )

       for tk in tk_l:
           hgt = get_theight( tk1d, cz, tk=tk )
           ax.hlines( y=hgt*zfac, xmin=xmin, xmax=xmax, linewidths=1.5, 
                      colors='gray', linestyle='dashed' ) 
           ax.text( xmax, hgt*zfac-0.1, r'{0:}$^\circ$C'.format( tk - tk0 ),
                   va='top', 
                   ha='right',
                   transform=ax.transData, 
                   color='gray', fontsize=11, )

       ax.plot( w1d*fac, cz*zfac, color=cu, linewidth=lw, label=labu )

       ax.vlines( x=xticks, ymin=ymin, ymax=ymax, linewidth=0.2, color='gray', 
                  linestyle='dashed' )

       ax.set_ylim( ymin, ymax )
       ax.set_xlim( xmin, xmax )

       ax.set_xticks( xticks )

       if i == 1: 
          ax.set_xlabel( xlab, fontsize=11 )
       if i == 0: 
          ax.set_ylabel( ylab, fontsize=11 )
          ax.legend( loc='upper right', fontsize=12 )

       tit = 'FT={0:.0f} min'.format( int( tlev*30/60 ) )
       ax.text( 0.5, 0.99, tit,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=12,
               bbox=bbox )

       if i == 2:
          ax.text( 1.0, 1.01, itime.strftime('%Y%m%d'),
                  va='bottom', 
                  ha='right',
                  transform=ax.transAxes,
                  color='k', fontsize=11,
                 )


    plt.show()
    sys.exit()


    ofig = "6p_fcst_crs_{0:}_{1:}_cll{2:.3f}_{3:}_{4:}.png".format(  itime.strftime('%m%d%H%M%S'), CRS, cll, nvar1, nvar2 )
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
#EXP = "20201117/D4_500m_CTRL_NOVR"

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






itime = datetime( 2019, 8, 24, 15, 30, 0 )
tlev_l = [ 0, 10, 20, 
           30, 40, 50 ]

lons = 139.35
lone = 139.9 
lats = 35.9
late = 36.22 

#itime = datetime( 2019, 8, 19, 13, 30, 0 )
#INFO["time0"] = datetime( 2019, 8, 19, 13, 0 )
#lons = 139.2
#late = 139.7
#lats = 35.7
#late = 36.1




main( INFO, itime=itime, tlev_l=tlev_l, lons=lons, lone=lone, lats=lats, late=late )

