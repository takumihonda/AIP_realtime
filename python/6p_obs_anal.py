import os
import sys
import numpy as np
from datetime import datetime
from tools_AIP import read_obs_grads, prep_proj_multi, read_nc_topo, read_mask, read_obs_grads_latlon, read_fcst_grads, read_nc_lonlat

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm

from scipy.interpolate import griddata

quick = True
#quick = False

def main( INFO, time_l=[], hgt=3000.0 ):

    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    mz1d, mlon2d, mlat2d = read_obs_grads_latlon()
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    mask = read_mask()
    mask2d = mask[mzidx,:,:]

    fig, (( ax1,ax2,ax3 ), (ax4,ax5,ax6) ) = plt.subplots( 2, 3, figsize=( 13, 9.0 ) )
    fig.subplots_adjust( left=0.03, bottom=0.03, right=0.97, top=0.97,
                         wspace=0.15, hspace=0.05)

    ax_l = [ ax1, ax2, ax3, ax4, ax5, ax6 ]

    if quick:
       res = "l"
    else:
       res = "f"


    lons = lon2d_4[0,0]
    lone = lon2d_4[-1,-1]

    lats = lat2d_4[0,0]
    late = lat2d_4[-1,-1]

    lon_0 = None
    lat_0 = None
    method = "merc"
    lon_r = 139.609
    lat_r = 35.861
    contc = "palegreen"
    contc = "burlywood"
    oc = "w"
    if quick:
       res = 'l'
    else:
       res = 'f'

    m_l = prep_proj_multi( method, ax_l, fs=7, res=res, lw=0.0, 
                           ll_lon=lons, ur_lon=lone, ll_lat=lats, ur_lat=late, 
                           pdlon=0.2, pdlat=0.2, blon=lon_r, blat=lat_0,
                           contc=contc, oc=oc )

    time = datetime( 2019, 8, 24, 15, 0, 30 )
    obs3d, olon2d, olat2d, oz1d = read_obs_grads( INFO, itime=time )
    ozidx = np.argmin( np.abs( oz1d - hgt ) )
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    imask2d = griddata( ( mlon2d.ravel(), mlat2d.ravel() ), mask[mzidx,:,:].ravel(),
                       (olon2d, olat2d),
                       method='cubic',
                      )

    # for pcolormesh
    olon2d -= np.abs( olon2d[1,0] - olon2d[0,0] )
    olat2d -= np.abs( olat2d[0,1] - olat2d[0,0] )
    mlon2d -= np.abs( mlon2d[1,0] - mlon2d[0,0] )
    mlat2d -= np.abs( mlat2d[0,1] - mlat2d[0,0] )


    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    fz1d = INFO["cz"]
    fzidx = np.argmin( np.abs( fz1d - hgt ) )

    ox2d, oy2d = m_l[0]( olon2d, olat2d )
    mx2d, my2d = m_l[0]( mlon2d, mlat2d )
    fx2d, fy2d = m_l[0]( flon2d, flat2d )

    x2d_, y2d_ = m_l[0]( lon2d_4, lat2d_4 )

    levs_dbz= np.array( [ 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 ] )
    cmap_dbz = mcolors.ListedColormap(['cyan', 'b', 'dodgerblue',
                                       'lime','yellow',
                                       'orange', 'red', 'firebrick', 'magenta',
                                       'purple'])
    cmap_dbz.set_over('gray', alpha=1.0)
    cmap_dbz.set_under('w', alpha=0.0)


    i1d = np.arange( lon2d_4.shape[0] ) + 1.0
    j1d = np.arange( lon2d_4.shape[1] ) + 1.0

    i1d -= np.mean( i1d )
    j1d -= np.mean( j1d )

    # 0.5km mesh
    j2d, i2d = np.meshgrid( i1d*0.5, j1d*0.5 )

    dist2d = np.sqrt( np.square(i2d) + np.square(j2d) )


    norm = BoundaryNorm( levs_dbz, ncolors=cmap_dbz.N, clip=False)

    x_r, y_r = m_l[0]( lon_r, lat_r )

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)" ]

    for i , ax in enumerate( ax_l ):
       itime = time_l[i]

       if i<= 2: 
          obs3d, _, _, _ = read_obs_grads( INFO, itime=itime )
          x2d = ox2d
          y2d = oy2d
          var2d = obs3d[ozidx,:,: ]
          #var2d = np.where( imask2d < 1.0, obs3d[ozidx,:,: ], 0.0 )
       else:
          fcst3d, _ = read_fcst_grads( INFO, itime=itime, tlev=0 , FT0=True, )
          x2d = fx2d
          y2d = fy2d
          var2d = fcst3d[fzidx,:,: ]

       SHADE = ax.pcolormesh( x2d, y2d, var2d, 
                       cmap=cmap_dbz, vmin=np.min(levs_dbz),
                       vmax=np.max(levs_dbz),
                       norm=norm, 
                       )

       CONT = ax.contour( x2d_, y2d_, dist2d, 
                          levels=[20, 40, 60], zorder=1,
                          colors='k', linewidths=0.5,
                          linestyles='dashed',
                          )
       ax.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
                   fontsize=8, fmt='%.0fkm', colors="k" )

       ax.plot( x_r, y_r, ms=4.0, marker='o', color='r',
                 markeredgecolor='w' )

       if i == 5:
          pos = ax.get_position()
          cb_width = 0.006
          cb_height = pos.height*1.0
          ax_cb = fig.add_axes( [ pos.x1+0.002, pos.y1-cb_height*0.5, 
                                  cb_width, cb_height] )
          cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                             ticks=levs_dbz[::1], extend='both' )
          cb.ax.tick_params( labelsize=8 )

          ax.text( 1.001, 1.51, "(dBZ)",
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

       if i <= 2:
          ax.text( 0.5, 1.01, itime.strftime('%H:%M:%S UTC') ,
                  va='bottom', 
                  ha='center',
                  transform=ax.transAxes,
                  color='k', fontsize=11, )

          if i == 2:
             ax.text( 0.9, 1.01, "Z={0:.0f}km".format( hgt/1000 ),
                     va='bottom', 
                     ha='left',
                     transform=ax.transAxes,
                     color='k', fontsize=10, )

#       if i == 0 or i == 3:
       if i <= 2:
          tit = "PAWR obs"
       else:
          tit = "Analysis"
   
       ax.text( 0.5, 0.99, tit,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=10, 
               bbox=bbox )

    ofig = "6p_obs_anal_" + itime.strftime('%m%d') + ".png"
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
EXP = "D4_500m_CTRL"

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/dafcst".format( EXP )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

fcst_zmax = 43

obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz )

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

time_l = [
          datetime( 2019, 8, 24, 15, 20),
          datetime( 2019, 8, 24, 15, 40),
          datetime( 2019, 8, 24, 16,  0),
          datetime( 2019, 8, 24, 15, 20),
          datetime( 2019, 8, 24, 15, 40),
          datetime( 2019, 8, 24, 16,  0),
         ]

hgt = 3000.0

main( INFO, time_l=time_l, hgt=hgt )

