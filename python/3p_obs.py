
import sys
import numpy as np
from datetime import datetime
from tools_AIP import read_obs_grads, prep_proj_multi, read_nc_topo, read_mask, read_obs_grads_latlon

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm

quick = True

def main( INFO, time_l=[], hgt=3000.0 ):

    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    mz1d, mlon2d, mlat2d = read_obs_grads_latlon()
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    mask = read_mask()
    mask2d = mask[mzidx,:,:]

    fig, (( ax1,ax2,ax3 )) = plt.subplots( 1, 3, figsize=( 13, 4.5 ) )
    fig.subplots_adjust( left=0.05, bottom=0.03, right=0.96, top=0.95,
                         wspace=0.15, hspace=0.02)

    ax_l = [ ax1, ax2, ax3, ]

    if quick:
       res = "l"
    else:
       res = "h"


    lons = lon2d_4[0,0]
    lone = lon2d_4[-1,-1]

    lats = lat2d_4[0,0]
    late = lat2d_4[-1,-1]

    lon_0 = None
    lat_0 = None
    method = "merc"
    lon_r = 139.609
    lat_r = 35.861
    res = 'l'
    contc = "palegreen"
    oc = "w"

    m_l = prep_proj_multi( method, ax_l, fs=7, res=res, lw=0.0, 
                           ll_lon=lons, ur_lon=lone, ll_lat=lats, ur_lat=late, 
                           pdlon=0.5, pdlat=0.5, blon=lon_r, blat=lat_0,
                           contc=contc, oc=oc )

    time = datetime( 2019, 8, 24, 15, 0, 30 )
    obs3d, olon2d, olat2d, oz1d = read_obs_grads( INFO, itime=time )
    ozidx = np.argmin( np.abs( oz1d - hgt ) )

    # for pcolormesh
    olon2d -= np.abs( olon2d[1,0] - olon2d[0,0] )
    olat2d -= np.abs( olat2d[0,1] - olat2d[0,0] )
    mlon2d -= np.abs( mlon2d[1,0] - mlon2d[0,0] )
    mlat2d -= np.abs( mlat2d[0,1] - mlat2d[0,0] )
    x2d, y2d = m_l[0]( olon2d, olat2d )
    mx2d, my2d = m_l[0]( mlon2d, mlat2d )

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

    for i , ax in enumerate( ax_l ):
       obs3d, _, _, _ = read_obs_grads( INFO, itime=time_l[i] )

       ax.pcolormesh( x2d, y2d, obs3d[ozidx,:,: ], 
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

#       ax.pcolormesh( mx2d, my2d, mask2d, ) 

    plt.show()

############3

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"
EXP = "D4_500m_CTRL"

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

INFO = { "TOP": TOP,
         "EXP": EXP,
         "time0": time0,
       }

time_l = [
          datetime( 2019, 8, 24, 15, 20),
          datetime( 2019, 8, 24, 15, 40),
          datetime( 2019, 8, 24, 16,  0),
         ]

hgt = 3000.0

main( INFO, time_l=time_l, hgt=hgt )

