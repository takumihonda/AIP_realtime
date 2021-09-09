import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from datetime import timedelta
import os
import sys

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize

from matplotlib.colors import BoundaryNorm
import matplotlib.colors as mcolors

from tools_AIP import read_nc_topo, get_GFS_grads, get_grads_JMA, draw_rec


quick = False
quick = True

res = 'i'
if quick:
  res = 'c'
  res = 'i'



lat2 = 40.0
blon = 138.0
blat = 36.3

def main( gtime=datetime( 2019, 8, 24, 12, 0 ), htime=datetime( 2019, 8, 24, 12, 0 ), 
          jtime=datetime(2019, 8, 24, 15, 0 ), FT=2 ):

    # get MSLP
    slp2d, glon2d, glat2d = get_GFS_grads( gtime, var="MSLETmsl",zdim=-1)

    lon2d, lat2d, _ = read_nc_topo( dom=1 )
    lon2d_2, lat2d_2, _ = read_nc_topo( dom=2 )
    lon2d_3, lat2d_3, _ = read_nc_topo( dom=3 )
    lon2d_4, lat2d_4, _ = read_nc_topo( dom=4 )


    rain2d, rlon2d, rlat2d = get_grads_JMA( jtime, FT=FT, ACUM=True )
    jtime2 = jtime + timedelta( hours=FT )


    fig, ((ax1,ax2)) = plt.subplots(1, 2, figsize=(10.0, 5.0) )
    fig.subplots_adjust( left=0.05, bottom=0.05, right=0.94, top=0.95, )

    ax_l = [ ax1, ax2 ]
 

    lons1, lone1 = lon2d[0,0], lon2d[-1,-1]
    lats1, late1 = lat2d[0,0], lat2d[-1,-1]

    lons1 = 124
    lone1 = 160
    lats1 = 20

    lons2, lone2 = lon2d_3[0,0], lon2d_3[-1,-1]
    lats2, late2 = lat2d_3[0,0], lat2d_3[-1,-1]

    pdlon1, pdlat1 = 5, 5
    pdlon2, pdlat2 = 0.5, 0.5

    lat2 = 40.0

    blon = 135.0
    blat = 35.0

    lons_l = [ lons1, lons2 ]
    lone_l = [ lone1, lone2 ]
    lats_l = [ lats1, lats2 ]
    late_l = [ late1, late2 ]

    pdlon_l = [ pdlon1, pdlon2 ]
    pdlat_l = [ pdlat1, pdlat2 ]

    lc = 'k'
    cc = 'k'

    lw = 0.1
    cont_alp = 0.3
    fs = 8
    zorder = 2
    contc = 'burlywood'

    clevs = np.arange( 800, 1200, 4 )
#    rlevs = np.array( [ 5, 10, 15, 20, 25, 30, 35, 40, 45, 50] )
    rlevs = np.array( [ 0.5, 1, 5, 10, 15, 20, 25, 30 ] )


    pro_l = [ 'merc', 'merc' ]

    res_l = [ 'i', 'f' ]
    if quick:
       res_l = [ 'c', 'c' ]

    m_l = []
    for i, ax in enumerate( ax_l ):
         m = Basemap( projection=pro_l[i],resolution = res_l[i],
                  llcrnrlon = lons_l[i], llcrnrlat = lats_l[i],
                  urcrnrlon = lone_l[i], urcrnrlat = late_l[i],
                  lat_0 = blat, lat_1 = lat2,
                  lat_2 = lat2, lon_0 = blon,
                  ax = ax )
         m_l.append(m)
  
         m.drawcoastlines( linewidth = 0.5, color=cc, zorder=zorder)
         m.fillcontinents( color=contc, lake_color='w', zorder=0, alpha=cont_alp)
         m.drawparallels( np.arange(0,70,pdlat_l[i]), labels=[1,0,0,0],
                  fontsize=fs, color=lc,linewidth=lw )
         m.drawmeridians( np.arange(0,180,pdlon_l[i]), labels=[0,0,0,1],
                  fontsize=fs, color=lc, linewidth=lw )


    x1, y1 = m_l[0]( glon2d, glat2d )
    CONT1 = ax1.contour( x1, y1, slp2d*0.01, levels=clevs, 
                 colors='k', linewidths=1.0 )
    ax1.clabel( CONT1, fontsize=6, fmt='%.0f' )

    lw = 3.0
    lc = 'k'
    draw_rec( ax1, m_l[0], lon2d_2, lat2d_2, lw=lw, lc=lc )
    draw_rec( ax1, m_l[0], lon2d_3, lat2d_3, lw=lw, lc=lc )
    draw_rec( ax2, m_l[1], lon2d_4, lat2d_4, lw=lw, lc=lc )

    x2, y2 = m_l[1]( rlon2d, rlat2d )
#    cmap = plt.cm.get_cmap("jet")
    rlevs= np.array([0.5, 1, 5, 10, 15, 20, 25, 30, ])
    cmap = mcolors.ListedColormap(['lightskyblue', 'dodgerblue', 'b',
                                    'yellow',
                                    'orange', 'red', 
                                    'purple'])
    cmap.set_under('w', alpha=0.0)
    cmap.set_over('gray', alpha=1.0)
    norm = BoundaryNorm( rlevs, ncolors=cmap.N, clip=False)

    SHADE2 = ax2.contourf( x2, y2, rain2d, levels=rlevs,
                       extend='max', cmap=cmap, norm=norm )


    x1d = np.arange( lon2d_3.shape[0]) + 1.0
    y1d = np.arange( lon2d_3.shape[1]) + 1.0
    
    x1d -= np.mean( x1d )
    y1d -= np.mean( y1d )

    # 1.5km mesh
    y2d, x2d = np.meshgrid( x1d*1.5, y1d*1.5 )

    dist2d = np.sqrt( np.square(x2d) + np.square(y2d) )

    print(x1d.shape, np.min(x1d), np.max(x1d) )
    x2_2, y2_2 = m_l[1]( lon2d_3, lat2d_3 )
    CONT = ax2.contour( x2_2, y2_2, dist2d, 
                        levels=[20, 40, 60], zorder=1,
                        colors='k', linewidths=0.5,
                        linestyles='dashed',
                       )
    ax2.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
                fontsize=8, fmt='%.0fkm', colors="k" )

    lon_r = 139.609
    lat_r = 35.861
    x_r, y_r = m( lon_r, lat_r )
    ax2.plot( x_r, y_r, ms=4.0, marker='o', color='r',
              markeredgecolor='w' )

    ptit_l = [ "GFS MSLP (hPa)".format( 13, htime.strftime('%m/%d %H%M UTC'),), 
               "JMA radar precipitation (mm) {0:}-{1:}".format( jtime.strftime('%m/%d %H%M'), jtime2.strftime('%H%M UTC')) ]
    
#    SHADE_l = [ SHADE1, SHADE2 ]
    pos = ax2.get_position()
    cb_width = 0.01
    cb_height = pos.height*1.0
    ax_cb = fig.add_axes( [pos.x1, pos.y0+0.01, cb_width, cb_height] )
    cb = plt.colorbar( SHADE2, cax=ax_cb, orientation = 'vertical', ) #ticks=levs[::2])
    cb.ax.tick_params( labelsize=8 )

    pnum_l = [ "(a)", "(b)" ]

    for i, ax in enumerate( ax_l ):
        ax2.text( 0.5, 1.01, ptit_l[i],
                va='bottom', 
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=11 )

        ax2.text( 0.0, 1.01, pnum_l[i],
                va='bottom', 
                ha='right',
                transform=ax.transAxes,
                color='k', fontsize=10 )


    ofig = "2p_overiew_" + gtime.strftime('%m%d') + ".png"
    print(ofig)

    if not quick:
       opath = "png"
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       plt.show()


    sys.exit()
#########


gtime = datetime(2019, 8, 24, 12, 0, 0 )
htime = datetime(2019, 8, 24, 12, 0, 0 )

jtime = datetime(2019, 8, 24, 15, 0, 0 )
FT = 1

gtime = datetime( 2019, 8, 19, 12, 0, 0 )
jtime = datetime( 2019, 8, 19, 12, 0, 0 )
FT = 2

htime = gtime
main( gtime=gtime, htime=htime, jtime=jtime, FT=FT )


