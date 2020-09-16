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

from tools_AIP import prep_proj_multi, read_nc, get_GFS_grads, get_grads_JMA, get_Him8_FLDK


quick = False
quick = True

res = 'i'
if quick:
  res = 'c'




lat2 = 40.0
blon = 138.0
blat = 36.3

def main( gtime=datetime( 2019, 8, 24, 12, 0 ), htime=datetime( 2019, 8, 24, 12, 0 ), 
          jtime=datetime(2019, 8, 24, 15, 0 ), FT=2 ):

    # get MSLP
    slp2d, glon2d, glat2d = get_GFS_grads( gtime, var="MSLETmsl",zdim=-1)

    lon2d, lat2d, _ = read_nc( dom=1 )
    lon2d3, lat2d3, _ = read_nc( dom=3 )


    rain2d, rlon2d, rlat2d = get_grads_JMA( jtime, FT=FT, ACUM=True )
    jtime2 = jtime + timedelta( hours=FT )

    tbb2d, hlon2d, hlat2d = get_Him8_FLDK( htime, band=13 )

    fig, ((ax1,ax2)) = plt.subplots(1, 2, figsize=(10.0, 5.0) )
    fig.subplots_adjust( left=0.1, bottom=0.1, right=0.9, top=0.9, )

    ax_l = [ ax1, ax2 ]
 
    res = "c"

    lons1, lone1 = lon2d[0,0], lon2d[-1,-1]
    lats1, late1 = lat2d[0,0], lat2d[-1,-1]
    lats1 = 26
    lone1 = 154
    lons1 = 122

    lons2, lone2 = lon2d3[0,0], lon2d3[-1,-1]
    lats2, late2 = lat2d3[0,0], lat2d3[-1,-1]

    pdlon1, pdlat1 = 5, 5
    pdlon2, pdlat2 = 1, 1

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
    cont_alp = 0.2
    fs = 10
    zorder = 2
    contc = 'burlywood'

    clevs = np.arange( 800, 1200, 4 )
    rlevs = np.array( [ 5, 10, 20, 30, 40, 50] )
#    hlevs = np.arange( 200, 300, 5 )
#
#    hcmap = mpl.colors.ListedColormap(['magenta','red', 'orange',
#                                      'gold','greenyellow', 'lime',
#                                      'turquoise','blue', 'cyan',
#                                      'white','whitesmoke','gainsboro',
#                                      'lightgray', 'silver','darkgray',
#                                      'gray','dimgray','black'])

    colors1 = plt.cm.jet_r(np.linspace(0, 1, 128))
    colors2 = plt.cm.binary(np.linspace(0., 1, 128)) # w/k
    colors = np.vstack((colors1, colors2))
    hcmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    hlevs = np.arange(200,304,4)

    pro_l = [ 'merc', 'merc' ]

    m_l = []
    for i, ax in enumerate( ax_l ):
         m = Basemap( projection=pro_l[i],resolution = res,
                  llcrnrlon = lons_l[i], llcrnrlat = lats_l[i],
                  urcrnrlon = lone_l[i], urcrnrlat = late_l[i],
                  lat_0 = blat, lat_1 = lat2,
                  lat_2 = lat2, lon_0 = blon,
                  ax = ax )
         m_l.append(m)
  
         m.drawcoastlines( linewidth = 0.5, color = cc, zorder=zorder)
         m.fillcontinents( color=contc,lake_color='w', zorder=0, alpha =cont_alp)
         m.drawparallels( np.arange(0,70,pdlat_l[i]), labels=[1,0,0,0],
                  fontsize=fs, color=lc,linewidth=lw )
         m.drawmeridians( np.arange(0,180,pdlon_l[i]), labels=[0,0,0,1],
                  fontsize=fs, color=lc, linewidth=lw )


    x1_, y1_ = m_l[0]( hlon2d, hlat2d )
    SHADE1 = ax1.contourf( x1_, y1_, tbb2d, levels=hlevs, 
                           cmap=hcmap, extend='max' )
    x1, y1 = m_l[0]( glon2d, glat2d )
    CONT1 = ax1.contour( x1, y1, slp2d*0.01, levels=clevs, 
                 colors='w', linewidths=1.0 )
    ax1.clabel( CONT1, fontsize=6, fmt='%.0f' )

    x2, y2 = m_l[1]( rlon2d, rlat2d )
    cmap = plt.cm.get_cmap("jet")
    SHADE2 = ax2.contourf( x2, y2, rain2d, levels=rlevs,
                       extend='max', cmap=cmap )

    
    ptit_l = [ "Himawari-8 B{0:0=2} (K) & GFS MSLP (hPa)".format( 13, htime.strftime('%m/%d %HUTC'),), 
               "JMA radar precipitation (mm) {0:}-{1:}".format( jtime.strftime('%m/%d %H'), jtime2.strftime('%HUTC')) ]
    
    SHADE_l = [ SHADE1, SHADE2 ]
    for i, ax in enumerate( ax_l ):
        pos = ax.get_position()
        cb_width = 0.01
        cb_height = pos.height*1.0
        ax_cb = fig.add_axes( [pos.x1, pos.y0+0.01, cb_width, cb_height] )
        cb = plt.colorbar( SHADE_l[i], cax=ax_cb, orientation = 'vertical', ) #ticks=levs[::2])
        cb.ax.tick_params( labelsize=8 )
    
        ax.text( 0.5, 1.01, ptit_l[i],
                verticalalignment='bottom', 
                horizontalalignment='center',
                transform=ax.transAxes,
                color='k', fontsize=12)

    plt.show()
    sys.exit()










  




    tit = r'Himawari-8 B13 (10.4 $\mu$m) Obs (K), GFS analysis MSLP (hPa), & JMA 48-h acum. precipitation (mm)'
    fig.suptitle(tit, fontsize=18)#, x = 0.08, y = 0.95, horizontalalignment = 'left')

    ofig = "5p_Him8_GFSANAL_JMA_" + itime.strftime('%m%d') + ".png"
    print(ofig)

    if not quick:
       opath = "png"
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       plt.show()


#########


gtime = datetime(2019, 8, 24, 12, 0, 0 )
#htime = datetime(2019, 8, 24, 12, 0, 0 )

jtime = datetime(2019, 8, 24, 14, 0, 0 )
FT = 2

htime = gtime
main( gtime=gtime, htime=htime, jtime=jtime, FT=FT )


