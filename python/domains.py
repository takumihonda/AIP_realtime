import numpy as np
import os
import sys
from netCDF4 import Dataset

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm

from tools_AIP import read_nc_topo, draw_rec, dist

data_path = "../../dat4figs_JAMES/Fig01"
os.makedirs( data_path, exist_ok=True )

USE_ARCH_DAT = True
USE_ARCH_DAT = False

quick = False
#quick = True

def prep_map( ax, method='merc',lon_0=139.609, lat_0=35.861, 
              ll_lon=1, ur_lon=2,
              ll_lat=1, ur_lat=2, res='c',
              lw=0.5 ):

    lat1 = None
    lat2 = None

    m = Basemap( projection=method, resolution=res,
                 llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                 urcrnrlon=ur_lon, urcrnrlat=ur_lat,
                 lat_0=lat_0, lat_1=lat2,
                 lat_2=lat2, lon_0=lon_0,
                 ax=ax)

    m.drawcoastlines( linewidth=lw )
    #m.bluemarble()

    return( m )





def main( dom=1, bar=False ):
 
    fn = '{0:}/topo1_dom{1:}.npz'.format( data_path, dom )

    if USE_ARCH_DAT:
       topo2d_1 = np.load( fn )['topo']
       lat2d_1 = np.load( fn )['lat']
       lon2d_1 = np.load( fn )['lon']
    else:
       lon2d_1, lat2d_1, topo2d_1 = read_nc_topo( dom=dom )
       np.savez( fn, lon=lon2d_1, lat=lat2d_1, topo=topo2d_1 )

    print( lon2d_1.shape )

    fig, ax1 = plt.subplots(1, 1, figsize=(8.0,8.0))
    fig.subplots_adjust( left=0.07, bottom=0.05, right=0.94, top=0.95,
                         wspace=0.2, hspace=0.01 )

    ll_lon = lon2d_1[0,0]
    ur_lon = lon2d_1[-1,-1]

    ll_lat = lat2d_1[0,0]
    ur_lat = lat2d_1[-1,-1]
 
    if dom == 1 or dom == 2:
       lon_0 = 135.0 
       lat_0 = 35.0
       method = "lcc"
       res = 'i'
       if dom == 2:
          res = 'h'
          method = "merc"
    else:
       lon_0 = None
       lat_0 = None
       method = "merc"
       lon_r=139.609
       lat_r=35.861
       res = 'f'

    m = prep_map( ax1, method=method, res=res,
                  ll_lon=ll_lon, ur_lon=ur_lon,
                  ll_lat=ll_lat, ur_lat=ur_lat,
                  lon_0=lon_0, lat_0=lat_0 )   

    fs = 24
    if dom == 1:
       pdlon = 10
       pdlat = 10
    elif dom == 2:
       pdlon = 3
       pdlat = 3
    else:
       pdlon = 0.5
       pdlat = 0.5
       fs = 18
    lw = 0.3
    lc = 'k'
    m.drawparallels( np.arange(0,60,pdlat), labels=[1,0,0,0],
                     fontsize=fs, color=lc, linewidth=lw )
    m.drawmeridians( np.arange(60,180,pdlon), labels=[0,0,0,1],
                     fontsize=fs, color=lc, linewidth=lw )

    x1, y1 = m( lon2d_1, lat2d_1 )

    #levs = np.arange(0, 4200, 200)
    levs = np.array( [ -20, -15, -10, -5, 1, 5, 10, 50, 100, 200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400, 2800, 3200, 3600, ] )

    #norm = None
    #cmap = plt.cm.get_cmap("Greens")
    cmap = plt.cm.get_cmap("terrain")
    norm = BoundaryNorm( levs, ncolors=cmap.N, clip=False )

    SHADE1 = m.contourf( x1, y1, topo2d_1, 
                         levels=levs, zorder=1,
                         cmap=cmap, norm=norm,
                         extend='max',
                        )

    if dom <= 3:

       fn = '{0:}/dom{1:}.npz'.format( data_path, dom )
       if USE_ARCH_DAT:
          data = np.load( fn )['data']
          lon2d_2 = np.load( fn )['lon']
          lat2d_2 = np.load( fn )['lat']
       else:
          lon2d_2, lat2d_2, topo2d_2 = read_nc_topo( dom=dom+1 )
          data = topo2d_2
          np.savez( fn, data=topo2d_2, lon=lon2d_2, lat=lat2d_2 )
       x2, y2 = m( lon2d_2, lat2d_2 )
       #SHADE2 = m.contourf( x2, y2, topo2d_2, 
       SHADE2 = m.contourf( x2, y2, data, 
                            levels=levs, zorder=1,
                            cmap=cmap, norm=norm,
                            extend='max',
                          )

       lw = 2.0
       lc = 'k'
       if dom == 3:
          lc = 'r'
          x_r, y_r = m( lon_r, lat_r )
          ax1.plot( x_r, y_r, ms=8.0, marker='o', color='r',
                    markeredgecolor='w' )
       draw_rec( ax1, m, lon2d_2, lat2d_2, lw=lw, lc=lc )


    if dom == 3 or dom == 4:
       if dom == 3:
          x1d = np.arange( lon2d_2.shape[0]) + 1.0
          y1d = np.arange( lon2d_2.shape[1]) + 1.0
          xd = x2
          yd = y2
       elif dom == 4:
          x1d = np.arange( lon2d_1.shape[0]) + 1.0
          y1d = np.arange( lon2d_1.shape[1]) + 1.0
          xd = x1
          yd = y1

       x1d -= np.mean( x1d )
       y1d -= np.mean( y1d )

       # 0.5km mesh
       y2d, x2d = np.meshgrid( x1d*0.5, y1d*0.5 )

       # distance based on SCALE
       #dist2d = np.sqrt( np.square(x2d) + np.square(y2d) ) 
       # distance in projection (great cicle)
       dist2d = dist( lon_r, lat_r, lon2d_2, lat2d_2 ) * 0.001

       print(x1d.shape, np.min(x1d), np.max(x1d) )
       CONT = ax1.contour( xd, yd, dist2d, 
                           levels=[20, 40, 60], zorder=1,
                           colors='k', linewidths=0.5,
                           linestyles='dashed',
                          )
       ax1.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
                   fontsize=16, fmt='%.0fkm', colors="k" )


    if bar:
       plt.colorbar( SHADE1, ticks=levs[1:] )


    ofig = 'dom{0:}.png'.format( dom )

    print( ofig )
    if quick:
       plt.show()
    else:
       odir = "png"
       os.makedirs( odir, exist_ok=True)
       plt.savefig( os.path.join(odir, ofig),
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')
    

####

bar = True
bar = False

dom = 3
dom = 2
#dom = 1
#dom = 4
main( dom=dom, bar=bar )

