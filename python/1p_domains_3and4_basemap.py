import numpy as np
import sys
import os

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from tools_AIP import read_nc_topo, dist

quick = True
#quick = False

data_path = "../../dat4figs_GRL/FigS01-2"
os.makedirs( data_path, exist_ok=True )

USE_ARCH_DAT = False
#USE_ARCH_DAT = True


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



def main():

    fn = '{0:}/topo.npz'.format( data_path, )

    if USE_ARCH_DAT:
       lat2d_ = np.load( fn )['lat2d_']
       lon2d_ = np.load( fn )['lon2d_']
       topo2d_ = np.load( fn )['topo2d_']
       
    else:
       lon2d_, lat2d_, topo2d_ = read_nc_topo( dom=3 )
       np.savez( fn, lat2d_=lat2d_, lon2d_=lon2d_, topo2d_=topo2d_ )



    lon2d_1 = lon2d_[2:-2,2:-2]
    lat2d_1 = lat2d_[2:-2,2:-2]
    topo2d_1 = topo2d_[2:-2,2:-2]


#    lons = lon2d_1[0,0] #np.min( lon2d_1 ) - 3
#    lone = lon2d_1[-1,-1] #np.max( lon2d_1 ) + 3
#    lats = lat2d_1[0,0] #np.min( lat2d_1 ) - 3
#    late = lat2d_1[-1,-1] #np.max( lat2d_1 ) + 3

    dx = 0.5
    dy = 0.5
    lons = np.min( lon2d_1 ) - dx
    lone = np.max( lon2d_1 ) + dx
    lats = np.min( lat2d_1 ) - dy
    late = np.max( lat2d_1 ) + dy

    ll_lon = lons
    ur_lon = lone

    ll_lat = lats
    ur_lat = late

#    fig = plt.figure( figsize=(11, 6) )
#    fig.subplots_adjust( left=0.05, bottom=0.03, right=0.98, top=0.97,
#                         wspace=0.05, hspace=0.01)

    fig, ax = plt.subplots( 1, 1, figsize=(9.0,8.0) )
    fig.subplots_adjust( left=0.07, bottom=0.05, right=0.98, top=0.95,
                         wspace=0.2, hspace=0.01 )

 
    lon_0 = 135.0 
    lat_0 = 35.0
    method = "lcc"
    method = "merc"
    res = 'i'
    if not quick:
       res = 'f'
    lon_r=139.609
    lat_r=35.861

    m = prep_map( ax, method=method, res=res,
                  ll_lon=ll_lon, ur_lon=ur_lon,
                  ll_lat=ll_lat, ur_lat=ur_lat,
                  lon_0=lon_0, lat_0=lat_0 )   

    fs = 14
    pdlon = 1
    pdlat = 1
#       pdlon = 3
#       pdlat = 3
#    else:
#       pdlon = 0.5
#       pdlat = 0.5
#       fs = 18
    lw = 1.0
    lc = 'k'
    m.drawparallels( np.arange(0,60,pdlat), labels=[1,0,0,0],
                     fontsize=fs, color=lc, linewidth=lw )
    m.drawmeridians( np.arange(60,180,pdlon), labels=[0,0,0,1],
                     fontsize=fs, color=lc, linewidth=lw )

    m.drawlsmask( land_color='gray', ocean_color='aqua', lakes=True, alpha=0.25,
                  grid=1.25, resolution=res )

    x1, y1 = m( lon2d_1, lat2d_1 )

    bbox = {'facecolor':'w', 'alpha':0.95, 'pad':1,
            'edgecolor':'w' }

    res_l = [ 18, 6, 1.5, 500 ]
    unit_l = [ "km", "km", "km", "m"]

    lw = 1.5
    lcr = 'b'

    fc = 'k'
    fs = 14
    for dom in range( 3, 5 ):
 
       fn = '{0:}/topo_{1:0=2}.npz'.format( data_path, dom )



       res = res_l[dom-1]
       unit = unit_l[dom-1]

       if USE_ARCH_DAT:
          lat2d_ = np.load( fn )['lat2d_']
          lon2d_ = np.load( fn )['lon2d_']
          topo2d_ = np.load( fn )['topo2d_']
       
       else:
          lon2d_, lat2d_, topo2d_ = read_nc_topo( dom=dom )
          np.savez( fn, lat2d_=lat2d_, lon2d_=lon2d_, topo2d_=topo2d_ )

       lon2d_ = lon2d_[2:-2,2:-2]
       lat2d_ = lat2d_[2:-2,2:-2]
       topo2d_ = topo2d_[2:-2,2:-2]

       xloc = np.max(lon2d_)
       yloc = np.max(lat2d_) + 0.1

       if dom == 4:
          fc = 'b'
          lc = 'b'
#          yloc = np.min(lat2d_) 

       x2, y2 = m( lon2d_, lat2d_ )

       ax.plot( x2[0,:],  y2[0,:],  color=lc, linewidth=lw, )
       ax.plot( x2[-1,:], y2[-1,:], color=lc, linewidth=lw, )

       ax.plot( x2[:,0],  y2[:,0],  color=lc, linewidth=lw, )
       ax.plot( x2[:,-1], y2[:,-1], color=lc, linewidth=lw, )


#       xtloc = xloc + 0
       tit_ = 'Domain {0:} ({1:} {2:})'.format( dom, res, unit ) 
       xlen_ = int( lon2d_.shape[1] // 2 )
       xloc = lon2d_[-1,xlen_]

       xt, yt = m( xloc, yloc )

       ax.text( xt, yt, tit_, 
                fontsize=fs, color=fc,
                #transform=data_crs,
                ha='center',
                va='bottom', 
                transform=ax.transData,
                bbox=bbox )

       if dom == 4:
          xr, yr = m( lon_r, lat_r )
          xtr, ytr = m( lon_r+0.05, lat_r+0.05 )
          ax.plot( xr,  yr, marker='o', color=lcr, linewidth=lw, )
 
          ax.text( xtr, ytr, "MP-PAWR", 
                   fontsize=12, color=fc,
                   ha='left',
                   va='bottom', 
                   transform=ax.transData, )

  
          dist2d = dist( lon_r, lat_r, lon2d_, lat2d_ ) * 0.001
          CONT = ax.contour( x2, y2, dist2d, 
                              levels=[20, 40, 60], zorder=1,
                              colors='k', linewidths=0.5,
                              linestyles='dashed',
                             )
          ax.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
                      fontsize=10, fmt='%.0fkm', colors="k" )
   


    ofig = 'doms_3and4_basemap.png'

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


###
main()
