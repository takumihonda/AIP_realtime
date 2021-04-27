import numpy as np
import sys
import os

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from tools_AIP import read_nc_topo

quick = True
quick = False

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

    lon2d_, lat2d_, topo2d_ = read_nc_topo( dom=1 )

    lon2d_1 = lon2d_[2:-2,2:-2]
    lat2d_1 = lat2d_[2:-2,2:-2]
    topo2d_1 = topo2d_[2:-2,2:-2]


#    lons = lon2d_1[0,0] #np.min( lon2d_1 ) - 3
#    lone = lon2d_1[-1,-1] #np.max( lon2d_1 ) + 3
#    lats = lat2d_1[0,0] #np.min( lat2d_1 ) - 3
#    late = lat2d_1[-1,-1] #np.max( lat2d_1 ) + 3

    lons = np.min( lon2d_1 ) - 3
    lone = np.max( lon2d_1 ) + 3
    lats = np.min( lat2d_1 ) - 3
    late = np.max( lat2d_1 ) + 3

    ll_lon = lons
    ur_lon = lone

    ll_lat = lats
    ur_lat = late

#    fig = plt.figure( figsize=(11, 6) )
#    fig.subplots_adjust( left=0.05, bottom=0.03, right=0.98, top=0.97,
#                         wspace=0.05, hspace=0.01)

    fig, ax = plt.subplots(1, 1, figsize=(11.0,8.0) )
    fig.subplots_adjust( left=0.07, bottom=0.05, right=0.98, top=0.95,
                         wspace=0.2, hspace=0.01 )

 
    lon_0 = 135.0 
    lat_0 = 35.0
    method = "lcc"
    method = "merc"
    res = 'i'
#    if dom == 2:
#       res = 'h'
#       method = "merc"
#       lon_r=139.609
#       lat_r=35.861
#       res = 'f'

    m = prep_map( ax, method=method, res=res,
                  ll_lon=ll_lon, ur_lon=ur_lon,
                  ll_lat=ll_lat, ur_lat=ur_lat,
                  lon_0=lon_0, lat_0=lat_0 )   

    fs = 12
    pdlon = 10
    pdlat = 10
#       pdlon = 3
#       pdlat = 3
#    else:
#       pdlon = 0.5
#       pdlat = 0.5
#       fs = 18
    lw = 1.5
    lc = 'k'
    m.drawparallels( np.arange(0,60,pdlat), labels=[1,0,0,0],
                     fontsize=fs, color=lc, linewidth=lw )
    m.drawmeridians( np.arange(60,180,pdlon), labels=[0,0,0,1],
                     fontsize=fs, color=lc, linewidth=lw )

    m.drawlsmask( land_color='gray', ocean_color='aqua', lakes=True, alpha=0.25 )

    x1, y1 = m( lon2d_1, lat2d_1 )

    bbox = {'facecolor':'w', 'alpha':0.95, 'pad':1,
            'edgecolor':'w' }

    res_l = [ 18, 6, 1.5, 500 ]
    unit_l = [ "km", "km", "km", "m"]

    fc = 'k'
    fs = 14
    for dom in range( 1, 5 ):
 
       res = res_l[dom-1]
       unit = unit_l[dom-1]

       lon2d_, lat2d_, topo2d_ = read_nc_topo( dom=dom )
       lon2d_ = lon2d_[2:-2,2:-2]
       lat2d_ = lat2d_[2:-2,2:-2]
       topo2d_ = topo2d_[2:-2,2:-2]

       xloc = np.max(lon2d_)
       yloc = np.max(lat2d_) + 1

       if dom == 4:
          fc = 'b'
          lc = 'b'
          yloc = np.min(lat2d_) - 3

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

    ofig = 'doms_1and2_basemap.png'

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
