import numpy as np

import matplotlib.pyplot as plt
from tools_AIP import prep_proj_multi_cartopy, setup_grids_cartopy, read_nc_topo

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

def main():

    lon2d_, lat2d_, topo2d_ = read_nc_topo( dom=1 )

    lon2d_1 = lon2d_[2:-2,2:-2]
    lat2d_1 = lat2d_[2:-2,2:-2]
    topo2d_1 = topo2d_[2:-2,2:-2]


    lons = np.min( lon2d_1 ) - 3
    lone = np.max( lon2d_1 ) + 3
    lats = np.min( lat2d_1 ) - 3
    late = np.max( lat2d_1 ) + 3

    fig = plt.figure( figsize=(11, 6) )
    fig.subplots_adjust( left=0.05, bottom=0.03, right=0.98, top=0.97,
                         wspace=0.05, hspace=0.01)
 
    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    xfig = 1
    yfig = 1
    ax_l = prep_proj_multi_cartopy( fig, xfig=xfig, yfig=yfig, proj='PlateCarree', )

    xticks = np.arange( 40, 205, 10 )
    yticks = np.arange( 0, 85, 10 )

    res = '50m'

    coast = cfeature.NaturalEarthFeature( 'physical', 'coastline', res,
                                         facecolor='none',
                                         edgecolor='k', )

    land = cfeature.NaturalEarthFeature( 'physical', 'land', res,
                                         edgecolor='face',
                                         facecolor=cfeature.COLORS['land'] )

    ocean = cfeature.NaturalEarthFeature( 'physical', 'ocean', res,
                                         edgecolor='face',
                                         facecolor=cfeature.COLORS['water'] )

    ax = ax_l[0]

    setup_grids_cartopy( ax, xticks=xticks, yticks=yticks,
                                fs=9, lw=0.25, color='k' )

    ax.set_extent([ lons, lone, lats, late ], crs=data_crs )

    ax.add_feature( coast, zorder=1 )

    ax.add_feature( land  )
    ax.add_feature( ocean )


    arrow_dict = dict(arrowstyle="->", color="k",
                  connectionstyle="angle3" )#"angle, angleA=0, angleB=10")

    bbox = {'facecolor':'w', 'alpha':0.95, 'pad':1,
            'edgecolor':'w' }

    lc = 'k'
    lw = 2.0

    res_l = [ 18, 6, 1.5, 500 ]
    unit_l = [ "km", "km", "km", "m"]

    fc = 'k'
    fs = 10
    transform = ccrs.PlateCarree()._as_mpl_transform(ax)
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


       ax.plot( lon2d_[0,:], lat2d_[0,:], color=lc, linewidth=lw,
                transform=data_crs )
       ax.plot( lon2d_[-1,:], lat2d_[-1,:], color=lc, linewidth=lw,
                transform=data_crs )

       ax.plot( lon2d_[:,0], lat2d_[:,0], color=lc, linewidth=lw,
                transform=data_crs )
       ax.plot( lon2d_[:,-1], lat2d_[:,-1], color=lc, linewidth=lw,
                transform=data_crs )


       xtloc = xloc + 1
       tit_ = 'Domain {0:} ({1:} {2:})'.format( dom, res, unit ) 

#       if dom >= 4:
#          xloc = np.min(lon2d_)
#          yloc = np.min(lat2d_)
#          if dom == 3:
#             ytloc = yloc + 1
#          else:
#             ytloc = yloc - 3
#             xtloc = xloc - 3
#
#          ax.annotate( tit_, xy=( xloc, yloc ), 
#                       xytext=( xtloc, ytloc ),
#                       size=fs, color=fc,
#                       xycoords=transform,
#                       arrowprops=arrow_dict,
##                       arrowprops=dict(facecolor='gray',
##                               arrowstyle="simple",
##                               #connectionstyle="arc3,rad=-0.2",
##                               alpha=0.5),
#                       )
#       else:
       xlen_ = int( lon2d_.shape[1] // 2 )
       xloc = lon2d_[-1,xlen_]
       ax.text( xloc, yloc, tit_, 
                fontsize=fs, color=fc,
                transform=data_crs,
                ha='center',
                va='bottom', 
                bbox=bbox )


#    ax2 = ax_l[1]
#
#    xticks = np.arange( 40, 205, 1 )
#    yticks = np.arange( 0, 85, 1 )
#
#    setup_grids_cartopy( ax2, xticks=xticks, yticks=yticks,
#                                fs=9, lw=0.25, color='k' )
#
#    ax2.add_feature( land  )
#    ax2.add_feature( ocean )
#
#    for dom in range( 3, 5 ):
#       lon2d_, lat2d_, topo2d_ = read_nc_topo( dom=dom )
#       lon2d_ = lon2d_[2:-2,2:-2]
#       lat2d_ = lat2d_[2:-2,2:-2]
#       topo2d_ = topo2d_[2:-2,2:-2]
#
#       if dom == 3:
#          lons = np.min( lon2d_ ) - 1
#          lone = np.max( lon2d_ ) + 1
#          lats = np.min( lat2d_ ) - 1
#          late = np.max( lat2d_ ) + 1
#
#          ax2.set_extent([ lons, lone, lats, late ], crs=data_crs )
#          fc = 'k'
#          fs = 10
#          lc = 'k'
#
#       elif dom == 4:
#          fc = 'b'
#          lc = 'b'
#
#       ax.plot( lon2d_[0,:], lat2d_[0,:], color=lc, linewidth=lw,
#                transform=data_crs )
#       ax.plot( lon2d_[-1,:], lat2d_[-1,:], color=lc, linewidth=lw,
#                transform=data_crs )
#
#       ax.plot( lon2d_[:,0], lat2d_[:,0], color=lc, linewidth=lw,
#                transform=data_crs )
#       ax.plot( lon2d_[:,-1], lat2d_[:,-1], color=lc, linewidth=lw,
#                transform=data_crs )

    plt.show()

###
main()
