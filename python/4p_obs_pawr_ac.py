from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import sys
import os

from tools_AIP import prep_proj_multi_cartopy, get_cfeature, setup_grids_cartopy, read_obs_nc, set_cmap_pawr_scat
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

quick = True
#quick = False


def main( INFO={}, vtime=datetime( 2019, 8, 24, 15, 10, 0), 
          height=3000, exp_l=[], clat=36.1, dll=0.05, clon=136.0,
          texp_l=[], odz=300, otyp='z' ):

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    zmax = 9.0

    # radar location
    lon_r = 139.609
    lat_r = 35.861

    if clon < 0.0:
       clon = lon_r

    odat_l = []
    olon_l = []
    olat_l = []
    olev_l = []

    for i in range( 2 ):

        for exp_ in exp_l:
            ofn = '{0:}/{1:}/pawr_grads/pawr_obs_{2:}.nc'.format( INFO["TOP"], exp_,
                         vtime.strftime('%Y%m%d-%H%M%S') )
            if i == 0:
               odat, olon, olat, olev, rlon, rlat = read_obs_nc( fn=ofn, dz=odz, otyp=otyp )
            elif i == 1:
               odat, olon, olat, olev, rlon, rlat = read_obs_nc( fn=ofn, dz=odz , VERT=True, clat=clat, clon=clon, dll=dll, otyp=otyp )
    
            odat_l.append( odat )
            olon_l.append( olon )
            olat_l.append( olat )
            olev_l.append( olev )

    lons = 138.90161
    lone = 140.31638
    lats = 35.285637
    late = 36.43222

    xticks = np.arange( 134.0, 142, 0.2 )
    yticks = np.arange( 30.0, 45, 0.2 )

    res = '10m'
    if quick:
       res = '50m'

    land = get_cfeature( typ='land', res=res )
    coast = get_cfeature( typ='coastline', res=res )

    # original data is lon/lat coordinate
    data_crs = ccrs.PlateCarree()

    fig = plt.figure( figsize=(10, 9) )
    fig.subplots_adjust( left=0.06, bottom=0.03, right=0.95, top=0.95,
                         wspace=0.1, hspace=0.2)

    projection = ccrs.Mercator( latitude_true_scale=lat_r, ) 
    ax_l = []
    xfig = 2
    yfig = 2
    for i in range( 1, 5 ):
       if i <= 2:
          ax_l.append( fig.add_subplot( yfig,xfig,i, projection=projection ) )
       else:
          ax_l.append( fig.add_subplot( yfig,xfig,i, projection=None ) )

#    ax_l = prep_proj_multi_cartopy( fig, xfig=2, yfig=1, proj='merc', 
#                         latitude_true_scale=lat_r )
# 
#    for i in range( 2 ):
#       ax_l.append( fig.add_subplot( 2,2,2+i, ) )


    unit, cmap, levels, norm, vmin, vmax, tvar = set_cmap_pawr_scat( otyp=otyp )


    for i, ax in enumerate( ax_l ):

        if i <= 1:

           ax.set_extent([ lons, lone, lats, late ] )
           ax.add_feature( coast, zorder=0 )
    
           setup_grids_cartopy( ax, xticks=xticks, yticks=yticks, 
                                fs=9, lw=0.0 )
    
           SHADE = ax.scatter( olon_l[i], olat_l[i], s=5.0, c=odat_l[i], 
                       norm=norm, 
                       cmap=cmap,
                       vmin=vmin,
                       vmax=vmax,
                       transform=data_crs ) 
           ax.plot( rlon, rlat, marker='o', color='k', markersize=10,
                    transform=data_crs )

           if clat > 0.0:
              ax.plot( [ lons, lone ], [ clat, clat ], transform=data_crs,
                       color='k', ls='dashed', lw=1.0 )
           elif clon > 0.0:
              ax.plot( [ clon, clon ], [ lats, late ], transform=data_crs,
                       color='k', ls='dashed', lw=1.0 )
 
           ax.text( 0.99, 0.01, r'Z={0:.1f}$\pm${1:.1f} km'.format( height*0.001, odz*0.001 ), 
                   va='bottom', 
                   ha='right',
                   bbox=bbox,
                   transform=ax.transAxes,
                   color='k', fontsize=10, )

           if i == 1:
   
              cvtime = vtime.strftime( '%H:%M:%S' )
              ax.text( 1.0, 1.01, 'Valid at {0:}'.format( cvtime ),
                      va='bottom', 
                      ha='right',
                      transform=ax.transAxes,
                      color='k', fontsize=10, )
   
              pos = ax.get_position()
              cb_width = 0.006
              cb_height = pos.height*0.9
              ax_cb = fig.add_axes( [ pos.x1+0.003, pos.y0, 
                                      cb_width, cb_height] )
              cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                                 ticks=levels[::1], extend='both' )
       
              ax.text( 1.02, 0.95, unit, 
                      va='top', 
                      ha='left',
                      transform=ax.transAxes,
                      color='k', fontsize=10, )
   
        else:
           ydata = olev_l[i]*0.001
           if clat > 0.0:
              ax.set_xlim( lons, lone )
              xdata = olon_l[i]
           elif clon > 0.0:
              ax.set_xlim( lats, late )
              xdata = olat_l[i]
           ax.set_ylim( 0, zmax )
#           odat_l[i] = np.where( odat_l[i] < -20, np.nan, odat_l[i])
           SHADE = ax.scatter( xdata, ydata, s=5.0, c=odat_l[i], 
                       norm=norm, 
                       cmap=cmap,
                       vmin=vmin,
                       vmax=vmax,
                             )
            
        if odat_l[i].shape[0] > 0:
           print( "Max:{0:.1f}, Min:{1:.1f}\n".format(
                  np.max( odat_l[i] ), np.min( odat_l[i] ), unit ) )
    


        ax.text( 0.5, 1.01, texp_l[i],
                va='bottom', 
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=12, )
 




    tit = 'MP-PAWR observations'
    fig.suptitle( tit, fontsize=14 )

    if clon > 0.0:
       ref_ = clon
    else:
       ref_ = clat

    ofig = "4p_obs_ac_h{0:}km_v{1:}_crs{2:.2f}.png".format( height*0.001, vtime.strftime('%Y%m%d%H%M%S'), ref_ )

#    fig.tight_layout()
#    fig.canvas.draw()
#    axpos1 = ax_l[0].get_position()
#    axpos2 = ax_l[2].get_position()
#    ax_l[2].set_position( [axpos2.x0, axpos2.y0, axpos1.width, axpos2.height] )


    print( ofig )
    if not quick:
       opath = "png/attenuation"
       os.makedirs( opath, exist_ok=True )
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
    else:
       plt.show()




INFO = { "TOP": "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test" }

vtime = datetime( 2019, 8, 24, 15, 30, 0 )
#vtime = datetime( 2019, 8, 24, 15, 40, 0 )
#vtime = datetime( 2019, 8, 24, 15, 50, 0 )
#vtime = datetime( 2019, 8, 24, 16,  0, 0 )

vtime = datetime( 2019, 8, 19, 13, 30, 0 )
#vtime = datetime( 2019, 8, 19, 13, 40, 0 )

exp_l = [ "d4", "d4_attenuation_corrected"]
texp_l = [ "CTRL", "TEST", "CTRL", "TEST" ]

height = 3000.0
height = 1000.0
odz = 300

clat = 36.15
clat = 36.2
clat = -1
clat =  35.96
clon = 139.37
dll = 0.01

otyp = 'cz'
otyp = 'z'
otyp = 'vr'

main( INFO=INFO, vtime=vtime, exp_l=exp_l,
      texp_l=texp_l, height=height, odz=odz, clat=clat, dll=dll, clon=clon,
      otyp=otyp )
