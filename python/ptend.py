import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, prep_proj_multi, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads_all, read_nc_lonlat, dist

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm

from scipy.interpolate import griddata



quick = True
#quick = False

def main( INFO, itime_l=[], EXP_l=[] ):

    ps3d = np.zeros( ( INFO["TMAX"], INFO["gy"], INFO["gx"] ) )
    ptend_l = np.zeros( ( len( EXP_l ), INFO["TMAX"]-2 ) )
    

    for i, EXP_ in enumerate( EXP_l ):

       INFO["FCST_DIR"] = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/dafcst".format( EXP_ )
       for itime in itime_l: 
          print( "initial", itime )
          dir_in = os.path.join( "dat_ptend", EXP_, )
          os.makedirs( dir_in, exist_ok=True)
          ofile = os.path.join( dir_in, "ptend2_abs_{0:}.npz".format(  itime.strftime('%Y%m%d%H%M%S') ) )
   
          try:
             data = np.load( ofile )
             ptend_ = data["ptend"]
          except:
             print( "No npz file ", ofile)
   
             for tlev in range( INFO["TMAX"] ):
                prs3d_ = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar="p" )
                ps3d[tlev,:,:] = prs3d_[0,:,:]
      
             ptend_ = np.average( np.abs( np.diff( ps3d, axis=0, n=2 ), ), axis=(1,2) ) / ( INFO["DT"]**2 )
             np.savez( ofile, ptend=ptend_ )
   
          ptend_l[i,:] += ptend_

    ptend_l = ptend_l / len( itime_l )

    print( "data plot")
    for i in range( len( EXP_l) ):
       plt.plot( ptend_l[i,:] )
    plt.show()
    sys.exit()



    print( prs3d.shape )
    plt.contourf( prs3d[0,:,:] )
    plt.colorbar()
    plt.show()
    plt.contourf( prs3d[10,:,:] )
    plt.colorbar()
    plt.show()
    sys.exit()

    lon2d_4, lat2d_4, topo2d_4 = read_nc_topo( dom=4 )
    flon2d = INFO["lon2d"]
    flat2d = INFO["lat2d"]
    mz1d, _, _ = read_obs_grads_latlon()
    mzidx = np.argmin( np.abs( mz1d - hgt ) )

    mask, mlon2d, mlat2d = read_mask_full()
    mask2d = mask[mzidx,:,:]

    print( mask2d.shape, mlon2d.shape, mlat2d.shape)

    fig, (( ax1,ax2,ax3 ), (ax4,ax5,ax6) ) = plt.subplots( 2, 3, figsize=( 13, 9.0 ) )
    fig.subplots_adjust( left=0.03, bottom=0.03, right=0.97, top=0.97,
                         wspace=0.15, hspace=0.05)

    ax_l = [ ax1, ax2, ax3, ax4,  ]

    if quick:
       res = "l"
    else:
       res = "f"


    lons = lon2d_4[0,0]
    lone = lon2d_4[-1,-1]

    lats = lat2d_4[0,0]
    late = lat2d_4[-1,-1]

    method = "merc"
    lon_r = 139.609
    lat_r = 35.861
    contc = "palegreen"
    contc = "burlywood"
    oc = "w"
    lon_0 = lon_r
    lat_0 = lat_r
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
                       (flon2d, flat2d),
                       #method='cubic',
                       method='nearest',
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
    cmap_dbz.set_bad( color='gray', alpha=0.5 )


    i1d = np.arange( lon2d_4.shape[0] ) + 1.0
    j1d = np.arange( lon2d_4.shape[1] ) + 1.0

    i1d -= np.mean( i1d )
    j1d -= np.mean( j1d )

    # 0.5km mesh
    j2d, i2d = np.meshgrid( i1d*0.5, j1d*0.5 )

    dist2d = np.sqrt( np.square(i2d) + np.square(j2d) )
    dist2d_ = dist( lon_r, lat_r, lon2d_4, lat2d_4 ) * 0.001

    norm = BoundaryNorm( levs_dbz, ncolors=cmap_dbz.N, clip=False)

    x_r, y_r = m_l[0]( lon_r, lat_r )

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    pnum_l = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)" ]

    for i , ax in enumerate( ax_l ):
       itime = time_l[i]
       tlev = tlev_l[i]

       if i<= 2: 
          obs3d, _, _, _ = read_obs_grads( INFO, itime=itime )
          obs2d_ = griddata( ( olon2d.ravel(), olat2d.ravel() ), 
                             obs3d[ozidx,:,:].ravel(),
                             (flon2d, flat2d),
                             #method='cubic',
                             method='nearest',
                            )

          x2d = fx2d
          y2d = fy2d
          #var2d = np.where( imask2d < 1.0, obs3d[ozidx,:,: ], np.nan )
          print( imask2d.shape, dist2d.shape  )
          #var2d = np.where( ( imask2d < 1.0 ) or ( dist2d > 60.0e3 ), obs2d_, np.nan )
          var2d = np.where( ( imask2d < 1.0 ) , obs2d_, np.nan )
       else:
          fcst3d, _ = read_fcst_grads( INFO, itime=itime, tlev=tlev , FT0=True, )
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
                   fontsize=8, fmt='%.0f km', colors="k" )

#       # DEBUG
#       CONT = ax.contour( x2d_, y2d_, dist2d_, 
#                          levels=[20, 40, 60], zorder=1,
#                          colors='r', linewidths=1.5,
#                          linestyles='dashed',
#                          )
#
#       ax.clabel( CONT, CONT.levels, inline=True, #inline_spacing=1, 
#                   fontsize=8, fmt='%.0fkm', colors="r" )




       ax.plot( x_r, y_r, ms=8.0, marker='o', color='r',
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

          ax.text( 1.002, 1.51, "(dBZ)",
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
          ax.text( 0.5, 1.01, itime.strftime('%H:%M:%S UTC %m/%d') ,
                  va='bottom', 
                  ha='center',
                  transform=ax.transAxes,
                  color='k', fontsize=11, )

          if i == 2:
             ax.text( 0.9, 1.01, "Z={0:.0f} km".format( hgt/1000 ),
                     va='bottom', 
                     ha='left',
                     transform=ax.transAxes,
                     color='k', fontsize=10, )

#       if i == 0 or i == 3:
       if i <= 2:
          tit = "PAWR obs"
       else:
          tit = "Forecast (FT={0:.0f} min)".format( tlev*0.5 )
   
       ax.text( 0.5, 0.99, tit,
               va='top', 
               ha='center',
               transform=ax.transAxes,
               color='k', fontsize=12, 
               bbox=bbox )

    ofig = "6p_obs_fcst_" + itime.strftime('%m%d') + ".png"
    print(ofig)

    if not quick:
       opath = "png"
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       plt.show()

#    plt.contourf( dist2d )
#    plt.colorbar()
#    plt.show()
#    sys.exit()

############3

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"
EXP = "20201117/D4_500m_CTRL"

EXP_l = [ 
#          "20201117/D4_500m_CTRL",
#          "20201117/D4_500m_H1V1",
          "20201117/D4_500m_H8V8",
          "20201117/D4_500m_H4V1",
          "20201117/D4_500m_H8V1",
        ]

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/dafcst".format( EXP )

# data should be stored in EXP/[time0]/dafcst
#time0 = datetime( 2019, 8, 24, 15, 0, 0 )
#
#stime = datetime( 2019, 8, 24, 15, 0, 30 )
#etime = datetime( 2019, 8, 24, 16, 0, 0 )

time0 = datetime( 2019, 8, 19, 13, 0, 0 )
stime = datetime( 2019, 8, 19, 13, 0, 30 )
etime = datetime( 2019, 8, 19, 14, 0, 0 )

fcst_zmax = 43

obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )

TMAX = 61
DT = 30.0
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
         "DT": DT,
         "TMAX": TMAX,
       }




itime_l = []
time = stime
while time <= etime:
   itime_l.append( time )
   time += timedelta( seconds=30 )

main( INFO, itime_l=itime_l, EXP_l=EXP_l )

