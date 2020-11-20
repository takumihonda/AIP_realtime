import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

from tools_AIP import read_obs, read_mask, read_nc_lonlat, read_fcst_grads, read_obs_grads_latlon

PLOT = True
quick = True
PLOT = False
quick = False


#def read_fcst_nc( INFO, itime=datetime(2019,9,3,2,0,0), tlev=0 ):
#
#    fn = os.path.join( INFO["TOP"], INFO["EXP"],
#                       INFO["time0"].strftime('%Y%m%d%H0000'),  "dafcst/fcst_ref3d_" +
#                       itime.strftime('%Y%m%d-%H%M%S.nc') )
#
#    print( "NetCDF:", fn )
#
#    try:
#       nc = Dataset( fn, "r", format="NETCDF4" )
#    except: 
#       print("Failed to open")
#       print( fn )
#       sys.exit()
#
#    ref3d = nc.variables['Reflectivity'][tlev,:,:,:]
#    cz = nc.variables['Height'][:]
#    nc.close()
#
#    return( ref3d.transpose((2,1,0)), cz )

def get_ts_bs( fcst2d, obs2d, thrs=10.0 ):

    fo = np.where( (fcst2d >= thrs) & (obs2d >= thrs), 1, 0 ) 
    fx = np.where( (fcst2d >= thrs) & (obs2d <  thrs), 1, 0 ) 
    xo = np.where( (fcst2d <  thrs) & (obs2d >= thrs), 1, 0 ) 

    return( np.sum(fo) / ( np.sum(fo)+np.sum(fx)+np.sum(xo) ), 
            ( np.sum(fo) + np.sum(fx) ) / ( np.sum(fo) + np.sum(xo) ) )

def vint_fcst( hgt3d, fcst3d, theight=3000.0 ):

    # get ratio for averaging
    tmp3d = np.where(hgt3d - theight >= 0.0, 0.0, hgt3d )
    hgt_below = np.where(tmp3d == np.max(tmp3d, axis=0), tmp3d, 0.0)

    tmp3d = np.where(hgt3d - theight < 0.0, 999999.0, hgt3d)
    hgt_above = np.where(tmp3d == np.min(tmp3d, axis=0), tmp3d, 0.0)

    rat2d = ( np.max(hgt_above, axis=0) - theight) / \
            ( np.max(hgt_above, axis=0) - np.max(hgt_below, axis=0) )

    fcst_below = np.max( np.where( hgt_below > 0.0, fcst3d, 0.0 ), axis=0 )
    fcst_above = np.max( np.where( hgt_above > 0.0, fcst3d, 0.0 ), axis=0 )

    return( fcst_below*rat2d + fcst_above*(np.ones(rat2d.shape)-rat2d) )

#############

def main( INFO, itime=datetime(2019,9,3,2,0), tlev=0, theight=3000, dbz_thrs_l=[], ):

    tlev = int( tlev )
    ftime = itime + timedelta( seconds=int(tlev)*30 )

    obs3d, ostat = read_obs( utime=ftime, mask=INFO["mask"] )

    if not ostat:
       return( None, None, ostat)

    fcst3d, fstat,= read_fcst_grads( INFO, itime=itime, tlev=tlev, )

    if not fstat:
       print( "No fcst" )
       return( None, None, fstat)

    ifcst2d = vint_fcst( INFO["hgt3d"], fcst3d, theight=theight )


    from scipy.interpolate import griddata
    # undef values
    obs3d[ obs3d < 10.0 ] = 5.0
    ok = np.argmin( np.abs( INFO["obsz"] - theight ) )

#    iobs2d = griddata( (olon2d.ravel(), olat2d.ravel()), obs3d[ok,:,:].ravel(), 
#                       (lon2d, lat2d), 
#                       method='cubic',
#                      )

    iobs2d = obs3d[ok,:,:]
    ifcst2d = griddata( ( INFO["lon2d"].ravel(), INFO["lat2d"].ravel()), ifcst2d.ravel(),
                       (INFO["olon2d"], INFO["olat2d"]),
                       method='cubic',
                      )


    if PLOT:
       olon2d = INFO["olon2d"]
       olat2d = INFO["olat2d"]
       plot( olon2d, olat2d,
             ifcst2d, itime=itime, tit="Fcst", ftsec=30*tlev )
       plot( olon2d, olat2d,
             iobs2d, itime=itime, tit="Obs", ftsec=30*tlev )

    ts_l = []
    bs_l = []
    for i, dbz in enumerate( dbz_thrs_l ):
        ts_, bs_ = get_ts_bs( ifcst2d, iobs2d,
                             thrs=dbz )
        ts_l.append( ts_ )
        bs_l.append( bs_ )

    return( ts_l, bs_l, True )


def plot( lon2d, lat2d, data2d, itime=datetime(2019, 6, 10, 8,0), tit="Fcst", ftsec=0 ):

    vtime = itime + timedelta(seconds=ftsec)

    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    import matplotlib.colors as mcolors


    levs_dbz= np.array([15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65])
    cmap_dbz = mcolors.ListedColormap(['cyan', 'b', 'dodgerblue',
                                       'lime','yellow',
                                       'orange', 'red', 'firebrick', 'magenta',
                                       'purple'])
#    cmap_dbz = plt.cm.get_cmap("jet")
    cmap_dbz.set_over('gray', alpha=1.0)
    cmap_dbz.set_under('w', alpha=1.0)

    cmap = cmap_dbz
    levs = levs_dbz 
    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)

#    data2d[ data2d < np.min(levs)] = np.nan


    fig, (ax1) = plt.subplots(1, 1, figsize=(7.5,9) )
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, )

    ax1.set_aspect('equal')

#    from mpl_toolkits.basemap import Basemap
#
#    res = "i"
#    ll_lon = lon2d[0,0]
#    ur_lon = lon2d[-1,-1]
#    ll_lat = lat2d[0,0]
#    ur_lat = lat2d[-1,-1]
#    lon_0=139.609
#    lat_0=35.861
#
#    m = Basemap( projection="merc", resolution=res,
#                 llcrnrlon=ll_lon, llcrnrlat=ll_lat,
#                 urcrnrlon=ur_lon, urcrnrlat=ur_lat,
#                 lat_0=lat_0, lon_0=lon_0,
#                 ax=ax1)
#
#    m.drawcoastlines()
#
#    x2d, y2d = m( lon2d, lat2d )
#    SHADE = m.contourf( x2d, y2d,

    x1d = np.arange( lon2d.shape[0] ) * 0.5 # km
    y1d = np.arange( lon2d.shape[1] ) * 0.5 # km

    x1d -= np.mean(x1d)
    y1d -= np.mean(y1d)

    x2d, y2d = np.meshgrid( x1d, y1d )


    SHADE = ax1.contourf(x2d, y2d,
                         data2d[:,:],
                         levels=levs, 
                         #vmin=np.min(levs), 
                         #vmax=np.max(levs), 
                         cmap=cmap, norm=norm,
                         extend='both',
                         )

    pos = ax1.get_position()
    cb_h = pos.height 
    ax_cb = fig.add_axes( [pos.x1+0.005, pos.y0, 0.01, 0.7] )
    cb = plt.colorbar( SHADE, cax=ax_cb, extend='both',
                       ticks=levs[:],  )

    tit_ = tit
    ctime = vtime.strftime('Valid: %H:%M:%S UTC %m/%d/%Y')
    ofig = tit + "_" + vtime.strftime('v%H%M%S_%Y%m%d') 
    if tit is "Fcst":
       tit_ = tit 
       ctime = "FT=" + str(ftsec).zfill(4) + "s\n" + \
               vtime.strftime('Valid: %H:%M:%S UTC %m/%d/%Y') + "\n" + \
               itime.strftime('Init: %H:%M:%S UTC %m/%d/%Y')
       ofig = tit + "_" + itime.strftime('i%H%M%S_%Y%m%d_') + \
              vtime.strftime('v%H%M%S_%Y%m%d') 

    ax1.text( 0.5, 1.03, tit_,
              fontsize=15, transform=ax1.transAxes,
              horizontalalignment='center',
              verticalalignment='bottom' )

    ax1.text( 1.05, 1.01, ctime,
              fontsize=11, transform=ax1.transAxes,
              horizontalalignment='right',
              verticalalignment='bottom' )

    #ax1.set_xlabel( r'Longitude (deg)', fontsize=14 )
    #ax1.set_ylabel( r'Latitude (deg)', fontsize=14 )
    ax1.set_xlabel( r'X (km)', fontsize=14 )
    ax1.set_ylabel( r'Y (km)', fontsize=14 )


    ofig += "_z{:.1f}km.png".format(theight/1000) 
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

###########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"




EXP = "D4_500m_CTRL"

theight = 3000.0
#theight = 6000.0


#etime = stime
stime = datetime( 2019, 8, 24, 15, 0, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )
#stime = datetime( 2019, 8, 19, 13, 0, 30 )
#etime = datetime( 2019, 8, 19, 14, 0, 0 )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )
#time0 = datetime( 2019, 8, 19, 13, 0, 0 )

tmin = 0
tmax = 61 # max time dimension does not include FT=0

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/20201117/{0:}/dafcst".format( EXP )
#OBS_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/{1:}/pawr_grads".format( EXP )

fcst_zmax = 43

obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )
mask = read_mask()


INFO = {
         "mask": mask,
         "EXP": EXP, 
#         "OBS_DIR": OBS_DIR,
         "FCST_DIR": FCST_DIR,
         "obsz": obsz,
         "olon2d": olon2d,
         "olat2d": olat2d,
         "lon2d": lon2d,
         "lat2d": lat2d,
         "hgt3d": hgt3d,
         "ohgt3d": ohgt3d,
         "cz": cz,
         "gz": fcst_zmax,
         "gy": lon2d.shape[0],
         "gx": lon2d.shape[1],
       }



dbz_thrs_l = [ 15.0, 30.0, ]

tskip = 1

itmax = int( ( etime - stime ).total_seconds() / 30 ) + 1

tlevs = np.arange( tmin, tmax, tskip, dtype=np.int32 )
tlev_max = np.shape( tlevs )[0]

ts_l = np.zeros( ( tlev_max, len( dbz_thrs_l) ) )
bs_l = np.zeros( ( tlev_max, len( dbz_thrs_l) ) )

ts_l[:] = np.nan
bs_l[:] = np.nan

odir = "ts_npz/" + INFO["EXP"]
os.makedirs( odir, exist_ok=True)



time = stime
while (time <= etime):
  print( "Initial time:", time )


  ftime_l = []
  for i, tlev in enumerate( tlevs ):
      ts_l_, bs_l_, stat = main( INFO, itime=time, tlev=tlev, theight=theight, 
                            dbz_thrs_l=dbz_thrs_l )

      ts_l[i,:] = ts_l_
      bs_l[i,:] = bs_l_

      print( ts_l_, i )
      ftime_l.append( time + timedelta(seconds=int( tlev*30 )) )
  
  for i, dbz in enumerate( dbz_thrs_l ):
      fn_ts = "TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
      np.savez( os.path.join(odir,fn_ts), ts=np.array(ts_l[:,i]), bs=np.array(bs_l[:,i]), 
             times=ftime_l )
  
  time += timedelta(seconds=30)

