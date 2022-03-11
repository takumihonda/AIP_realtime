import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

from tools_AIP import read_obs, read_mask, read_nc_lonlat, read_obs_grads_latlon, rain2dbz, read_obs_nc, calc_fss

PLOT = True
PLOT = False
quick = True
quick = False


def read_fcst_nc( fn='', fstime=datetime( 2019, 8, 24, 15, 10, 0 ),
             vtime=datetime( 2019, 8, 24, 15, 10, 0 ), 
             fdtsec=30, height=3000.0, ):

    ftsec_ = ( vtime - fstime ).total_seconds()

#    tlev = int( ftsec_ / fdtsec )
#    print( tlev )

    nc = Dataset( fn, "r", format="NETCDF4" )
    var4d = nc.variables['Reflectivity'][:]
    lon1d = nc.variables['Longitude'][:]
    lat1d = nc.variables['Latitude'][:]
    z1d = nc.variables['Height'][:]
    ftimes = nc.variables['time'][:]
    nc.close()

    tlev = np.argmin( np.abs( ftimes - ftsec_ ) )

    zlev = np.argmin( np.abs( z1d - height ) )

    return( var4d[tlev,:,:,zlev], lon1d, lat1d )


def get_ts_bs( fcst2d, obs2d, thrs=10.0 ):

    fo = np.where( (fcst2d >= thrs) & (obs2d >= thrs), 1, 0 ) 
    fx = np.where( (fcst2d >= thrs) & (obs2d <  thrs), 1, 0 ) 
    xo = np.where( (fcst2d <  thrs) & (obs2d >= thrs), 1, 0 ) 

    noobs_ = np.where( obs2d < thrs, 1, 0 ) 
    
    # Note
    # this function calculates the "false alarm ratio," not the "false alarm rate."
    # 10.1175/2009WAF2222300.1 

    return( np.sum(fo) / ( np.sum(fo)+np.sum(fx)+np.sum(xo) ), 
            ( np.sum(fo) + np.sum(fx) ) / ( np.sum(fo) + np.sum(xo) ),
            np.sum(fx) / ( np.sum(fo)+np.sum(fx) ) )

#############

def main( INFO, itime=datetime(2019,9,3,2,0), tlev=0, theight=3000, dbz_thrs_l=[], 
          lons=0.0, lone=0.0, lats=0.0, late=0.0, AC=False, odz=300.0,
          ng_l=[] ):

    tlev = int( tlev )
    ftime = itime + timedelta( seconds=int(tlev)*30 )

    from scipy.interpolate import griddata

    USE_NOWCAST_DATA = True
#    USE_NOWCAST_DATA = False

    if USE_NOWCAST_DATA:
       obs3d, ostat = read_obs( utime=ftime, AC=AC )
       if not ostat:
          return( None, None, None, None, False)
       obs3d = rain2dbz( obs3d )
   
       # undef values
       obs3d[ obs3d < 10.0 ] = 5.0
       ok = np.argmin( np.abs( INFO["obsz"] - theight ) )
   
#       iobs2d = griddata( (olon2d.ravel(), olat2d.ravel()), obs3d[ok,:,:].ravel(), 
#                          (lon2d, lat2d), 
#                          method='cubic',
#                         )
       iobs2d = obs3d[ok,:,:]
       ostat = True
    else:

       ofn = '{0:}/{1:}/pawr_grads/pawr_obs_{2:}.nc'.format( INFO["TOP"], INFO["EXP"],
                    ftime.strftime('%Y%m%d-%H%M%S') )
       odat, olon, olat, olev, rlon, rlat = read_obs_nc( fn=ofn, dz=odz )
       iobs2d = griddata( (olon.ravel(), olat.ravel()), odat.ravel(), 
                          (INFO["olon2d"], INFO["olat2d"]), 
                          method='cubic',
                         )
       ostat = True
   


    if not ostat:
       return( None, None, None, None, ostat)

    ffn = '{0:}/{1:}.nc'.format( INFO["FCST_DIR"],
                    itime.strftime('%Y%m%d-%H%M%S') )
    try:
       ifcst2d, flon1d, flat1d = read_fcst_nc( fn=ffn, fstime=itime, vtime=ftime )
       fstat = True
    except:
       fstat = False

    if not fstat:
       print( "No fcst" )
       return( None, None, None, fstat)



    ifcst2d = griddata( ( INFO["lon2d"].ravel(), INFO["lat2d"].ravel()), ifcst2d.ravel(),
                       (INFO["olon2d"], INFO["olat2d"]),
                       method='cubic',
                      )




    olon2d = INFO["olon2d"]
    olat2d = INFO["olat2d"]
    if PLOT:
       plot( olon2d, olat2d,
             ifcst2d, itime=itime, tit="Fcst", ftsec=30*tlev )
       plot( olon2d, olat2d,
             iobs2d, itime=itime, tit="Obs", ftsec=30*tlev )


    ts_l = []
    bs_l = []
    fa_l = []
    for i, dbz in enumerate( dbz_thrs_l ):
#        if REGION:
#           ifcst2d_ = ifcst2d[ ( olon2d >= lons ) & ( olon2d <= lone ) &
#                              ( olat2d >= lats ) & ( olat2d <= late ) ]
#           iobs2d_ = iobs2d[ ( olon2d >= lons ) & ( olon2d <= lone ) & 
#                            ( olat2d >= lats ) & ( olat2d <= late ) ]
#        else:
        ifcst2d_ = ifcst2d
        iobs2d_ = iobs2d
        ts_, bs_, fa_ = get_ts_bs( ifcst2d_, iobs2d_,
                             thrs=dbz )
        ts_l.append( ts_ )
        bs_l.append( bs_ )
        fa_l.append( fa_ )

    fss_l = np.zeros( ( len( ng_l), len( dbz_thrs_l) ) )
    fss_l[:,:] = np.nan

    for i, ng in enumerate( ng_l ):
       for j, dbz in enumerate( dbz_thrs_l ):
          fss_l[i,j] = calc_fss( ng=ng, thrs=dbz, fcst2d=ifcst2d, obs2d=iobs2d )

    return( ts_l, bs_l, fa_l, fss_l, True )

###########


def plot( lon2d, lat2d, data2d, itime=datetime(2019, 6, 10, 8,0), tit="Fcst", ftsec=0 ):

    print( 'Max: {0:.2f}, Min:{1:.2f}'.format( np.max( data2d ), np.min( data2d ) ) )
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
    cmap_dbz.set_over('k', alpha=1.0)
    cmap_dbz.set_under('w', alpha=1.0)

    cmap = cmap_dbz
    levs = levs_dbz 
    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)

    cmap.set_bad( color='gray', alpha=0.3 )

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


    #SHADE = ax1.contourf(x2d, y2d,
    SHADE = ax1.pcolormesh(x2d, y2d,
                         data2d[:,:],
                         #levels=levs, 
                         vmin=np.min(levs), 
                         vmax=np.max(levs), 
                         cmap=cmap, norm=norm,
                         #extend='both',
                         )

    pos = ax1.get_position()
    cb_h = pos.height 
    ax_cb = fig.add_axes( [pos.x1+0.005, pos.y0, 0.01, 0.7] )
    cb = plt.colorbar( SHADE, cax=ax_cb, extend='both',
                       ticks=levs[:],  )

    tit_ = tit
    ctime = vtime.strftime('Valid: %H:%M:%S UTC %m/%d/%Y')
    ofig = tit + "_" + vtime.strftime('v%H%M%S_%Y%m%d') 
    if tit == "Fcst":
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




theight = 3000.0
#theight = 3500.0
#theight = 6000.0


#etime = stime
stime = datetime( 2019, 8, 24, 15, 0, 30 )
#stime = datetime( 2019, 8, 24, 15, 24, 30 )
#stime = datetime( 2019, 8, 24, 15, 26, 0 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )
#etime = stime

stime = datetime( 2019, 8, 19, 13, 0, 30 )
etime = datetime( 2019, 8, 19, 14, 0, 0 )

EXP = "d4"
#EXP = "d4_500m_ac"
#EXP = "d4_500m_ac_vr"
#EXP = "d4_500m_ac_z"

# attenuation correction
AC = True
#AC = False

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )
time0 = datetime( 2019, 8, 19, 13, 0, 0 )

tmin = 0
tmax = 61 # max time dimension does not include FT=0

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test/{0:}/dafcst_nc".format( EXP )
TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test"

fcst_zmax = 43


data_path = "../../dat4figs_JAMES/info"
#os.makedirs( data_path, exist_ok=True )
fn_info = '{0:}/data.npz'.format( data_path, )


#obsz, olon2d, olat2d = read_obs_grads_latlon()
#lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )
obsz = np.load( fn_info )['obsz']
olon2d = np.load( fn_info )['olon2d']
olat2d = np.load( fn_info )['olat2d']
lon2d = np.load( fn_info )['lon2d']
lat2d = np.load( fn_info )['lat2d']
hgt3d = np.load( fn_info )['hgt3d']
cz = np.load( fn_info )['cz']
ohgt3d = np.load( fn_info )['ohgt3d']

INFO = {
         "EXP": EXP, 
         "TOP": TOP, 
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
fa_l = np.zeros( ( tlev_max, len( dbz_thrs_l) ) )

ts_l[:] = np.nan
bs_l[:] = np.nan
fa_l[:] = np.nan

ng_l = [ 0, 1, 2, 3, 4, 6, ]
fss_l = np.zeros( ( tlev_max, len( ng_l), len( dbz_thrs_l ) ) )
fss_l[:] = np.nan



if AC:
   odir = "ts_npz/{0:}_acObs".format( INFO["EXP"] )
else:
   odir = "ts_npz/{0:}_noacObs".format( INFO["EXP"] )
os.makedirs( odir, exist_ok=True)


lats = 36.05
late = 36.1

lons = 139.7
lone = 139.8

time = stime
while (time <= etime):
  print( "Initial time:", time )


  ftime_l = []
  for i, tlev in enumerate( tlevs ):
      print( "start ", tlev)
      ts_l_, bs_l_, fa_l_, fss_l_, stat = main( INFO, itime=time, tlev=tlev, theight=theight, 
                            dbz_thrs_l=dbz_thrs_l, lons=lons, lone=lone, lats=lats, late=late,
                            AC=AC, ng_l=ng_l )
      print( "end ", tlev, stat )

      if stat:
         ts_l[i,:] = ts_l_
         bs_l[i,:] = bs_l_
         fa_l[i,:] = fa_l_
         fss_l[i,:,:] = fss_l_
      else:
         ts_l[i,:] = np.nan
         bs_l[i,:] = np.nan
         fa_l[i,:] = np.nan
         fss_l[i,:,:] = np.nan

      print( ts_l_, i )
      ftime_l.append( time + timedelta(seconds=int( tlev*30 )) )
  
  for i, dbz in enumerate( dbz_thrs_l ):
      fn_ts = "20220216_ac_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
      print( "Store", fn_ts, "\n")
      np.savez( os.path.join(odir,fn_ts), ts=np.array(ts_l[:,i]), bs=np.array(bs_l[:,i]), 
             fa=np.array(fa_l[:,i]),
             fss=np.array(fss_l[:,:,i]),
             ng_l=np.array(ng_l),
             times=ftime_l )
  
  time += timedelta(seconds=30)

