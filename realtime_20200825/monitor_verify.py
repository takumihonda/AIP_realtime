import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

PLOT = True
PLOT = False
quick = False
#quick = True

def read_mask():
    fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_new/shadow/mask.nc"
    nc = Dataset( fn, "r", format="NETCDF4" )
    mask = nc.variables['mask'][0,:,:,:]

    nc.close()

    return( mask )

def read_obs( utime=datetime(2019,9,3,2,0,0) ):

    jtime = utime + timedelta(hours=9)

    OBS_EXIST = False
    for sec in range( 0, 30, 1 ):
        jtime2 = jtime + timedelta( seconds=sec )
        fn = os.path.join("/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_new/",
                          jtime2.strftime('%Y/%m/%d/%H/%M/%S'),
                          "rain_cart_0002.nc")

        if os.path.isfile( fn ):
           OBS_EXIST = True
           break

    if not OBS_EXIST:
       print( "Not found OBS ", fn)
       return( None, False )

    try:
       nc = Dataset( fn, "r", format="NETCDF4" )
    except:
       print("Failed to open")
       print( fn )
       sys.exit()

    obs = nc.variables['rain'][0,:,:,:]
    nc.close()

    obs = np.where( mask > 0.0, np.nan, obs )

    return( obs, True )

def read_obs_grads_latlon( ):

    # tlev starts from "0" (not from "1")
    gx = 161
    gy = 161
    gz = 29

    count = gx*gy

    lonlat = {"lon":None, "lat":None}

    for tvar in ["lon", "lat"]:
        fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_new/" + tvar + ".bin"

        try:
           infile = open(fn)
        except:
           print("Failed to open")
           print( fn )
           sys.exit()

        tmp2d = np.fromfile(infile, dtype=np.dtype('<f8'), count=count)  # big endian   
        input2d = np.reshape( tmp2d, (gy,gx) )

        lonlat[tvar] = input2d

    dz = 500.0
    oz = np.arange( 0.0, dz*gz, dz )

    return( oz, lonlat["lon"], lonlat["lat"] )

def get_obs_grid( ):

    gx = 241
    gy = 241

    lons = 138.94319715
    lats = 35.32193244
    dlon = 0.00554812
    dlat = 0.00449640

    olon1d = np.arange( lons, lons+dlon*gx, dlon )
    olat1d = np.arange( lats, lats+dlat*(gy-0.1), dlat )

    olon2d, olat2d = np.meshgrid( olon1d, olat1d )

    # obs contains 22 levels every 500 m
    gz = 22
    obs_dz = 500.0

    obsz = np.arange( 0, gz*obs_dz, obs_dz ) 

    return( obsz, olon2d, olat2d ) 


def read_CZ( ):
    fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/topo/init.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" )
    cz = nc.variables['CZ'][HALO:-HALO]
    nc.close()

    return( cz )

def read_nc_lonlat( fcst_zmax=43, obsz=np.arange(1)):

    #fn = "/work/jh200062/share/honda/SCALE-LETKF/monitor/topo/topo.d4.nc"
    fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/topo/topo.d4.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" ) 
    lon = nc.variables['lon'][HALO:-HALO,HALO:-HALO]
    lat = nc.variables['lat'][HALO:-HALO,HALO:-HALO]
    topo = nc.variables['TOPO'][HALO:-HALO,HALO:-HALO]
    nc.close()

    cz = read_CZ( ) 

    hgt3d = np.zeros( ( fcst_zmax, topo.shape[0],topo.shape[1]) )
    for k in range( fcst_zmax ):
        hgt3d[k,:,:] = topo[:,:] + cz[k]

    ohgt3d = np.zeros( ( obsz.shape[0], topo.shape[0],topo.shape[1]) )
    for k in range( obsz.shape[0] ):
        ohgt3d[k,:,:] = topo[:,:] + cz[k]

    return( lon, lat, hgt3d, cz, ohgt3d )

def read_obs_grads( INFO, itime=datetime(2019,9,10,9), ):

    fn = os.path.join( INFO["OBS_DIR"],
                       itime.strftime('pawr_ref3d_%Y%m%d-%H%M%S.grd') )

    print( fn )
    
    try:
       infile = open(fn)
    except: 
       print("Failed to open")
       print( fn )
       return( None, None, None, None, False )
    # obs does include Vr
    nv = 2
  
    gx = 241
    gy = 241


    # obs contains 22 levels every 500 m
    gz = 22
    rec3d = gx*gy*gz

    rec = 0

    infile.seek(rec*4)
    tmp3d = np.fromfile(infile, dtype=np.dtype('>f4'), count=rec3d)  # big endian   
    obs3d = np.reshape( tmp3d, (gz,gy,gx) )

    return( obs3d[:,:,:], True ) 

def read_fcst_grads( INFO, itime=datetime(2019,9,3,2,0,0), tlev=0 , FT0=True, ):

    fn = os.path.join( INFO["FCST_DIR"],
                       itime.strftime('%Y%m%d-%H%M%S.grd') )

    print( fn )
    try:
       infile = open(fn)
    except: 
       print("Failed to open")
       print( fn )
       return( None, False )
       sys.exit()
    
    # tlev starts from "0" (not from "1")
    #gz = 43 # only below 15km is stored in fcst
    gz = INFO["gz"]
    gx = INFO["gx"]
    gy = INFO["gy"]

    rec3d = gx*gy*gz

    nv = 1
 
    if FT0:
      rec = rec3d * nv * tlev
    else:
      rec = rec3d * nv * ( tlev - 1 )

    try:
       infile.seek(rec*4)
       tmp3d = np.fromfile(infile, dtype=np.dtype('>f4'), count=rec3d)  # big endian   
       input3d = np.reshape( tmp3d, (gz,gy,gx) )
       fstat = True
    except:
       input3d = None
       fstat = False

    return( input3d, fstat ) 

def get_ts_bs( fcst2d, obs2d, thrs=10.0 ):

    fo = np.where( (fcst2d >= thrs) & (obs2d >= thrs), 1, 0 ) 
    fx = np.where( (fcst2d >= thrs) & (obs2d <  thrs), 1, 0 ) 
    xo = np.where( (fcst2d <  thrs) & (obs2d >= thrs), 1, 0 ) 

    return( np.sum(fo) / ( np.sum(fo)+np.sum(fx)+np.sum(xo) ), 
            ( np.sum(fo) + np.sum(fx) ) / ( np.sum(fo) + np.sum(xo) ) )

def vint_fcst( hgt3d, fcst3d, theight=3000.0 ):

    # get ratio for averaging
    tmp3d = np.where( hgt3d - theight >= 0.0, 0.0, hgt3d )
    hgt_below = np.where( tmp3d == np.max(tmp3d, axis=0), tmp3d, 0.0 )

    tmp3d = np.where( hgt3d - theight < 0.0, 999999.0, hgt3d )
    hgt_above = np.where( tmp3d == np.min(tmp3d, axis=0), tmp3d, 0.0 )

    rat2d = ( np.max(hgt_above, axis=0) - theight) / \
            ( np.max(hgt_above, axis=0) - np.max(hgt_below, axis=0) )

    fcst_below = np.max( np.where( hgt_below > 0.0, fcst3d, 0.0 ), axis=0 )
    fcst_above = np.max( np.where( hgt_above > 0.0, fcst3d, 0.0 ), axis=0 )

    return( fcst_below*rat2d + fcst_above*(np.ones(rat2d.shape)-rat2d) )

#############

def main( INFO, itime=datetime(2019,9,3,2,0), tlev=0, theight=3000, dbz_thrs_l=[], ):

    tlev = int( tlev )
    ftime = itime + timedelta( seconds=int(tlev)*30 )

    #obs3d, ostat  = read_obs_grads( INFO, itime=ftime, )
    obs3d, ostat = read_obs( utime=ftime )

    if not ostat:
       return( None, None, ostat)

    fcst3d, fstat,= read_fcst_grads( INFO, itime=itime, tlev=tlev, )

    if not fstat:
       print( "No fcst" )
       return( None, None, fstat)

    ifcst2d = vint_fcst( hgt3d, fcst3d, theight=theight )
    #iobs2d = vint_fcst( ohgt3d, obs3d, theight=theight )

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


    # preview
#    levs = np.arange( 0, 50, 5)
#    import matplotlib.pyplot as plt
#    fig, (ax1) = plt.subplots(1, 1, figsize=(7.5,7) )
#    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, )
#    ax1.set_aspect('equal')
#    ax1.contourf( lon2d, lat2d, ifcst2d[:,:], levels=levs )
#    ax1.set_xlim(139, 140.2)
#    ax1.set_ylim(35.3, 36.3)
#    plt.show()
#    sys.exit()

#    halo = 0
#    imin = halo
#    imax = -halo
#    jmin = halo
#    jmax = -halo

    if PLOT:
       olon2d = INFO["olon2d"]
       olat2d = INFO["olat2d"]
       plot( olon2d, olat2d, 
             ifcst2d, itime=itime, tit="Fcst", ftsec=30*tlev )
       plot( olon2d, olat2d, 
             iobs2d, itime=itime, tit="Obs", ftsec=30*tlev )
       #plot( olon2d[jmin:jmax,imin:imax], olat2d[jmin:jmax,imin:imax], 
       #      ifcst2d[jmin:jmax,imin:imax], itime=itime, tit="Fcst", ftsec=30*tlev )
       #plot( olon2d[jmin:jmax,imin:imax], olat2d[jmin:jmax,imin:imax], 
       #      iobs2d[jmin:jmax,imin:imax], itime=itime, tit="Obs", ftsec=30*tlev )

    
    ts_l = []
    bs_l = []
    for i, dbz in enumerate( dbz_thrs_l ):
        #ts_, bs_ = get_ts_bs( ifcst2d[jmin:jmax,imin:imax], iobs2d[jmin:jmax,imin:imax], 
        ts_, bs_ = get_ts_bs( ifcst2d, iobs2d, 
                             thrs=dbz )
        ts_l.append( ts_ )
        bs_l.append( bs_ )

#    print( "TS: {0:.4f} BS: {1:.4f} TS: {2:.4f} BS: {3:.4f} FT: {4:0=4} s".format(ts1, bs1, ts2, bs2, tlev*30) )
    return( ts_l, bs_l, True )


def plot( lon2d, lat2d, data2d, itime=datetime(2019, 6, 10, 8,0), tit="Fcst", ftsec=0 ):

    vtime = itime + timedelta(seconds=ftsec)

    from matplotlib.colors import BoundaryNorm
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt


    levs_dbz= np.array([15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65])
    cmap_dbz = mcolors.ListedColormap(['cyan', 'b', 'dodgerblue',
                                       'lime','yellow',
                                       'orange', 'red', 'firebrick', 'magenta',
                                       'purple'])
#    cmap_dbz = plt.cm.get_cmap("jet")
    cmap_dbz.set_over('gray', alpha=1.0)
    cmap_dbz.set_under('k', alpha=1.0)

    cmap = cmap_dbz
    levs = levs_dbz 
    norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)

#    data2d[ data2d < np.min(levs)] = np.nan


    fig, (ax1) = plt.subplots(1, 1, figsize=(7.5,7.5) )
    fig.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.93, )

    ax1.set_aspect('equal')

    x1d = np.arange( lon2d.shape[0] ) * 0.5 # km
    y1d = np.arange( lon2d.shape[1] ) * 0.5 # km

    x1d -= np.mean(x1d)
    y1d -= np.mean(y1d)

    #y2d, x2d = np.meshgrid( y1d, x1d )
    y2d, x2d = lon2d, lat2d

    SHADE = ax1.contourf(y2d, x2d,
                         data2d[:,:],
                         levels=levs, 
                         #vmin=np.min(levs), 
                         #vmax=np.max(levs), 
                         cmap=cmap, norm=norm,
                         extend='both',
                         )

    pos = ax1.get_position()
    cb_h = pos.height 
    ax_cb = fig.add_axes( [pos.x1+0.005, pos.y0, 0.01, 0.5] )
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
       odir = "png/2d"
       os.makedirs( odir, exist_ok=True)
       plt.savefig( os.path.join(odir, ofig), 
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')

###########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"
# /work/jh200062/share/SCALE-LETKF-rt/result/ope/d4_500m/dafcst

OBS_DIR = os.path.join( TOP, "monitor_dev/obs" ) 
FCST_DIR = os.path.join( TOP, "monitor_dev/fcst" ) 

FCST_DIR = "/work/jh200062/share/SCALE-LETKF-rt/result/ope/d4_500m/dafcst"
OBS_DIR = "/work/jh200062/share/honda/SCALE-LETKF/monitor/obs"

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/fcst"
OBS_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/obs/pawr_obs"

theight = 3000.0
theight = 2000.0

fcst_zmax = 43

stime = datetime( 2020, 8, 27, 7, 11, 30 )
etime = datetime( 2020, 8, 27, 23, 0, 0 )
etime = datetime( 2020, 8, 28,  9, 0, 0 )

stime = datetime( 2020, 8, 27, 8, 0, 0 )
stime = datetime( 2020, 8, 27, 9, 0, 0 )
etime = stime

stime = datetime( 2020, 8, 29, 9, 0, 0 )
etime = datetime( 2020, 8, 31, 3, 0, 0 )

stime = datetime( 2020, 8, 31, 3, 0, 0 )
etime = datetime( 2020, 8, 31, 10, 0, 0 )

#stime = datetime( 2020, 9, 1, 0, 0, 0 )
#etime = datetime( 2020, 9, 2, 0, 0, 0 )

stime = datetime( 2020, 9, 5, 22, 42, 0 )
etime = datetime( 2020, 9, 6, 0, 0, 0 )

stime = datetime( 2020, 9, 1, 0, 0, 0 )
etime = datetime( 2020, 9, 2, 0, 0, 0 )
stime = datetime( 2020, 9, 2, 0, 0, 30 )
etime = datetime( 2020, 9, 3, 0, 0, 0 )

stime = datetime( 2020, 9, 3, 0, 0, 30 )
etime = datetime( 2020, 9, 4, 0, 0, 0 )

stime = datetime( 2020, 8, 31, 0, 0, 30 )
etime = datetime( 2020, 9, 1, 0, 0, 0 )

stime = datetime( 2020, 8, 30, 0, 0, 30 )
etime = datetime( 2020, 8, 31, 0, 0, 0 )

stime = datetime( 2020, 8, 29, 0, 0, 30 )
etime = datetime( 2020, 8, 30, 0, 0, 0 )

stime = datetime( 2020, 8, 28, 0, 0, 30 )
etime = datetime( 2020, 8, 29, 0, 0, 0 )

stime = datetime( 2020, 8, 27, 0, 0, 30 )
etime = datetime( 2020, 8, 28, 0, 0, 0 )

stime = datetime( 2020, 8, 26, 0, 0, 30 )
etime = datetime( 2020, 8, 27, 0, 0, 0 )

stime = datetime( 2020, 8, 25, 2, 26, 0 )

stime = datetime( 2020, 8, 24, 15, 36, 0 )
etime = datetime( 2020, 8, 25, 0, 0, 0 )
etime = datetime( 2020, 9, 1, 0, 0, 0 )

stime = datetime( 2020, 9, 1, 0, 0, 0 )
etime = datetime( 2020, 9, 7, 0, 0, 0 )

#etime = stime

tmin = 0
tmax = 61 # max time dimension includes FT=0

#obsz, olon2d, olat2d = get_obs_grid()
obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz )

mask = read_mask()

INFO = { 
         "mask": mask,
         "OBS_DIR": OBS_DIR,
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

tskip = 20

itmax = int( ( etime - stime ).total_seconds() / 30 ) + 1

tlevs = np.arange( tmin, tmax, tskip, dtype=np.int32 )
tlev_max = np.shape( tlevs )[0]

ts_l = np.zeros( ( tlev_max, len( dbz_thrs_l) ) )
bs_l = np.zeros( ( tlev_max, len( dbz_thrs_l) ) )

odir = "ts_npz/realtime_score"
os.makedirs( odir, exist_ok=True)

OVERW = False
OVERW = True

time = stime
while (time <= etime):
  time += timedelta(seconds=30)
  print( "Initial time:", time )

  ts_l[:] = np.nan
  bs_l[:] = np.nan

  exist = False
  for j, dbz in enumerate( dbz_thrs_l ):
     fn_ts = "TS_thrs{0:.1f}dbz_z{1:.1f}_tskip{2:0=3}_i{3:}.npz".format( dbz, theight, tskip, time.strftime('%H%M%S_%Y%m%d') )
     exist = os.path.exists( os.path.join(odir,fn_ts) )
     
  if exist and not OVERW:
     print( "Not overwrite" )
     continue

  ftime_l = []
  for i, tlev in enumerate( tlevs ):
      tlev = int( tlev )


      ts_l_, bs_l_, stat = main( INFO, itime=time, tlev=tlev, theight=theight, 
                               dbz_thrs_l=dbz_thrs_l, )
      
      if stat:
         print( "completed tlev={0:0=3}, itime:{1:}".format( tlev, time ) )
      else:
         print( "failed ")
         continue

      ts_l[i,:] = ts_l_
      bs_l[i,:] = bs_l_
      ftime_l.append( time + timedelta(seconds=int(tlev)*30) )
  
  for j, dbz in enumerate( dbz_thrs_l ):
     fn_ts = "TS_thrs{0:.1f}dbz_z{1:.1f}_tskip{2:0=3}_i{3:}.npz".format( dbz, theight, tskip, time.strftime('%H%M%S_%Y%m%d') )

     np.savez( os.path.join(odir,fn_ts), ts_l=ts_l[:,j], bs_l=bs_l[:,j], times=ftime_l )


sys.exit()

