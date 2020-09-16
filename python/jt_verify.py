import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

PLOT = True
PLOT = False
quick = False
#quick = True

FT0 = False
#FT0 = True

def read_CZ( Kobe=False ):
    if Kobe:
       fn = "/data_honda01/honda/SCALE-LETKF/AIP_SAFE/DEBUG256p/D4_500m20191206_np256/20190610080000/anal/0001/init_20190610-080000.000.pe000000.nc"
    else:
       fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/TEST_INPUT/D4/20190910130000/anal/mean/init_20190910-130000.000.pe000000.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" )
    cz = nc.variables['CZ'][HALO:-HALO]
    nc.close()

    return( cz )

def read_nc_lonlat( Kobe=False ) :
    if Kobe:
       fn = "/data_honda01/honda/SCALE-LETKF/AIP_SAFE/DEBUG256p/D4_500m20191206_np256/const/topo_sno_np00001/topo.pe000000.nc"
    else:
       fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/domains/topo.d4.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" ) 
    lon = nc.variables['lon'][HALO:-HALO,HALO:-HALO]
    lat = nc.variables['lat'][HALO:-HALO,HALO:-HALO]
    topo = nc.variables['TOPO'][HALO:-HALO,HALO:-HALO]
    nc.close()

    cz = read_CZ( Kobe=Kobe ) 

    hgt3d = np.zeros( ( cz.shape[0], topo.shape[0],topo.shape[1]) )
    for k in range( cz.shape[0] ):
        hgt3d[k,:,:] = topo[:,:] + cz[k]

    return( lon, lat, hgt3d, cz )

def read_obs( utime=datetime(2019,9,3,2,0,0), Kobe=False ):

    jtime = utime + timedelta(hours=9)

    OBS_EXIST = False
    for sec in range( -15, 16, 1 ):
        jtime2 = jtime + timedelta( seconds=sec )
        if Kobe:
           fn = os.path.join("/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/test_kobe/obs/500m",
                             jtime2.strftime('%Y/%m/%d/%H/%M/%S'),
                             "rain_cart_0002.nc")
        else:
           fn = os.path.join("/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m/",
                             jtime2.strftime('%Y/%m/%d/%H/%M/%S'),
                             "rain_cart_0002.nc")

        if os.path.isfile( fn ):
           OBS_EXIST = True
           break
    
    if not OBS_EXIST:
       print( "Not found OBS ", fn)
       sys.exit()

    try:
       nc = Dataset( fn, "r", format="NETCDF4" ) 
    except: 
       print("Failed to open")
       print( fn )
       sys.exit()
    
    obs = nc.variables['rain'][0,:,:,:]

    return( obs )

def read_fcst_grads( INFO, itime=datetime(2019,9,3,2,0,0), tlev=0 ):

    #fn = os.path.join( "/data15/honda/SCALE-LETKF/AIP/verify/d4_500m_small",
    fn = os.path.join( INFO["TOP"], INFO["EXP"],
                       itime.strftime('%Y%m%d%H0000'),  "dafcst/fcst_ref3d_" +
                       itime.strftime('%Y%m%d-%H%M%S.grd') )
    
    try:
       infile = open(fn)
    except: 
       print("Failed to open")
       print( fn )
       sys.exit()
    
    # tlev starts from "0" (not from "1")
    gx = 256
    gy = 256
    gz = 60
    rec3d = gx*gy*gz

    nv = 1
 
#    # grads file starts from FT > 0s
#    gtlev = tlev - 1
#    rec = rec3d * nv * gtlev

    if FT0:
      rec = rec3d * nv * tlev
    else:
      rec = rec3d * nv * ( tlev - 1 )

    infile.seek(rec*4)
    tmp3d = np.fromfile(infile, dtype=np.dtype('>f4'), count=rec3d)  # big endian   
    input3d = np.reshape( tmp3d, (gz,gy,gx) )

    return( input3d ) 

def read_obs_grads_latlon( kmax=1, Kobe=False ):

    # tlev starts from "0" (not from "1")
    gx = 161
    gy = 161

    count = gx*gy

    lonlat = {"lon":None, "lat":None}

    for tvar in ["lon", "lat"]:
        if Kobe:
           fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/test_kobe/obs/500m/" + tvar + ".bin"
        else:
           fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m/" + tvar + ".bin"

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
    oz = np.arange( 0.0, dz*kmax, dz )

    return( lonlat["lon"], lonlat["lat"], oz )

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

def main( INFO, itime=datetime(2019,9,3,2,0), tlev=0, theight=3000, thrs_dbz=15.0, Kobe=False ):

    dbz_min = 15.0

    lon2d, lat2d, hgt3d, cz =  read_nc_lonlat( Kobe=Kobe )
    ##print(lon2d.shape, lat2d.shape, hgt3d.shape)

    #cz = read_CZ( Kobe=Kobe ) 
    ##print( cz.shape )
  
    ftime = itime + timedelta(seconds=tlev*30)
    fcst3d = read_fcst_grads( INFO, itime=itime, tlev=tlev )
    fcst3d[ fcst3d < dbz_min ] = 0.0
    ##print( fcst3d.shape )

    # vertical interpolation for fcst3d
    ifcst2d = vint_fcst( hgt3d, fcst3d, theight=theight )

    obs3d = read_obs( utime=ftime, Kobe=Kobe )
    ##print( obs3d.shape, np.min(obs3d), np.max(obs3d), ftime )

    olon2d, olat2d, ocz = read_obs_grads_latlon( kmax=obs3d.shape[0], Kobe=Kobe )

    # horizontal interpolation
    k = np.argmin( np.abs(ocz - theight) )
    obs3d[ obs3d < dbz_min ] = 0.0
    from scipy.interpolate import griddata
    iobs2d = griddata( (olon2d.ravel(), olat2d.ravel()), obs3d[k,:,:].ravel(), 
                       (lon2d, lat2d), 
                       method='cubic',
                      )

#    print(ifcst2d.shape, lon2d.shape)
    halo = 30
    imin = halo
    imax = -halo
    jmin = halo
    jmax = -halo

    ts, bs = get_ts_bs( ifcst2d[jmin:jmax,imin:imax], iobs2d[jmin:jmax,imin:imax], 
                        thrs=thrs_dbz)
    print( "TS: {0:.4f} BS: {1:.4f} FT: {2:0=4} s".format(ts, bs, tlev*30) )

    if PLOT:
       print( lon2d[jmin:jmax,imin:imax].shape )
       plot( lon2d[jmin:jmax,imin:imax], lat2d[jmin:jmax,imin:imax], 
             ifcst2d[jmin:jmax,imin:imax], itime=itime, tit="Fcst", ftsec=30*tlev )
       plot( lon2d[jmin:jmax,imin:imax], lat2d[jmin:jmax,imin:imax], 
             iobs2d[jmin:jmax,imin:imax], itime=itime, tit="Obs", ftsec=30*tlev )

    return( ts, bs )

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

    y2d, x2d = np.meshgrid( x1d, y1d )


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
EXP = "D4_500m_TEST_DEFAULT"

INFO = { "TOP": TOP,
         "EXP": EXP,
       }

theight = 3000.0
#theight = 6000.0
thrs_dbz = 15.0


#time = datetime(2019,9,3,2,50,0)
time = datetime(2019, 6, 10, 8, 10, 0)
time = datetime(2019, 6, 10, 8, 20, 0)

#stime = datetime(2019, 6, 10, 8, 11, 0)
etime = datetime(2019, 6, 10, 8, 29, 30)


stime = datetime(2019, 8, 24, 19, 0, 30)
stime = datetime(2019, 8, 24, 19, 5, 30)
stime = datetime(2019, 8, 24, 19, 15, 30)

stime = datetime( 2019, 8, 24, 15, 0, 30 )
stime = datetime( 2019, 8, 24, 15, 30, 30 )
etime = stime

if FT0:
  tmin = 0
  tmax = 61 # max time dimension includes FT=0
else:
  tmin = 1
  tmax = 61 # max time dimension does not include FT=0

tskip = 1

Kobe = False

time = stime
while (time <= etime):
  print( "Initial time:", time )

  ts_l = []
  bs_l = []

  ftime_l = []
  fn_ts = "TS_" + "thrs{:.1f}dbz_".format(thrs_dbz) + \
          time.strftime('i%H%M%S_%Y%m%d.npz')
  for tlev in range( tmin, tmax, tskip ):
  #for tlev in range(1,2):
      ts, bs = main( INFO, itime=time, tlev=tlev, theight=theight, thrs_dbz=thrs_dbz, Kobe=Kobe )
      ts_l.append( ts )
      bs_l.append( bs )
      ftime_l.append( time + timedelta(seconds=tlev*30) )
  
  print(fn_ts)
  
  odir = "ts_npz"
  os.makedirs( odir, exist_ok=True)
  np.savez(os.path.join(odir,fn_ts), ts=np.array(ts_l), bs=np.array(bs_l), times=ftime_l )

  time += timedelta(seconds=30)

sys.exit()

