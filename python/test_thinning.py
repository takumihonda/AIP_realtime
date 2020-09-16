import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

#PLOT = True
PLOT = False
quick = False
#quick = True

#FT0 = False
FT0 = True

NETCDF = True
NETCDF = False

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

def read_obs_grads( INFO, itime=datetime(2019,9,10,9) ):
    fn = os.path.join( INFO["TOP"], INFO["EXP"],
                       itime.strftime('%Y%m%d%H0000'),  "pawr_grads/pawr_ref3d_" +
                       itime.strftime('%Y%m%d-%H%M%S.grd') )
    
    try:
       infile = open(fn)
    except: 
       print("Failed to open")
       print( fn )
       sys.exit()

    print( fn )
    
    gx = 241
    gy = 241
    gz = 22
    rec3d = gx*gy*gz

    nv = 2
    rec = 0

    infile.seek(rec*4)
    tmp3d = np.fromfile(infile, dtype=np.dtype('>f4'), count=rec3d)  # big endian   
    input3d = np.reshape( tmp3d, (gz,gy,gx) )

    lons = 138.94319715
    lats = 35.32193244
    dlon = 0.00554812
    dlat = 0.00449640
    lon1d = np.arange( lons, lons+dlon*gx, dlon )
    lat1d = np.arange( lats, lats+dlat*gx, dlat )
    z1d = np.arange( 0.0, gz*500.0, 500.0 )

    lon2d, lat2d = np.meshgrid( lon1d, lat1d )

    return( input3d, lon2d, lat2d, z1d  ) 

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

#############

def main( INFO, itime=datetime(2019,9,3,2,0), Kobe=False ):

#    obs3d = read_obs( utime=itime, Kobe=Kobe )
#    olon2d, olat2d, ocz = read_obs_grads_latlon( kmax=obs3d.shape[0], Kobe=Kobe )
    obs3d, olon2d, olat2d, ocz = read_obs_grads( INFO, itime=itime )
    print( obs3d.shape, np.nanmin(obs3d), np.nanmax(obs3d), itime )

    i1d = np.arange( 0, olon2d.shape[1])
    j1d = np.arange( 0, olon2d.shape[0])
    i2d, j2d = np.meshgrid( i1d, j1d )

    print( "ok" )

    zlev = 5 
    zlev = 6 

    import matplotlib.pyplot as plt

    fig, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=( 10, 3 ) )
    fig.subplots_adjust( wspace=0.3, )

    cmap = plt.cm.get_cmap("jet")
    levs = np.arange( 10, 55, 5 )
 
    var = obs3d[zlev,:,:]
    var[ var <= -999000000.0 ] = np.nan

    ii, jj = 60, 40
    ii, jj = 69, 50


    ng = 4
    ng2 = 1
    i2d_ = i2d - ii
    j2d_ = j2d - jj
#    var1_ = np.where( ( i2d_ % ng == 0 ) & ( j2d_ % ng == 0 ), 
#                      var, np.nan )
    var1_ = np.where( ( ( i2d_ % ng == 0 ) & ( j2d_ % ng == 0 ) ) |
                        ( np.abs(i2d_) <= ng2) & ( np.abs(j2d_) <= ng2) , 
                      var, np.nan )

    ng2 = 3
#    i2d_ = i2d - ii + 1
#    j2d_ = j2d - jj + 1
    var2_ = np.where( ( ( i2d_ % ng == 0 ) & ( j2d_ % ng == 0 ) ) |
                        ( np.abs(i2d_) <= ng2) & ( np.abs(j2d_) <= ng2) , 
                      var, np.nan )

    var_l = [ var, var1_, var2_ ]

#    var = np.where( var > 0.0, 1.0, 0.0 )

#    SHADE = g.contourf( olon2d, olat2d, var,
#                        extend='both',
#                        cmap=cmap, levels=levs )

    ax_l = [ ax1, ax2, ax3 ]

    xmin_ = 139.2
    xmax_ = 139.4
    ymin_ = 35.4
    ymax_ = 35.6

    xmin_ = 139.3
    xmax_ = 139.4
    ymin_ = 35.5
    ymax_ = 35.6

    for i in range( len(ax_l) ):
        ax = ax_l[i]

        ax.pcolormesh( olon2d, olat2d, var_l[i], 
                            cmap=cmap, 
                            vmax=np.max(levs), vmin=np.min(levs))
        ax.set_aspect('equal')
    
        ax.set_xlim( xmin_, xmax_ )
    
        ax.set_ylim( ymin_, ymax_ )
    
        print( olon2d[jj,ii], olat2d[jj,ii] )
        ax.plot( olon2d[jj,ii], olat2d[jj,ii], marker='x', color='r', ms=4.0)

#    cb = plt.colorbar( SHADE )
    plt.show()

    sys.exit()



###########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"

EXP = "D4_500m_TEST_DEFAULT_0604_MEAN"

theight = 3000.0
theight = 6000.0
thrs_dbz1 = 15.0
thrs_dbz2 = 30.0


stime = datetime( 2019, 8, 24, 15, 0, 30 )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

if FT0:
  tmin = 0
  tmax = 61 # max time dimension includes FT=0
else:
  tmin = 1
  tmax = 61 # max time dimension does not include FT=0

INFO = { "TOP": TOP,
         "EXP": EXP,
         "time0": time0,
       }


tskip = 1

Kobe = False


main( INFO, itime=stime, Kobe=Kobe )
sys.exit()
