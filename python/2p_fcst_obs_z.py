from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
import sys

def read_obs_nc( fn='', height=3000.0, dz=100.0, otyp=4001 ):
    nc = Dataset( fn, "r", format="NETCDF4" )
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    lev = nc.variables['lev'][:]
    dat = nc.variables['dat'][:]
    elm = nc.variables['elm'][:]
    nc.close()

    print( dat.shape )
    #dat = dat[ ( elm == 4003 ) * ( lev == height ) ]
    dat = dat[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]
    lon = lon[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]
    lat = lat[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]

    lev = lev[ ( np.abs( lev - height ) < dz ) * ( elm == otyp ) ]
    #print( lev )

    return( dat, lon, lat )

def read_nc( fn='', fstime=datetime( 2019, 8, 24, 15, 10, 0 ),
             vtime=datetime( 2019, 8, 24, 15, 10, 0 ), 
             fdtsec=60 ):

    ftsec_ = ( vtime - fstime ).total_seconds()

    tlev = int( ftsec_ / fdtsec )
    print( tlev )

    nc = Dataset( fn, "r", format="NETCDF4" )
    var4d = nc.variables['Reflectivity'][:]
    lon1d = nc.variables['Longitude'][:]
    lat1d = nc.variables['Latitude'][:]
    z3d = nc.variables['Height'][:]
    ftimes = nc.variables['time'][:]
    nc.close()

    return( var4d[tlev,:,:,:], lon1d, lat1d, z3d )

def main( INFO={}, vtime=datetime( 2019, 8, 24, 15, 10, 0), 
          fstime=datetime( 2019, 8, 24, 15, 10, 0 ) ):


    ofn = '{0:}/pawr_grads/pawr_obs_{1:}.nc'.format( INFO["TOP"],
                 vtime.strftime('%Y%m%d-%H%M%S') )
    odat, olon, olat = read_obs_nc( fn=ofn )
    print( odat.shape )

    ffn = '{0:}/dafcst_nc/{1:}.nc'.format( INFO["TOP"],
                 fstime.strftime('%Y%m%d-%H%M%S') )
    print( ffn )
    fdat3d, flon1d, flat1d, fz3d = read_nc( fn=ffn, fstime=fstime, vtime=vtime )
 
    flon2d, flat2d = np.meshgrid( flon1d, flat1d )
    #sys.exit()

    print( fdat3d.shape, fz3d.shape)
    import matplotlib.pyplot as plt

    #plt.scatter( olon, olat, s=5.0, c=odat ) 
    levels = np.arange( 15, 65, 5 )
#    levels = np.arange( 0, 6500, 500 )
    plt.contourf( flon2d, flat2d, fdat3d[:,:,10], levels=levels )
    plt.colorbar()
    plt.show()

    sys.exit()
    ref4d = read_nc( fn=fn )

    print( ref4d.shape )
    import matplotlib.pyplot as plt
    plt.contourf( ref4d[10,:,:,18], levels=[15, 20, 30 ] )
    plt.show()

INFO = { "TOP": "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test/d4" }

vtime = datetime( 2019, 8, 24, 15, 20, 0 )
fstime = datetime( 2019, 8, 24, 15, 10, 0 )


main( INFO=INFO, vtime=vtime, fstime=fstime )
