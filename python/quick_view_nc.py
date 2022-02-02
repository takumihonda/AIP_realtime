from netCDF4 import Dataset
import numpy as np
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

    return( dat, lon, lat )

def read_nc( fn='' ):

    nc = Dataset( fn, "r", format="NETCDF4" )
    print( nc )
    var4d = nc.variables['Reflectivity'][:]
    lon1d = nc.variables['Longitude'][:]
    z1d = nc.variables['Height'][:]
    nc.close()

    print( lon1d )
    print( z1d )

    return( var4d )

def main():
    #fn = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test/d4/dafcst_nc/20190825-001000.nc'

#    fn = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test/d4/dafcst_nc/20190824-150700.nc'

    odz = 400

    fn = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/tmp/pawr_obs_20190824-151200.nc'
    odat, olon, olat = read_obs_nc( fn=fn, dz=odz )
    import matplotlib.pyplot as plt

    plt.scatter( olon, olat, s=1.0, c=odat )
    plt.show()

    sys.exit()
    ref4d = read_nc( fn=fn )

    print( ref4d.shape )
    import matplotlib.pyplot as plt
    plt.contourf( ref4d[10,:,:,18], levels=[15, 20, 30 ] )
    plt.show()

main()
