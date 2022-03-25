
import numpy as np
import sys
from datetime import datetime, timedelta
from tools_AIP import get_grads_JMA_snap, get_JMA_lonlat, read_nc_lonlat, read_CZ

lon2d, lat2d = get_JMA_lonlat()

#mlon1d, mlat1d, _, _, _ = read_nc_lonlat()
#mlon2d, mlat2d = np.meshgrid( mlon1d, mlat1d )

#lons = np.min( mlon1d )
#lone = np.max( mlon1d )

#lats = np.min( mlat1d )
#late = np.max( mlat1d )

late = 36.331504821777344
lats = 35.3876838684082

lone = 140.1912841796875
lons = 139.02670288085938


# UTC
stime = datetime( 2020, 8, 24, 3 )
etime = datetime( 2020, 9,  7, 3 )

stime = datetime( 2020, 7, 31, 0 )
etime = datetime( 2020, 8,  7, 0 )

time_l = []
rmax_l = []
rarea_l = []

rmin = 30.0 # mm/h
rmin = 20.0 # mm/h
rmin = 1.0 # mm/h
#rmin = 5.0 # mm/h

#of = "JMA_rmax_rmin{0:.0f}.npz".format( rmin )    
of = "JMA_rmax_rmin{0:.0f}_{1:}_{2:}.npz".format( rmin, 
                                        stime.strftime( '%Y%m%d' ),    
                                        etime.strftime( '%Y%m%d' ),    
                                      )
print( of )
#sys.exit

time = stime   
while time <= etime:

    data_ = get_grads_JMA_snap( time )
#
#    print( data_.shape )

    data_[ ( lon2d < lons ) | ( lon2d > lone ) | ( lat2d < lats) | ( lat2d > late) ] = np.nan
#    import matplotlib.pyplot as plt
#    plt.contourf( lon2d, lat2d, data_ )
#    plt.show()
    rmax = np.nanmax( data_ )

    data2_ = np.where( data_ >= rmin, 1.0, 0.0 )
    rarea = np.nansum( data2_ )

    time_l.append( time )
    rmax_l.append( rmax )
    rarea_l.append( rarea )
    print( time )
#    sys.exit()

    time += timedelta( minutes=10 )

print( rmax_l )
print( time_l )

np.savez( of, rmax_l=rmax_l, time_l=time_l, rarea_l=rarea_l, rmin=rmin )

