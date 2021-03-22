import os
import sys
from datetime import datetime, timedelta

import numpy as np

def mon2mm( month='Mar' ):
    month_ = datetime.strptime( month, "%b")
    return( int( month_.month ) )

def ncepobs_time( fn='', dtime_max=2000.0 ):

    times = []

    print( fn )

    with open( fn ) as f:
         lines = f.readlines()

         for l in lines:
             data = l.split()
             # get time ( JST )
             gtime = datetime( year=2021, month=mon2mm( data[5] ), 
                        day=int( data[6] ), hour=int( data[7][0:2] ), 
                        minute=int( data[7][3:5] ) )
             print( gtime )
             # data time ( JST )
             dtime = datetime( year=int( data[8][0:4] ), 
                               month=int( data[8][4:6] ),
                               day=int( data[8][6:8] ), 
                               hour=int( data[8][8:10] )  ) + timedelta( hours=9)
             print( dtime )
             lag = gtime - dtime
             times.append( lag.total_seconds() )
             

    return( times )


####
fn = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/log_from_amemiya/ncepobs_gdas_letkf/log'

typ = 'ncepobs'

ctimes = ncepobs_time( fn=fn )
print( '{0:} average: {1:}'.format( typ, np.mean( ctimes ) ) )


