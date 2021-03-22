import os
import sys
from datetime import datetime, timedelta

import numpy as np

def d3_computation_time( top='', dtime_max=2000.0 ):
    dirs = [ f.name for f in os.scandir( top ) ] #if f.is_file() ]

    times = []
    path_l = []

    # Prepare file path list
    for dir_ in dirs:
        fname = [ f.name for f in os.scandir( os.path.join( top, dir_ ) ) ] #if f.is_file() ]
        path_l.append( os.path.join( top, dir_, fname[0] ) )
 
    # Get computation time
    for path in path_l:
          
        if not os.path.isfile( path ):
           break

        with open( path ) as f:
             lines = f.readlines()

             for l in lines:
                 if 'Finish' in l:
                     ftime = datetime( year=int( l[1:5] ), month=int( l[6:8] ), 
                                       day=int( l[9:12] ), hour=int( l[12:14] ), 
                                       minute=int( l[15:17] ), second=int( l[18:20] ) )
                 elif 'Start' in l:
                     stime = datetime( year=int( l[1:5] ), month=int( l[6:8] ), 
                                       day=int( l[9:12] ), hour=int( l[12:14] ), 
                                       minute=int( l[15:17] ), second=int( l[18:20] ) )
                 elif 'Error' in l:
                      break

        dtime =  ( ftime - stime ).total_seconds() 
        if np.abs( dtime ) > dtime_max:
           print( "    Exception ", dtime, path )
        else:
           times.append( dtime )
 
    return( times )


####
top = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/log_from_amemiya/d3/exp'

dtime_max = 10000

typ = 'fcst'
ctimes = d3_computation_time( top=top, dtime_max=dtime_max )
print( '{0:} average: {1:}'.format( typ, np.mean( ctimes ) ) )
print( ctimes )


