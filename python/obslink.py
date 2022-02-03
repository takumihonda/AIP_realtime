#import numpy as np
from datetime import datetime, timedelta
import os
import sys

def main( stime=datetime( 2019, 8, 24, 15, 0, 0 ), 
          etime=datetime( 2019, 8, 24, 15, 0, 0 ), 
          dt=timedelta(seconds=30), 
          max_dif_sec=29 ):

   #path = os.getcwd()
   path = "/vol0003/hp150019/u01161/SCALE-LETKF/scale-5.4.3/SCALE-LETKF-rt/obs_attenuation_corrected"

   jstime = stime + timedelta( hours=9 ) 
   # original data directory
   odir = os.path.join( path, 'org_JST', 
                        jstime.strftime('%Y%m%d') )

   # target data directory
   tdir = os.path.join( path, '../obs' )

   jtime_l = []

   h2_list = [ h.name for h in os.scandir( os.path.join( odir ) ) ]
   for h2 in h2_list:
      m2_list = [ m.name for m in os.scandir( os.path.join( odir, h2 ) ) ]

      for m2 in m2_list:
         s2_list = [ s.name for s in os.scandir( os.path.join( odir, h2, m2 ) ) ]
         for s2 in s2_list:
            jtime_ = datetime( year=jstime.year, month=jstime.month, day=jstime.day,
                               hour=int( h2 ), minute=int( m2 ), second=int( s2 ) )
            jtime_l.append( jtime_ )

   time = stime
   while time <= etime:
      #print( len( jtime_l ) )

      jtime = time + timedelta( hours=9 )

      dif_l = []
      for jtime_ in jtime_l:
         dif_l.append( ( abs( jtime_ - jtime ).seconds ) )
      idx = dif_l.index( min( dif_l ) )

      org_jtime = jtime_l[idx]
      dif_ = dif_l[idx]

      tfile = os.path.join( tdir, 'obs.pawr.ze-ac_{0:}.dat'.format( time.strftime('%Y%m%d-%H%M%S') ) )
      if os.path.exists( tfile ):
         os.remove( tfile )

      if dif_ <= max_dif_sec:
         #print( org_jtime, jtime_ )
         # make a symbolic link

         ofile = os.path.join( odir, org_jtime.strftime('%H/%M/%S'), 'corrected_zh_polar_00000.bin' )
         #print( ofile )
         #print( tfile )
         os.symlink( ofile, tfile )

         jtime_l.remove( org_jtime )
      else:
         print( "No link ", "t:", jtime, "o:", org_jtime, dif_ )
         

      time += dt 


##########
stime = datetime( 2019, 8, 24, 15, 0, 0 )
etime = datetime( 2019, 8, 24, 16, 1, 0 )

max_dif_sec = 29
max_dif_sec = 15
main( stime=stime, etime=etime, max_dif_sec=max_dif_sec )
