#import numpy as np
from datetime import datetime, timedelta
import os
import sys

def main( stime=datetime( 2019, 8, 24, 15, 0, 0 ), 
          etime=datetime( 2019, 8, 24, 15, 0, 0 ), 
          dt=timedelta(seconds=30), 
          max_dif_sec=29 ):

   path = os.getcwd()

   #path = "/vol0003/hp150019/u01161/SCALE-LETKF/scale-5.4.3/SCALE-LETKF-rt/obs/obs_attenuation_corrected"

   jstime = stime + timedelta( hours=9 ) 
   # original data directory
   odir = os.path.join( path, 'attenuation_not_corrected_JST', 
                        jstime.strftime('%Y%m%d'), 'raw' )

   h2_list = [ h.name for h in os.scandir( os.path.join( odir ) ) ]

   print( odir, h2_list )
   for h2 in h2_list:
       raw_dir = os.path.join( odir, h2 )
       print( raw_dir )
       file_list = [ f.name for f in os.scandir( os.path.join( raw_dir ) ) ]
       print( file_list )

       for file_ in file_list:
           y4_ = file_[0:4] 
           m2_ = file_[4:6] 
           d2_ = file_[6:8] 
           h2_ = file_[9:11] 
           n2_ = file_[11:13] 
           s2_ = file_[13:15] 

           print( y4_, m2_, d2_, h2_, n2_, s2_ )

           dir_ = os.path.join( path, 'attenuation_not_corrected_JST', 
                                jstime.strftime('%Y%m%d'), h2_, n2_, s2_ )
           print( dir_, file_ )
           os.makedirs( dir_, exist_ok=True )

           ofile_ = os.path.join( raw_dir, file_ )
           tfile_ = os.path.join( dir_, file_ )

           print( ofile_ )
           print( tfile_ )
           os.symlink( ofile_, tfile_ )

   sys.exit()

   # target data directory
   tdir = os.path.join( path, 'link' )
   os.makedirs( tdir, exist_ok=True )

#   print( tdir )
#   sys.exit()

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
         print( ofile )
         print( tfile )
#         os.symlink( ofile, tfile )

         jtime_l.remove( org_jtime )
      else:
         print( "No link ", "t:", jtime, "o:", org_jtime, dif_ )
         

      time += dt 


##########
stime = datetime( 2019, 8, 24, 15, 0, 0 )
etime = datetime( 2019, 8, 24, 16, 1, 0 )

#max_dif_sec = 29
max_dif_sec = 15
main( stime=stime, etime=etime, max_dif_sec=max_dif_sec )
