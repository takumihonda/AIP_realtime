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
   odir_corrected = os.path.join( path, 'attenuation_corrected_JST', 
                        jstime.strftime('%Y%m%d')  )

   odir_ncorrected = os.path.join( path, 'attenuation_not_corrected_JST', 
                        jstime.strftime('%Y%m%d')  )

   # target data directory
   tdir = os.path.join( path, 'link' )
   os.makedirs( tdir, exist_ok=True )

   jtime_l = []

   h2_list = [ h.name for h in os.scandir( os.path.join( odir_corrected ) ) ]
   for h2 in h2_list:
      m2_list = [ m.name for m in os.scandir( os.path.join( odir_corrected, h2 ) ) ]

      for m2 in m2_list:
         s2_list = [ s.name for s in os.scandir( os.path.join( odir_corrected, h2, m2 ) ) ]
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

      tfile_ac = os.path.join( tdir, 'obs.pawr.ze-ac_{0:}.dat'.format( time.strftime('%Y%m%d-%H%M%S') ) )
      if os.path.exists( tfile_ac ):
         os.remove( tfile_ac )

      if dif_ <= max_dif_sec:
         #print( org_jtime, jtime_ )
         # make a symbolic link

         ofile_ac = os.path.join( odir_corrected, org_jtime.strftime('%H/%M/%S'), 'corrected_zh_polar_00000.bin' )
         print( ofile_ac )
         print( tfile_ac )
         os.symlink( ofile_ac, tfile_ac )


         for typ_, typ2_ in zip( ['ZH', 'VH'], ['ze', 'vr'] ):
             ofile_nac = os.path.join( odir_ncorrected, org_jtime.strftime('%H/%M/%S'), '{0:}.00-00-PPI.RAW-{1:}_MTI.AUTO-02-NORMAL.saitama.dat'.format( org_jtime.strftime('%Y%m%d_%H%M%S'), typ_ ) )
             tfile_nac = os.path.join( tdir, 'obs.pawr.{0:}_{1:}.dat'.format( typ2_, time.strftime('%Y%m%d-%H%M%S') ) )

             if os.path.exists( tfile_nac ):
                os.remove( tfile_nac )
             os.symlink( ofile_nac, tfile_nac )

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
