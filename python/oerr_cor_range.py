import numpy as np
import sys
import os

from datetime import datetime
from datetime import timedelta

from tools_AIP import read_dat

quick = True
#quick = False

def take_ave1d_range( INFO, otyp=4002, amax=10, range_in=10.0, elv_in=5.0,
                dr=10.0, de=5.0, ):

   tfoot = INFO["time"].strftime('%Y%m%d-%H%M%S')
   fn = os.path.join( INFO["TOP"], INFO["EXP"], INFO["stime0"].strftime('%Y%m%d%H%M%S'),
                      "obsdep",
                      "obsdep_" + tfoot + ".dat")
   
   print( "Get data:", tfoot )
   elms, lons, lats, levs, dats, ombs, qcs, omas = read_dat( fn )
 
   # radar reflectivity
   xs = lons[ elms == otyp ]
   ys = lats[ elms == otyp ]
   zs = levs[ elms == otyp ]
   obs = ombs[ elms == otyp ]
   dats = dats[ elms == otyp ]
   qcs = qcs[ elms == otyp ]
   oas = ombs[ elms == otyp ]

   azms = np.where( ys != 0.0, np.arctan2(xs,ys)*180.0/np.pi, 0.0 )
   rdist2d = np.sqrt( np.square( xs ) + np.square( ys ) )
   rdist3d = np.sqrt( np.square( xs ) + np.square( ys ) + np.square( zs ) )
   elvs = np.where( rdist2d > 0.0, np.arctan(zs/rdist2d)*180.0/np.pi, 0.0 )


   azm_i = np.round( azms )
   elv_i = np.round( elvs )
   range_i = np.round( rdist3d / INFO["DR"] )


   range_il = [ range_in ]
   elv_il = [ elv_in ]

   ba1d = np.zeros( amax*2+1)
   ab1d = np.zeros( amax*2+1)

   cnt1d = np.zeros( amax*2+1)
   cor1d = np.zeros( amax*2+1)

   azm_max = 1

   for ie in elv_il:
#       print( ir, ie )
       obs_ = obs[ ( elv_i >= ie ) &
                   ( elv_i < ( ie + de ) ) ]
       oas_ = oas[ ( elv_i >= ie ) &
                   ( elv_i < ( ie + de ) ) ]
       azm_i_ = azm_i[ ( elv_i >= ie ) &
                       ( elv_i < ( ie + de ) ) ]
       range_i_ = range_i[ ( elv_i >= ie ) &
                       ( elv_i < ( ie + de ) ) ]

#       print( obs_.shape )

       obs_ -= np.mean( obs_ )
       oas_ -= np.mean( oas_ )
       for i in range( len( obs_ ) ):
           # shift such that it is consistent with az1d
           da_ = ( range_i_ - range_i_[i] ).astype('int32') + amax 
           da2_ = np.abs( azm_i_ - azm_i_[i] ).astype('int32') 
           da_[ da_ > 2*amax ] = 2*amax
           da_[ da_ < 0] = 0
 
           da_[ da2_ > azm_max ] = 0

           ba1d[da_] +=  obs_[i]*oas_ 
           ab1d[da_] +=  oas_[i]*obs_ 
           cnt1d[da_] += 1

   cor1d = np.where( cnt1d > 0.0, 0.5 * ( ba1d + ab1d ),
                    0.0 ) 

   return( cor1d, cnt1d )








def main_range( stime=datetime(2019,9,10,15,0), etime=datetime(2019,9,10,15,0), otyp=4002,
          stime0=datetime(2019,9,10,15,0), 
          exp="D4_500m_20190910_THIN", range_in=1.0, elv_in=5.0, dr=1.0, de=1.0 ):

   DR = 500.0
   DZ = 500.0


   nr = 10
   r1d = np.arange( 0, nr*DR, DR)
   r1dh = np.arange( DR*0.5, (nr+0.5)*DR, DR)


   INFO = { "NR":nr, "DR":DR, "DZ":DZ,
            "R1D":r1d, "R1DH":r1dh, "EXP":exp, "stime":stime, "stime0":stime0,
            "TOP":"/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"  }


   amax = 10
   cnt1d = np.zeros( amax*2+1)
   cor1d = np.zeros( amax*2+1)


   ofile = "dat/oerr_range_e{0:}_{1:}_r{2:}_{3:}.npz".format( elv_in, de, range_in, dr)
   file_ext = False



   try:
      print( "Found npz file ")
      data = np.load( ofile )
      file_ext = True
      cor1d = data["cor1d"]
      cnt1d = data["cnt1d"]
      oerr = data["oerr"]
   except:
      print( "No npz file ")
 
      tcnt = 0
      time = stime
      while( time <= etime):
          INFO["time"] = time
          try:
             cor1d_, cnt1d_ = take_ave1d_range( INFO, otyp=otyp, amax=amax, 
                                  range_in=range_in, elv_in=elv_in, 
                                  dr=dr, de=de,  )
             print( tcnt, np.sum(cnt1d_), "\n" )
      
             cor1d += cor1d_
             cnt1d += cnt1d_
          except:
             print(" Failed to opne ")
   
          time += timedelta(seconds=30)
          tcnt += 1
   
   
      print( cor1d )
      print( cnt1d )
   
      cor1d = cor1d[:] / cnt1d[:]
   
      oerr = np.sqrt( cor1d[amax] )
      cor1d = cor1d / cor1d[amax]
      print( cor1d, cor1d[amax] )
 
   # 
   amin = -amax
   az1d = np.zeros(  amax*2+1 )
   az1d[0:amax+1] = np.arange( 0, amax+1 )
   az1d[amax+1:amax*2+1] = np.arange( amin, 0 )

   az1d = np.arange( amin, amax+1 )

   if np.sum(cnt1d) > 100:
      plot_oerr( INFO, cor1d, cnt1d, az1d, oerr=oerr, range_in=range_in, dr=dr, 
                 elv_in=elv_in, de=de )

   if not file_ext:
      np.savez( ofile, cor1d=cor1d, cnt1d=cnt1d, oerr=oerr )

def plot_oerr( INFO, cor1d, cnt1d, az1d, oerr=1.0, range_in=5.0, dr=1.0, 
               elv_in=1.0, de=1.0 ):

   import matplotlib.pyplot as plt
   fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
   fig.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.8, )

   ardist = ( range_in + dr*0.5 ) * INFO["DR"]
   print( ardist )
   dist1d = az1d * INFO["DR"] * 0.001 # km

   #ax1.plot( az1d, cor1d, )
   #ax1.set_xlim( -7, 7  )
   #xlab1 = 'Azimuth angle ($^\circ$)'
   ax1.plot( dist1d, cor1d, lw=3.0, color='k' )
   xmin_ = -3.0
   xmax_ = 3.0
   xlab1 = "Separation distance in the range direction (km)"

   ax1.set_xlim( xmin_, xmax_  )


#   ax1.set_xticks( np.arange( -10, 11, 1 ) )

   ylab1 = "Correlations"

   ymin_ = 0.0
   ymax_ = 1.0
   ax1.set_ylim( ymin_, ymax_ )
   ax1.set_xlabel( xlab1, fontsize=12 )
   ax1.set_ylabel( ylab1, fontsize=12 )

   ax1.vlines( x=np.arange( xmin_, xmax_+0.5, 0.5), ymin=ymin_, ymax=ymax_,
               ls='dashed', lw=0.5, color='k')

   ax1.hlines( y=[ 0.2 ], xmin=xmin_, xmax=xmax_,
               ls='dashed', lw=0.5, color='k')

   bbox = {'facecolor':'w', 'alpha':0.95, 'pad':2,
           'edgecolor':'w' }

   note = 'Elevation: {0:}-{1:}$^\circ$'.format ( elv_in, elv_in+de )  #\n\
#           Range: {4:}-{5:}km\n\
#           $\sigma_o$={6:.1f}dBZ\n\
#           N={7:.0f}'.format(
#            elv_in, elv_in+de, 
#            np.sin( np.deg2rad(elv_in+de) )*dist_max*0.001,
#            dist_min*0.001, dist_max*0.001,
#            oerr,
#            cnt1d[0] )

   ax1.text( 0.99, 0.99, note,
             fontsize=10, transform=ax1.transAxes,
             horizontalalignment='right',
             verticalalignment='top',
             bbox=bbox )
   
   if otyp == 4001:
      nvar = "PAWR X-band Radar Reflectivity\n"
   tit = nvar + "Estimated Observation Error Correlations"
   ax1.text( 0.5, 1.13, tit,
             fontsize=14, transform=ax1.transAxes,
             horizontalalignment='center',
             verticalalignment='bottom' )


   plt.show()
   sys.exit()





   print( np.sin( np.deg2rad(elv_in) )*dist_min )
   print( np.sin( np.deg2rad(elv_in+de) )*dist_max )


   ofig = "oerr_az_e{0:}_{1:}_r{2:}_{3:}.png".format( elv_in, de, range_in, dr)

   print( ofig )
   if quick:
      plt.show()
   else:
      odir = "png/oerr"
      os.makedirs( odir, exist_ok=True)
      plt.savefig( os.path.join(odir, ofig), 
                   bbox_inches="tight", pad_inches = 0.1)
      plt.clf()
      plt.close('all')




##################

otyp = 4002 # vr
#otyp = 4004 # zero Z
otyp = 4001 # Z


exp = "D4_500m_TEST_DEFAULT_0515_MEAN"

exp = "D4_500m_TEST_DEFAULT_0604_MEAN"

stime0 = datetime( 2019, 8, 24, 15, 0, 0)


stime = datetime( 2019, 8, 24, 15, 30, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

stime = datetime( 2019, 8, 24, 16, 0, 0 )
etime = stime # DEBUG



dr = 5.0
de = 1.0

range_inl = np.arange( 10, 100+dr, dr )
elv_inl = np.arange( 1, 20+de, de )

for elv_in in elv_inl:
   for range_in in range_inl:

      main_range( stime=stime, etime=etime, stime0=stime0, otyp=otyp, exp=exp, 
               range_in=range_in, elv_in=elv_in, dr=dr, de=de )  

