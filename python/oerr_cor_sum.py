import numpy as np
import sys
import os

from datetime import datetime
from datetime import timedelta

from tools_AIP import read_dat

quick = True
quick = False


def main( stime=datetime(2019,9,10,15,0), etime=datetime(2019,9,10,15,0), otyp=4002,
          stime0=datetime(2019,9,10,15,0), 
          exp="D4_500m_20190910_THIN", ):

   DR = 500.0
   DZ = 500.0

   dr = 5.0
   de = 1.0
   
   range_inl = np.arange( 10, 100+dr, dr )
   elv_inl = np.arange( 1, 20+de, de )

   nr = len( range_inl )
   ne = len( elv_inl )
   print( nr )



   INFO = { "DR":DR, "DZ":DZ,
            "EXP":exp, "stime":stime, "stime0":stime0,
            "TOP":"/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"  }


   amax = 10
   cnt1d = np.zeros( (ne, nr, amax*2+1) )
   cor1d = np.zeros( (ne, nr, amax*2+1) )

   cnt1d[:] = np.nan
   cor1d[:] = np.nan

   for e, elv_in in enumerate( elv_inl ):
      for r, range_in in enumerate( range_inl ):

         ofile = "dat/oerr_az_e{0:}_{1:}_r{2:}_{3:}.npz".format( elv_in, de, range_in, dr)
         file_ext = False
         try:
#            print( "Found npz file ")
            data = np.load( ofile )
            file_ext = True
            cor1d[e,r,:] = data["cor1d"]
            cnt1d[e,r,:] = data["cnt1d"]
            oerr = data["oerr"]
         except:
            print( "No npz file ")
            print( ofile )
            sys.exit()
 
   # 
   amin = -amax
   az1d = np.zeros(  amax*2+1 )
   az1d[0:amax+1] = np.arange( 0, amax+1 )
   az1d[amax+1:amax*2+1] = np.arange( amin, 0 )

   az1d = np.arange( amin, amax+1 )

   econst = 10
   rconst = None

   econst = None
   rconst = 10

   plot_oerr( INFO, cor1d, cnt1d, az1d, oerr=oerr, range_in=range_in, dr=dr, 
                 elv_in=elv_in, de=de, range_inl=range_inl, elv_inl=elv_inl,
                 econst=econst, rconst=rconst )


def plot_oerr( INFO, cor1d, cnt1d, az1d, oerr=1.0, range_in=5.0, dr=1.0, 
               elv_in=1.0, de=1.0, range_inl=[], elv_inl=[], econst=10, rconst=10 ):

   import matplotlib.pyplot as plt
   import matplotlib.cm as cm

   fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
   fig.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.8, )

   ardist = ( range_in + dr*0.5 ) * INFO["DR"]
   print( ardist )
   dist1d = ardist * np.sin( np.deg2rad(az1d) ) * 0.001 # km
   dist_min = range_in*INFO["DR"]
   dist_max = ( range_in + dr )*INFO["DR"]

   #ax1.plot( az1d, cor1d, )
   #ax1.set_xlim( -7, 7  )
   #xlab1 = 'Azimuth angle ($^\circ$)'

   itmax = len( range_inl ) 

#   range_inl = [ 20.0 ]
#   itmax = len( elv_inl )

   if econst is not None:
      note = 'Elevation: {0:}-{1:}$^\circ$'.format( elv_inl[e], elv_inl[e] + de )
      e = econst
      for r in range( len( range_inl) ):
         ref = r 
         dist1_ = ( range_inl[r]  ) * INFO["DR"] * 0.001
         dist2_ = ( range_inl[r] + dr ) * INFO["DR"] * 0.001
   
         if np.any( np.isnan(cor1d[e,r,:] ) ):
            label = None
         else:
            label = "Range: {0:}-{1:}km".format( dist1_, dist2_ ) 
   
         ardist = ( range_inl[r] + dr*0.5 ) * INFO["DR"]
         dist1d = ardist * np.sin( np.deg2rad(az1d) ) * 0.001 # km
         ax1.plot( dist1d, cor1d[e,r,:], lw=3.0, color=cm.jet(ref/itmax),
                   label=label )

   elif rconst is not None:
      r = rconst
      dist1_ = ( range_inl[r]  ) * INFO["DR"] * 0.001
      dist2_ = ( range_inl[r] + dr ) * INFO["DR"] * 0.001
      note = 'Range: {0:}-{1:}km'.format( dist1_, dist2_ )
      for e in range( len( elv_inl) ):
         ref = e 

         if np.any( np.isnan(cor1d[e,r,:] ) ):
            label = None
         else:
            label = "Elevation: {0:}-{1:}".format( elv_inl[e], elv_inl[e]+de ) 

         ardist = ( range_inl[r] + dr*0.5 ) * INFO["DR"]
         dist1d = ardist * np.sin( np.deg2rad(az1d) ) * 0.001 # km
         ax1.plot( dist1d, cor1d[e,r,:], lw=3.0, color=cm.jet(ref/itmax),
                   label=label )


   ax1.legend( loc='upper left', fontsize=8 )

   ylab1 = "Correlations"

   ymin_ = 0.0
   ymax_ = 1.0
   ax1.set_ylim( ymin_, ymax_ )
   ax1.set_ylabel( ylab1, fontsize=12 )

   xlab1 = 'Azimuth angle ($^\circ$)'
   xlab1 = 'Separation distance in azimuth (km)'

   xmin_ = -4.0
   xmax_ = 4.0
   ax1.set_xlabel( xlab1, fontsize=12 )
   ax1.set_xlim( xmin_, xmax_  )

   ax1.hlines( y=[ 0.2 ], xmin=xmin_, xmax=xmax_,
               ls='dashed', lw=0.5, color='k')

   bbox = {'facecolor':'w', 'alpha':0.95, 'pad':2,
           'edgecolor':'w' }


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





   xlab1 = "Horizontal distance (km)"
#   ax1.set_xticks( np.arange( -10, 11, 1 ) )


   ax1.vlines( x=np.arange( xmin_, xmax_+0.5, 0.5), ymin=ymin_, ymax=ymax_,
               ls='dashed', lw=0.5, color='k')


   ax2 = ax1.twiny()
   c = 'b'
   xlab2 = 'Azimuth angle ($^\circ$)'
   ax2.plot( az1d , cor1d, color=c, lw=0.0 )
   ax2.set_xticks( az1d ) #np.arange( -4, 4.5, 0.5 ) )
   ax2.set_xlabel( xlab2 )

   xmin2_ = np.rad2deg( np.arcsin( xmin_ / ardist * 1000.0 ) )
   xmax2_ = np.rad2deg( np.arcsin( xmax_ / ardist * 1000.0 ) )
   print( xmin2_, xmax2_ )
   ax2.set_xlim( xmin2_, xmax2_ )

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
exp = "20201117/D4_500m_H1V1"

stime0 = datetime( 2019, 8, 24, 15, 0, 0)


stime = datetime( 2019, 8, 24, 15, 30, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

#stime = datetime( 2019, 8, 24, 16, 0, 0 )
#etime = stime # DEBUG




main( stime=stime, etime=etime, stime0=stime0, otyp=otyp, exp=exp, )
