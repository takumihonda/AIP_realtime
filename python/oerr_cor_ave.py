import numpy as np
import sys
import os

from datetime import datetime
from datetime import timedelta

from tools_AIP import read_dat, read_oerr_npz

quick = True
#quick = False

AVE = False



###################

def main( otyp=4002,
          exp="D4_500m_20190910_THIN", range_inl=[1.0], elv_inl=[5.0], dr=1.0, de=1.0, mode="az",
          azm_inl=[10.0], da=1.0, 
          rskip=5 ):


   DR = 500.0

   import matplotlib.pyplot as plt
   import matplotlib.cm as cm
   fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
   fig.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.8, )

   ax1.set_ylim( 0.0, 1.0 )

   if AVE:

      cor1d, cnt1d, oerr, x1d = read_oerr_npz( range_inl=range_inl, elv_inl=elv_inl, dr=dr, de=de, mode=mode,
                                azm_inl=azm_inl, da=da, DR=DR, exp=exp, oytp=otyp )
   
   
      print( cor1d )
      ax1.plot( x1d, cor1d, )
      plt.show()

   else:

      pcnt = 0
      pmax = len( range_inl ) // rskip

      xmin = -5.0
      xmax = -xmin

      range0 = range_inl[0]*DR/1000.0
      #xlab = r'Azimuth angle separation ($^\circ$)'
      xlab = r'Azimuth separation distance (km)'
      ylab = 'Correlation'
      for i, range_in_ in enumerate( range_inl ):
          ii = i + 1
          range_inl_ = [ range_in_ ]

          cor1d_, cnt1d_, oerr_, x1d_ = read_oerr_npz( range_inl=range_inl_, elv_inl=elv_inl, dr=dr, de=de, mode=mode,
                                     azm_inl=azm_inl, da=da, DR=DR, 
                                     exp=exp, otyp=otyp )
          print( cor1d_[10], i )

          if ii == 1:
             cor1d = np.copy( cor1d_ )
             cnt1d = np.copy( cnt1d_ )
             x1d = np.copy( x1d_ )
             oerr = np.copy( oerr_ )
             cnt = 1
             len2 = len( cnt1d ) // 2
          else:
             cnt += 1
             cor1d += cor1d_ 
             cnt1d += cnt1d_ 
             oerr += oerr_ 
             x1d += x1d_ 

             if ii % rskip == 0:
                print( "plot", i )
                range1 = range_in_ * DR / 1000.0

                cor1d = cor1d / cnt
                # dgree
                x1d = x1d / cnt
 
                # distance
                x1d = np.sin( np.deg2rad( x1d ) ) * range1 

                # plot
                if not np.isnan( cor1d ).any(): # and cnt1d[len2] > 20000: 
                   pcnt += 1 

                   label = '{0:.1f}-{1:.1f}km, err:{2:.1f}, cnt:{3:.0f}'.format( range0, 
                                                                range1,
                                                                oerr/cnt, cnt1d[len2] )
                   print( cor1d )
                   print( cnt1d )
                   ax1.plot( x1d, cor1d, label=label, 
                             color=cm.jet( pcnt/pmax ) )


                # store distance for the next plot
                range0 = range1 + DR / 1000.0

                # reset count
                cor1d = 0.0
                cnt1d = 0.0
                oerr = 0.0
                cnt = 0


      ax1.legend()
      ax1.set_xlabel( xlab, fontsize=12 )
      ax1.set_ylabel( ylab, fontsize=12 )
      ax1.set_xlim( xmin, xmax )
      plt.show()

   sys.exit()


def plot_oerr( INFO, cor1d, cnt1d, az1d, oerr=1.0, range_in=5.0, dr=1.0, 
               elv_in=1.0, de=1.0 ):

   import matplotlib.pyplot as plt
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
   ax1.plot( dist1d, cor1d, lw=3.0, color='k' )
   xmin_ = -3.0
   xmax_ = 3.0
   xlab1 = "Horizontal distance (km)"

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

   bbox = {'facecolor':'w', 'alpha':0.95, 'pad':2,
           'edgecolor':'w' }

   print( np.sin( np.deg2rad(elv_in) )*dist_min )
   print( np.sin( np.deg2rad(elv_in+de) )*dist_max )

   note = 'Elevation: {0:}-{1:}$^\circ$\n\
           (Height: {2:.1f}-{3:.1f}km)\n\
           Range: {4:}-{5:}km\n\
           $\sigma_o$={6:.1f}dBZ\n\
           N={7:.0f}'.format(
            elv_in, elv_in+de, 
            np.sin( np.deg2rad(elv_in) )*dist_min*0.001,
            np.sin( np.deg2rad(elv_in+de) )*dist_max*0.001,
            dist_min*0.001, dist_max*0.001,
            oerr,
            cnt1d[0] )

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
exp = "20201117/D4_500m_H1V1_Z"

stime0 = datetime( 2019, 8, 24, 15, 0, 0)


stime = datetime( 2019, 8, 24, 15, 30, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

#stime = datetime( 2019, 8, 24, 16, 0, 0 )
#etime = stime # DEBUG


de = 1
da = 1

dr = 1
rskip = 10

dr = 1
rskip = 8
rskip = 16


range_inl = np.arange( 40, 120+dr, dr )
#range_inl = np.arange( 40, 100+dr, dr )
elv_inl = np.arange( 1, 60+de, de )
azm_inl = np.arange( 0, 360, da )


mode = "az"
#mode = "el"
#mode = "ra"


main( otyp=otyp, exp=exp, 
      range_inl=range_inl, elv_inl=elv_inl, azm_inl=azm_inl, dr=dr, de=de, da=da, mode=mode,
      rskip=rskip,  )  
   

