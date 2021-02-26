import numpy as np
import sys
import os

from datetime import datetime
from datetime import timedelta

from tools_AIP import read_dat, read_oerr_npz

quick = True
quick = False

AVE = False



###################

def main( otyp=4002,
          exp="D4_500m_20190910_THIN", range_inl=[1.0], elv_inl=[5.0], dr=1.0, de=1.0, mode="az",
          azm_inl=[10.0], da=1.0, 
          rskip=5 ):


   DR = 500.0

   import matplotlib.pyplot as plt
   import matplotlib.cm as cm
   fig, (ax1) = plt.subplots(1, 1, figsize=(7,4))
   fig.subplots_adjust(left=0.1, bottom=0.12, right=0.98, top=0.9, )

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

      xmin = -4.0
      xmax = -xmin

      # plot distance array
      dx = 500.0 # m # >= DR
      np_ = 40
      min_ = -dx * np_
      max_ = dx * np_ + dx
      pdist1d = np.arange( min_, max_, dx )

      cor1d = np.zeros( pdist1d.shape )
      cnt1d = np.zeros( pdist1d.shape )
      oerr = 0.0

      range0 = range_inl[0]*DR/1000.0
      #xlab = r'Azimuth angle separation ($^\circ$)'

      if mode == "az":
         xlab = r'Azimuth separation distance (km)'
      elif mode == "el":
         xlab = r'Elevation separation distance (km)'
 
      tit = "Estimated observation-error correlations for\n"
      if otyp == 4001:   # Z
         tit += "MP-PAWR radar reflectivity"
         otyp_ = "Z"
         unit_ = "dBZ"
      elif otyp == 4001: # vr
         tit += "MP-PAWR Doppler velocity"
         otyp_ = "VR"
         unit_ = "m/s"

     
      ylab = 'Correlation'
      lw = 4.0
      for i, range_in_ in enumerate( range_inl ):
          ii = i + 1
          range_inl_ = [ range_in_ ]

          cor1d_, cnt1d_, oerr_, x1d_ = read_oerr_npz( range_inl=range_inl_, elv_inl=elv_inl, dr=dr, de=de, mode=mode,
                                     azm_inl=azm_inl, da=da, DR=DR, 
                                     exp=exp, otyp=otyp )

          x1d_ = np.sin( np.deg2rad( x1d_ ) ) * range_in_ * DR # [m]
          idx1d_ = ( np.round( x1d_ / dx ) + np_ ).astype( int ) # index for plot distance array


          len2 = len( cor1d_ ) // 2 
          cor1d[ idx1d_] += cor1d_ * cnt1d_
          cnt1d[ idx1d_ ] += cnt1d_
          oerr += oerr_ * cnt1d_[ len2 ]
          range1 = range_in_ * DR * 0.001 # km

          if ii % rskip == 0 or ii == len( range_inl ):
             print( "plot", i )
             range1 = range_in_ * DR / 1000.0

             cor1d = cor1d / cnt1d

#             # dgree
#             x1d = x1d / cnt
# 
#             # distance
#             x1d = np.sin( np.deg2rad( x1d ) ) * range1 

             # plot
             #if not np.isnan( cor1d ).any(): # and cnt1d[len2] > 20000: 
             pcnt += 1 

             label = r'Range: {0:.1f}-{1:.1f} km, $\sigma_o$:{2:.1f} {3:}'.format( range0, 
                                                        range1, 
                                                        oerr/cnt1d[np_],
                                                        unit_ )
             #print( "cor" , cor1d )
             #print( "cnt", cnt1d )

             # Remove NaN
             pcor1d = cor1d[ ~np.isnan(cor1d) ]
             pdist1d_ = pdist1d[ ~np.isnan(cor1d) ] * 0.001

             if pmax > 0:
               lc = cm.jet( pcnt/pmax )
             else:
               lc = 'k'

             ax1.plot( pdist1d_, pcor1d, label=label, 
                       color=lc, lw=lw )

             # store distance for the next plot
             range0 = range1 + DR * 0.001

             # reset count
             cor1d[:] = 0.0
             cnt1d[:] = 0
             oerr = 0.0
             cnt = 0


      ax1.legend()
      ax1.set_xlabel( xlab, fontsize=11 )
      ax1.set_ylabel( ylab, fontsize=11 )
      ax1.set_xlim( xmin, xmax )
      xlabs = np.arange( xmin, xmax+0.5, 0.5 )
      ax1.set_xticks( xlabs )
      ax1.grid( lw=0.5, ls='dashed' )


      ax1.text( 0.5, 1.01, tit,
                fontsize=12, transform=ax1.transAxes,
                ha='center',
                va='bottom',
              )

   
      odir = "png/oerr/{0:}_{1:}".format( exp, otyp_ )
      ofig = "oerr_{0:}_e{1:03}_r{2:03}_a{3:03}_dr{4:03}_skip{5:03}.png".format( mode, de, dr, da, int( DR ), rskip ) 
   
      print( odir, ofig )
      if quick:
         plt.show()
      else:
         os.makedirs( odir, exist_ok=True)
         plt.savefig( os.path.join(odir, ofig), 
                      bbox_inches="tight", pad_inches = 0.1)
         plt.clf()
         plt.close('all')
   



##################

otyp = 4002 # vr
#otyp = 4004 # zero Z
otyp = 4001 # Z


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
rskip = 32
rskip = 64
rskip = 128


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
   

