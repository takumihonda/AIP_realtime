import numpy as np
import sys
import os

from datetime import datetime
from datetime import timedelta

from tools_AIP import read_dat

import matplotlib.pyplot as plt

quick = True
#quick = False

OVERW = True

def read_dat_idxs( INFO, otyp=4002 ):

   tfoot = INFO["time"].strftime('%Y%m%d-%H%M%S')
   fn = os.path.join( INFO["TOP"], INFO["EXP"], INFO["stime0"].strftime('%Y%m%d%H%M%S'),
                      "obsdep",
                      "obsdep_" + tfoot + ".dat")
   
   print( "Get data:", tfoot )
   elms, lons, lats, levs, dats, ombs, qc_l, omas = read_dat( fn )
 
   # radar reflectivity
   x_l = lons[ ( elms == otyp ) & ( qc_l == 0.0 ) ]
   y_l = lats[ ( elms == otyp ) & ( qc_l == 0.0 ) ]
   z_l = levs[ ( elms == otyp ) & ( qc_l == 0.0 ) ]
   ob_l = ombs[ ( elms == otyp ) & ( qc_l == 0.0 ) ]
   dat_l = dats[ ( elms == otyp ) & ( qc_l == 0.0 ) ]
   oa_l = ombs[ ( elms == otyp ) & ( qc_l == 0.0 )]

   azms = np.where( y_l != 0.0, np.arctan2(x_l,y_l)*180.0/np.pi, 0.0 )         # degrees
   rdist2d = np.sqrt( np.square( x_l ) + np.square( y_l ) )
   rdist3d = np.sqrt( np.square( x_l ) + np.square( y_l ) + np.square( z_l ) )
   elvs = np.where( rdist2d > 0.0, np.arctan(z_l/rdist2d)*180.0/np.pi, 0.0 ) # degrees

   azm_i = np.round( azms )
   elv_i = np.round( elvs )
   range_i = np.round( rdist3d / INFO["DR"] )

   return( z_l, ob_l, oa_l, azm_i, elv_i, range_i )

def get_oerr_elv( INFO, otyp=4002, eidx_l=np.arange(0, 60) ):

   z_l, ob_l, oa_l, iazm_l, ielv_l, irange_l = read_dat_idxs( INFO, otyp=otyp )

   oerr_l = np.zeros( eidx_l.shape )
   oerr_l[:] = np.nan

   cnt_l = np.zeros( eidx_l.shape )
   cnt_l[:] = np.nan

   for i, e in enumerate( eidx_l ):
      ob_ = ob_l[ ( ielv_l >= e ) & ( ielv_l < (e + 1) ) ]
      oa_ = oa_l[ ( ielv_l >= e ) & ( ielv_l < (e + 1) ) ]
      oerr_l[i] = np.sqrt( np.mean( ob_ * oa_ ) )
      if len( ob_ ) > 0: 
         cnt_l[i] = len( ob_ )

   return( oerr_l, cnt_l )


def test( INFO, otyp=4002, eidx_l=np.arange(0, 60), oerr_l=np.ones(1), 
          elv=30, r=30, azm=20,
          mode='az' ):

   z_l, ob_l, oa_l, iazm_l, ielv_l, irange_l = read_dat_idxs( INFO, otyp=otyp )

   de = INFO["de"]
   dr = INFO["dr"]
   da = INFO["da"]


   if mode == "az":
      cmax_ = INFO["amax"]
      azm = -999.0
      da = 99999.0

      imax = INFO["rmax"]
      imin = INFO["rmin"]
      di = dr

   elif mode == "el" or mode == "range":
      cmax_ = INFO["emax"]
      imin = 0
      imax = 360

      di = da

      if mode == "el":
         elv = -999.0
         de = 99999.0
      elif mode == "range":
         r = -999.0
         dr = 99999.0

   cor_l = np.zeros( cmax_ + 1 )
   cnt_l = np.zeros( cmax_ + 1 )


   for i in range( imin, imax+di, di ):

      if mode == "az":
         r = i
      elif mode == "el" or mode == "range":
         azm = i

      ob_l_ = ob_l[ ( ielv_l >= elv ) & ( ielv_l < ( elv + de ) ) &
                    ( irange_l >= r ) & ( irange_l < ( r + dr ) ) &
                    ( iazm_l >= azm ) & ( iazm_l < ( azm + da ) ) ]

      oa_l_ = oa_l[ ( ielv_l >= elv ) & ( ielv_l < ( elv + de ) ) &
                    ( irange_l >= r ) & ( irange_l < ( r + dr ) ) &
                    ( iazm_l >= azm ) & ( iazm_l < ( azm + da ) ) ]

      irange_l_ = irange_l[ ( ielv_l >= elv ) & ( ielv_l < ( elv + de ) ) &
                            ( irange_l >= r ) & ( irange_l < ( r + dr ) ) &
                            ( iazm_l >= azm ) & ( iazm_l < ( azm + da ) ) ]

      iazm_l_ = iazm_l[ ( ielv_l >= elv ) & ( ielv_l < ( elv + de ) ) &
                        ( irange_l >= r ) & ( irange_l < ( r + dr ) ) &
                        ( iazm_l >= azm ) & ( iazm_l < ( azm + da ) ) ]
   
      ielv_l_ = ielv_l[ ( ielv_l >= elv ) & ( ielv_l < ( elv + de ) ) &
                        ( irange_l >= r ) & ( irange_l < ( r + dr ) ) &
                        ( iazm_l >= azm ) & ( iazm_l < ( azm + da ) ) ]

      ob_l_ -= np.mean( ob_l_ )
      oa_l_ -= np.mean( oa_l_ )
   
      for i in range( len( ob_l_ ) ):
          if mode == "az":
             didx_l_ =  ( np.abs( iazm_l_ - iazm_l_[i] ) ).astype('int32')
          elif mode == "el":
             didx_l_ =  ( np.abs( ielv_l_ - ielv_l_[i] ) ).astype('int32')
          elif mode == "range":
             didx_l_ =  ( np.abs( irange_l_ - irange_l_[i] ) ).astype('int32')

          didx_l_[ np.abs( didx_l_ ) > cmax_ ] = cmax_
   
          cor_l[didx_l_] += ( ob_l_[i] * oa_l_ + oa_l_[i] * ob_l_ ) * 0.5 
          cnt_l[didx_l_] += 1
   

   return( cor_l, cnt_l )


def take_ave1d( INFO, otyp=4002, amax=10, range_in=10.0, elv_in=5.0, azm_in=10.0,
                dr=10.0, de=5.0, mode='az' ):

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

   azms = np.where( ys != 0.0, np.arctan2(xs,ys)*180.0/np.pi, 0.0 )         # degrees
   rdist2d = np.sqrt( np.square( xs ) + np.square( ys ) )
   rdist3d = np.sqrt( np.square( xs ) + np.square( ys ) + np.square( zs ) )
   elvs = np.where( rdist2d > 0.0, np.arctan(zs/rdist2d)*180.0/np.pi, 0.0 ) # degrees

   azm_i = np.round( azms )
   elv_i = np.round( elvs )
   range_i = np.round( rdist3d / INFO["DR"] )

#   print( "debug" )
#   print( np.max( elv_i), np.min(elv_i) )
#   emax = int( np.max( elv_i ) ) 
#   emin = int( np.min( elv_i ) ) 
#   for e in range( emin, emax ):
#      ob_ = obs[ ( elv_i >= e ) & ( elv_i < (e + 1) ) ]
#      oa_ = oas[ ( elv_i >= e ) & ( elv_i < (e + 1) ) ]
#      oestd_ = np.sqrt( np.mean( ob_ * oa_ ) )
#      print( e, len( ob_ ), oestd_  )
   sys.exit()
  


   range_il = [ range_in ]
   elv_il = [ elv_in ]
   azm_il = [ azm_in ]

   ba1d = np.zeros( amax*2+1)
   ab1d = np.zeros( amax*2+1)

   cnt1d = np.zeros( amax*2+1)
   cor1d = np.zeros( amax*2+1)

   print( "mode", mode )

   if mode == "az":
      il1 = range_il
      il2 = elv_il

      idx1 = range_i
      idx2 = elv_i
      idx3 = azm_i

   elif mode == "el":
      il1 = range_il
      il2 = azm_il

      idx1 = range_i
      idx2 = azm_i
      idx3 = elv_i

   for ir, ie in zip( il1, il2 ):
#       print( ir, ie )
       obs_ = obs[ ( idx1 >= ir ) & 
                   ( idx1 < ( ir + dr ) ) &
                   ( idx2 >= ie ) &
                   ( idx2 < ( ie + de ) ) ]

       oas_ = oas[ ( idx1 >= ir ) & 
                   ( idx1 < ( ir + dr ) ) &
                   ( idx2 >= ie ) &
                   ( idx2 < ( ie + de ) ) ]

       idx3_ = idx3 [ ( idx1 >= ir ) & 
                       ( idx1 < ( ir + dr ) ) &
                       ( idx2 >= ie ) &
                       ( idx2 < ( ie + de ) ) ]
#       print( obs_.shape )

       obs_ -= np.mean( obs_ )
       oas_ -= np.mean( oas_ )
       for i in range( len( obs_ ) ):
           # shift such that it is consistent with az1d
           da_ = ( idx3_ - idx3_[i] ).astype('int32') + amax 
           da_[ da_ > 2*amax ] = 2*amax
           da_[ da_ < 0] = 0
#           print( da_ )
           
           ba1d[da_] +=  obs_[i]*oas_ 
           ab1d[da_] +=  oas_[i]*obs_ 
           cnt1d[da_] += 1
   
   cor1d = np.where( cnt1d > 0.0, 0.5 * ( ba1d + ab1d ),
                    0.0 ) 

   return( cor1d, cnt1d )


 




def plot_oerr( INFO, cor1d, cnt1d, az1d, oerr=1.0, range_in=5.0, dr=1.0, 
               elv_in=1.0, de=1.0 ):

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



###################

def main( stime=datetime(2019,9,10,15,0), etime=datetime(2019,9,10,15,0), otyp=4002,
          stime0=datetime(2019,9,10,15,0), 
          exp="D4_500m_20190910_THIN", range_in=1, elv_in=5, dr=1, de=1, mode="az",
          azm_in=10.0, da=1 ):

   DR = 500.0
   DZ = 500.0
   DA = 2.0 # degres


   nr = 10
   r1d = np.arange( 0, nr*DR, DR)
   r1dh = np.arange( DR*0.5, (nr+0.5)*DR, DR)

   emax = 60
   amax = 10

   rmin = 20
   rmax = 30

   INFO = { "NR":nr, "DR":DR, "DZ":DZ, "DA": DA, "dr":dr, "de":de, "da":da, 
            "emax": emax, "amax":amax,
            "rmax": rmax, "rmin": rmin,
            "R1D":r1d, "R1DH":r1dh, "EXP":exp, "stime":stime, "stime0":stime0,
            "TOP":"/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"  }


   cnt1d = np.zeros( amax*2+1)
   cor1d = np.zeros( amax*2+1)


   eidx_l = np.arange( 0, emax+de, de )
   oerr_l = np.zeros( eidx_l.shape )
   oerr_l[:] = np.nan

   dir_in = "dat/{0:}/{1:05}".format( exp, otyp )
   ofile = "oerr_{0:}_e{1:03}_{2:03}_r{3:03}_{4:03}_a{5:03}_{6:03}_{7:03}.npz".format( mode, elv_in, de, range_in, dr, azm_in, da, int( DR ) )

   os.makedirs( dir_in, exist_ok=True)
   print( dir_in, ofile )

#   file_ext = False
#   try:
#      print( "Found npz file ")
#      data = np.load( os.path.join( dir_in, ofile ) )
#      file_ext = True
#   except:
#      print( "No npz file ")
# 
#   if file_ext and not OVERW:
#      cor1d = data["cor1d"]
#      cnt1d = data["cnt1d"]
#      oerr = data["oerr"]
#   else:
#
#      tcnt = 0
#      time = stime
#      while( time <= etime):
#          INFO["time"] = time
#          try:
#             oerr_l_, cnt_l_ = get_oerr_elv( INFO, otyp=otyp, eidx_l=eidx_l )
#          except:
#             time += timedelta(seconds=30)
#             continue
#
#          if time == stime:
#             oerr_l = np.copy( oerr_l_*cnt_l_ )
#             cnt_l = np.copy( cnt_l_ )
#          else:
#
#             oerr_l = np.where( np.isnan( oerr_l_ ), oerr_l, oerr_l+oerr_l_*cnt_l_ )
#             cnt_l = np.where( np.isnan( cnt_l_ ), cnt_l, cnt_l+cnt_l_ )
#
#          time += timedelta(seconds=30)
#          tcnt += 1
#   
#      print( "debug2" )
#      oerr_l =  np.sqrt( oerr_l/cnt_l ) 

   time = stime
   while( time <= etime):
       INFO["time"] = time
       try:
          cor_l_, cnt_l_ = test( INFO, otyp=otyp, eidx_l=eidx_l, #oerr_l=oerr_l, 
                               mode=mode )
       except:
          time += timedelta(seconds=30)
          continue

       print( cnt_l_ )

       if time == stime:
          cor_l = np.copy( cor_l_ )
          cnt_l = np.copy( cnt_l_ )
       else:
          cor_l += cor_l_
          cnt_l += cnt_l_
       time += timedelta(seconds=30)

   cor_l = np.where( cnt_l > 0, cor_l / cnt_l, np.nan )
   print( cnt_l )
   fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
   ax1.set_ylim( -0.5, 1.0 )
   ax1.plot( cor_l/cor_l[0] )
   ax1.set_xlim( 0, 9 )
   plt.show()

   sys.exit()  


#      print( cor1d )
#      print( cnt1d )
#   
#      cor1d = cor1d[:] / cnt1d[:]
#   
#      oerr = np.sqrt( cor1d[amax] )
#      cor1d = cor1d / cor1d[amax]
#      print( cor1d, cor1d[amax] )
#
# 
#   # 
#   amin = -amax
#   az1d = np.zeros(  amax*2+1 )
#   az1d[0:amax+1] = np.arange( 0, amax+1 )
#   az1d[amax+1:amax*2+1] = np.arange( amin, 0 )
#
#   az1d = np.arange( amin, amax+1 )
#
##   if np.sum(cnt1d) > 00:
##      plot_oerr( INFO, cor1d, cnt1d, az1d, oerr=oerr, range_in=range_in, dr=dr, 
##                 elv_in=elv_in, de=de )
#
#   if not file_ext:
#      np.savez( os.path.join( dir_in, ofile ), cor1d=cor1d, cnt1d=cnt1d, oerr=oerr )

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
#otyp = 4001 # Z


exp = "D4_500m_TEST_DEFAULT_0515_MEAN"

exp = "D4_500m_TEST_DEFAULT_0604_MEAN"
exp = "20201117/D4_500m_H1V1_VR"
exp = "20201117/D4_500m_H1V1"
#exp = "20201117/D4_500m_H1V1_Z"

stime0 = datetime( 2019, 8, 24, 15, 0, 0)


stime = datetime( 2019, 8, 24, 15, 30, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

stime = datetime( 2019, 8, 24, 16, 0, 0 )
etime = stime # DEBUG

stime = datetime( 2019, 8, 24, 15, 40, 30 )
etime = datetime( 2019, 8, 24, 15, 59, 0 )


dr = 1 # [xDR m]
de = 1 # [ degrees ]
da = 1 # [ degrees ]


mode = "az"
mode = "range"
mode = "el"


main( stime=stime, etime=etime, stime0=stime0, otyp=otyp, exp=exp, 
      dr=dr, de=de, da=da, mode=mode )  


