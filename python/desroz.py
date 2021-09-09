import numpy as np
import sys
import os

from datetime import datetime
from datetime import timedelta

from tools_AIP import read_dat, get_nbin, desroz_diag_R

quick = True
#quick = False

HIST = False
HIST = True

def get_dat( INFO, otyp=4002 ):

   tfoot = INFO["time"].strftime('%Y%m%d-%H%M%S')
   fn = os.path.join( INFO["TOP"], INFO["EXP"], INFO["stime0"].strftime('%Y%m%d%H%M%S'),
                      "obsdep",
                      "obsdep_" + tfoot + ".dat")
   
   print("Get data:",tfoot)
   elms, lons, lats, levs, dats, ombs, qcs, omas = read_dat( fn )
 
   r_dxs = lons[ elms == otyp ]
   r_dys = lats[ elms == otyp ]
   r_dzs = levs[ elms == otyp ]
   r_ombs = ombs[ elms == otyp ]
   r_omas = omas[ elms == otyp ]
   r_dats = dats[ elms == otyp ]
   r_qcs = qcs[ elms == otyp ]


   return( r_dxs.tolist(), r_dys.tolist(), r_dzs.tolist(), 
           r_ombs.tolist(), r_omas.tolist(), r_dats.tolist(), r_qcs.tolist(),)
 

def main( stime=datetime(2019,9,10,15,0), etime=datetime(2019,9,10,15,0), otyp=4002,
          stime0=datetime(2019,9,10,15,0), TYP="DIST",
          exp="D4_500m_20190910_THIN", tit="", zmin=None, zmax=None ):

   INFO = {
            "EXP":exp, "stime":stime, "stime0":stime0,
            "TOP":"/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"  }

   xs_l = [] 
   ys_l = [] 
   zs_l = [] 
   ombs_l = [] 
   omas_l = [] 
   dats_l = [] 
   qcs_l = [] 

   tcnt = 0
   time = stime
   while( time <= etime):
       INFO["time"] = time
       x_, y_, z_, omb_, oma_, dat_, qc_ = get_dat( INFO, otyp=4001 )

       xs_l.extend(x_)
       ys_l.extend(y_)
       zs_l.extend(z_)

       ombs_l.extend(omb_)
       omas_l.extend(oma_)
       dats_l.extend(dat_)
       qcs_l.extend(qc_)

#       x_, y_, z_, omb_, oma_, dat_, qc_ = get_dat( INFO, otyp=4004 )
#
#       xs_l.extend(x_)
#       ys_l.extend(y_)
#       zs_l.extend(z_)
#
#       ombs_l.extend(omb_)
#       omas_l.extend(oma_)
#       dats_l.extend(dat_)
#       qcs_l.extend(qc_)


       time += timedelta(seconds=30)
       tcnt += 1
   print( len(xs_l) )


   R_diag = desroz_diag_R( ombs_l, omas_l )
   print( "Obs err std:", R_diag )


   ombs_l = np.array( ombs_l )
   omas_l = np.array( omas_l )
   zs_l = np.array( zs_l )

   print( np.min(zs_l), np.max(zs_l), np.size(ombs_l) )
 
   if zmin is not None and zmax is not None: 
      ombs_l = ombs_l[ ( zs_l < zmax ) & ( zs_l >= zmin ) ]
      omas_l = omas_l[ ( zs_l < zmax ) & ( zs_l >= zmin ) ]
      zs_l = zs_l[ ( zs_l < zmax ) & ( zs_l >= zmin ) ]

      tit += " (zmin: {0:.1f}km, zmax: {1:.1f}km)".format( zmin/1000, zmax/1000)

      if np.size(ombs_l) < 1:
         print("No samples", zmin, zmax, np.size(ombs_l))
         sys.exit()
   

   if HIST:
   
      dat_b = ombs_l
      dat_a = omas_l

  
      vmax = 30
      vmin = -30
   
      nbin_b = get_nbin( vmin, vmax, dat_b )
      nbin_a = get_nbin( vmin, vmax, dat_a )
   
      import matplotlib.pyplot as plt
   
      fig, ( (ax) ) = \
      plt.subplots(1, 1, figsize=(8,6))
   
      yy0, xx0, _ = ax.hist( dat_b, range=(vmin,vmax), bins=nbin_b, density=True, 
                             color="k", alpha=0.2, label="O-B, $\sigma$={0:.2f}".format( np.std(dat_b, ddof=1) ) )
      yy1, xx1, _ = ax.hist( dat_a, range=(vmin,vmax), bins=nbin_a, density=True, 
                             color="r", alpha=0.2, label="O-A, $\sigma$={0:.2f}".format( np.std(dat_a, ddof=1) ) )
      ax.legend( fontsize=9 )
   
       
      ymax = yy0.max()
      if ymax < yy1.max():
         ymax = yy1.max()
      ax.vlines( x=np.mean(dat_b), ymin=0.0, ymax=ymax, ls='dashed', color='k', alpha=0.4, lw=2.5 )
      ax.vlines( x=np.mean(dat_a), ymin=0.0, ymax=ymax, ls='dashed', color='r', alpha=0.4, lw=2.5 )
   
#      ax.vlines( x=xx0[np.argmax(yy0)], ymin=0.0, ymax=ymax, ls='dashed', color='k', alpha=0.2 ) # mode

      #ax.vlines( x=0.0, ymin=0.0, ymax=ymax, ls='solid', color='k', alpha=0.4, lw=1.5 )

      ax.set_ylim( 0, ymax )
      ax.set_xlim( vmin, vmax )
 
      from scipy.stats import norm
      mu, std = norm.fit( dat_b )
      x_ = np.linspace( vmin, vmax, 100)
      p = norm.pdf(x_, mu, std)
      plt.plot(x_, p, 'k', linewidth=2, ls='dashed', alpha=0.4 )
  
      ax.text( 0.5, 1.05, "{0:}, N={1:}".format( tit, len(ombs_l) ),
               fontsize=13, transform=ax.transAxes,
               horizontalalignment='center',
               verticalalignment='center', )
   

   
      ofig = "hist_250m_{0:}km_{1:}km.png".format( zmin*0.001, zmax*0.001 )
      print( ofig )
      if quick:
         plt.show()
      else:
         odir = "png"
         os.makedirs( odir, exist_ok=True)
         plt.savefig( os.path.join(odir, ofig),
                      bbox_inches="tight", pad_inches = 0.1)
         plt.clf()
         plt.close('all')
   


##################

otyp = 4002 # vr
#otyp = 4004 # zero Z
otyp = 4001 # Z

exp = "D4_500m_20190910_THIN_NOMEMQC"
exp = "D4_500m_20190910_THIN"
exp = "D4_500m_20190910"
#exp = "D4_500m_20190910_RMEM20"

exp = "D4_250m_20190910_0324_NP0256_NOV_GROSS5_OERR3_HT2_VT2"
tit = "250-m mesh"

#exp = "D4_20190910_TEST0409_bias20"
#tit = "Bias (+20dBZ)"

#exp = "D4_20190910_TEST0409_bias10"
#tit = "Bias (+10dBZ)"

#exp = "D4_20190910_TEST0409_bias00"
#tit = "Bias (+0dBZ)"

stime0 = datetime(2019,9,10,12,0,0)
stime = datetime(2019,9,10,12,10,30)
#stime = datetime(2019,9,10,12,15,30)
etime = datetime(2019,9,10,12,20,0)

stime = datetime(2019,9,10,12,20,30)
etime = datetime(2019,9,10,12,25,0)
#etime = stime

#stime0 = datetime(2019,9,10,15,0,0)
#stime = datetime(2019,9,10,15,5,30)
#etime = datetime(2019,9,10,15,9,0)

zmin = 0
zmin = 2000
zmin = 4000
zmin = 6000
zmin = 8000
zmin = 10000
zmax = zmin + 2000


main( stime=stime, etime=etime, stime0=stime0, otyp=otyp, exp=exp, tit=tit,
      zmin=zmin, zmax=zmax )  
