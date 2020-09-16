import numpy as np
import sys
import os

from datetime import datetime
from datetime import timedelta

from tools_AIP import read_dat

def take_ave2d( rz2d, rz2d_cnt, INFO, otyp=4002, TYP="DIST", REF="all" ):

   tfoot = INFO["time"].strftime('%Y%m%d-%H%M%S')
   fn = os.path.join( INFO["TOP"], INFO["EXP"], INFO["stime0"].strftime('%Y%m%d%H%M%S'),
                      "obsdep",
                      "obsdep_" + tfoot + ".dat")
#   fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/D4_500m_20190910_THIN/20190910150000/obsdep/obsdep_" + tfoot + ".dat"
   
   print("Get data:",tfoot)
   elms, lons, lats, levs, dats, ombs, qcs, _ = read_dat( fn )
 
  
   # radar reflectivity
   r_dxs = lons[ elms == otyp ]
   r_dys = lats[ elms == otyp ]
   r_dzs = levs[ elms == otyp ]
   r_ombs = ombs[ elms == otyp ]
   r_dats = dats[ elms == otyp ]
   r_qcs = qcs[ elms == otyp ]

   print( "Num:", len(r_dats), np.min(r_dats), np.max(r_dats) )

   if REF == "all" and ( otyp == 4001 or otyp == 4004):
      r_dxs = lons[ (elms == 4001) | (elms == 4004) ]
      r_dys = lats[ (elms == 4001) | (elms == 4004) ]
      r_dzs = levs[ (elms == 4001) | (elms == 4004) ]
      r_ombs = ombs[ (elms == 4001) | (elms == 4004) ]
      r_dats = dats[ (elms == 4001) | (elms == 4004) ]
      r_qcs = qcs[ (elms == 4001) | (elms == 4004) ]


   hdist = np.sqrt( np.square( r_dxs ) + np.square( r_dys ) )

   zidxs = np.rint( r_dzs[:]/INFO["DZ"] ).astype(np.int32) 

   if TYP is "DIST":
      ridxs = np.rint( ( hdist[:] - np.min(INFO["r1d"]) ) / INFO["DR"] ).astype(np.int32) 
   elif TYP is "RAW":
      ridxs = np.rint( ( r_dats[:] - np.min(INFO["r1d"]) ) / INFO["DR"] ).astype(np.int32) 


   for z in range(INFO["nz"]):
       for r in range(INFO["nr"]):
           rz2d[r,z] += np.sum( r_ombs[ (ridxs == r) & (zidxs == z)] )

           if len( r_ombs[ (ridxs == r) & (zidxs == z)] ) > 0:
              rz2d_cnt[r,z] += len( r_ombs[ (ridxs == r) & (zidxs == z) ] )

   return( rz2d, rz2d_cnt )


def main( stime=datetime(2019,9,10,15,0), etime=datetime(2019,9,10,15,0), otyp=4002,
          stime0=datetime(2019,9,10,15,0), TYP="DIST",
          exp="D4_500m_20190910_THIN", REF="all" ):

   DR = 500.0
   DZ = 500.0

   if TYP is "DIST":
      r1d = np.arange(0, 60.0e+3+DR, DR)
      zfac = 0.001
      rfac = 0.001
   elif TYP is "RAW":
      DR = 5.0
      r1d = np.arange(0, 70.0+DR, DR)
      if otyp == 4002:
         r1d = np.arange(-DR*10, DR*11, DR)
      zfac = 0.001
      rfac = 1.0
   nr = len(r1d)
   z1d = np.arange(0, 12.0e+3+DZ, DZ)
   nz = len(z1d)

   rz2d = np.zeros( (nr,nz) )
   rz2d_cnt = np.zeros( (nr,nz) )

   INFO = { "nr":nr, "nz":nz, "DR":DR, "DZ":DZ, 
            "z1d":z1d, "r1d":r1d, "EXP":exp, "stime":stime, "stime0":stime0,
            "TOP":"/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"  }


   tcnt = 0
   time = stime
   while( time <= etime):
       INFO["time"] = time
       rz2d, rz2d_cnt = take_ave2d( rz2d, rz2d_cnt, INFO, otyp=otyp, TYP=TYP, REF=REF )

       time += timedelta(seconds=30)
       tcnt += 1
   import matplotlib.pyplot as plt

   cmap = plt.cm.get_cmap( "RdBu_r" )
   levs = np.arange(-5,5.5,0.5)
   #levs = np.arange(-10,10.5,0.5)
   #levs = np.arange(-3,3.5,0.5)
   cmap.set_under('gray', alpha=1.0)
   cmap.set_over('k', alpha=1.0)

   y2d, x2d = np.meshgrid( z1d*zfac, (r1d + DR*0.5)*rfac )

   fig, ( (ax) ) = \
   plt.subplots(1, 1, figsize=(8,6))

   fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.92,
                       wspace=0.2, hspace=0.2)


   rz2d = rz2d / rz2d_cnt
   rz2d[ rz2d_cnt < tcnt*10 ] = np.nan
   print(x2d.shape, y2d.shape, rz2d.shape, np.max(rz2d_cnt))
   SHADE = ax.pcolormesh( x2d, y2d, rz2d,
                          cmap=cmap, 
                          vmax=np.max(levs),
                          vmin=np.min(levs),
                          )
 
   pos = ax.get_position()
   cb_h = pos.height
   cb_w = 0.01
   ax_cb = fig.add_axes( [pos.x1+0.01, pos.y0, cb_w, cb_h] )
   cb = plt.colorbar( SHADE, cax=ax_cb, orientation = 'vertical',
                      ticks=levs[::2], extend='both' )

   if otyp == 4001: 
      tvar = "radar reflectivity"
   elif otyp == 4004: 
      tvar = "\"zero\" radar reflectivity"
   elif otyp == 4002: 
      tvar = "radial velocity"
   else:
      tvar = ""

   if TYP is "DIST":
      xlabel = "Distance from Saitama PAWR (km)"
   elif TYP is "RAW":
      xlabel = "Observed " + tvar + " (dBZ)"
   ylabel = "Height relative to Saitama PAWR (km)"

   ax.set_xlabel( xlabel, fontsize=14 )
   ax.set_ylabel( ylabel, fontsize=14 )


   ax.text( 0.5, 1.05, "Mean of " + tvar + " O-B",
           fontsize=16, transform=ax.transAxes,
           horizontalalignment='center',
           verticalalignment='top', 
           )


   plt.show()
   sys.exit()

   print(rz2d.shape)
   print(hdist.shape, np.min(hdist), np.max(hdist))
#   print(r_lons)
   print(elms.shape)

   print( np.min(dats) )

##################

otyp = 4002 # vr
#otyp = 4004 # zero Z
#otyp = 4001 # Z

exp = "D4_500m_20190910_THIN_NOMEMQC"
exp = "D4_500m_20190910_THIN"
exp = "D4_500m_20190910"
#exp = "D4_500m_20190910_RMEM20"

stime = datetime(2019,9,10,15,5,0)
etime = datetime(2019,9,10,15,10,30)

stime0 = datetime(2019,9,10,9,0,0)

stime = datetime(2019,9,10,9,10,0)
etime = datetime(2019,9,10,9,40,0)
#etime = datetime(2019,9,10,9,10,0)

stime = datetime(2019,9,10,9,10,0)
etime = datetime(2019,9,10,9,40,0)
etime = datetime(2019,9,10,9,10,0)

etime = stime # DEBUG

TYP = "DIST"
#TYP = "RAW"

REF = "all"
#REF = "only"

main( stime=stime, etime=etime, stime0=stime0, otyp=otyp, TYP=TYP, exp=exp, REF=REF )  
