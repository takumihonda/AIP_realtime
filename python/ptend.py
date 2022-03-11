import os
import sys
import numpy as np
from datetime import datetime, timedelta
from tools_AIP import read_obs_grads, prep_proj_multi, read_nc_topo, read_mask_full, read_obs_grads_latlon, read_fcst_grads_all, read_nc_lonlat, dist

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm

from scipy.interpolate import griddata



quick = True
#quick = False

USE_ARCH_DAT = True
#USE_ARCH_DAT = False

def main( INFO, itime_l=[], EXP_l=[], tit_l=[] ):

    data_path = "../../dat4figs_JAMES/Fig18"
    ofig = "Fig18.pdf"
    os.makedirs( data_path, exist_ok=True )
    fn = '{0:}/data.npz'.format( data_path, )


    nd = 2 # second derivative

    ps3d = np.zeros( ( INFO["TMAX"], INFO["gy"], INFO["gx"] ) )
    ptend_l = np.zeros( ( len( EXP_l ), INFO["TMAX"]-nd ) )
    

    if not USE_ARCH_DAT:

       for i, EXP_ in enumerate( EXP_l ):
   
          INFO["FCST_DIR"] = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/dafcst".format( EXP_ )
          for itime in itime_l: 
             print( "initial", itime )
             dir_in = os.path.join( "dat_ptend", EXP_, )
             os.makedirs( dir_in, exist_ok=True)
             ofile = os.path.join( dir_in, "ptend2_abs_{0:}.npz".format(  itime.strftime('%Y%m%d%H%M%S') ) )
      
             try:
                data = np.load( ofile )
                ptend_ = data["ptend"]
             except:
                print( "No npz file ", ofile)
      
                for tlev in range( INFO["TMAX"] ):
                   prs3d_ = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar="p" )
                   ps3d[tlev,:,:] = prs3d_[0,:,:]
         
                ptend_ = np.average( np.abs( np.diff( ps3d, axis=0, n=nd ), ), axis=(1,2) ) / ( INFO["DT"]**2 )
                np.savez( ofile, ptend=ptend_ )
      
             ptend_l[i,:] += ptend_
   
       ptend_l = ptend_l / len( itime_l )
 
       np.savez( fn, ptend_l=ptend_l )   

    else:

       ptend_l = np.load( fn )['ptend_l']

    fig, (ax1) = plt.subplots(1, 1, figsize= (6,4 ))
    fig.subplots_adjust(left=0.15, bottom=0.12, right=0.98, top=0.9, )


    xlab = "Forecast time (min)"
    #ylab = r'Second derivative of the lowest-level pressure ($\partial^2/\partial t^2, Pa s^{-2}$)'
    if nd == 2:
       ylab = r'$\partial^2p/\partial t^2$ (Pa s$^{-2})$'
       fts = np.arange( len( ptend_l[0,:] ) ) * 0.5 + 0.5 # every 30 seconds
    print( fts )

    xmin = 0
    xmax = 30

    ymin = 0.0
    ymax = 0.2

    ax1.set_xlabel( xlab, fontsize=11 )
    ax1.set_ylabel( ylab, fontsize=11 )
    ax1.set_xlim( xmin, xmax )
    ax1.set_ylim( ymin, ymax )
    ax1.grid( lw=0.5, ls='dashed' )


    dy = 0.05
    ylabs = np.arange( ymin, ymax+dy, dy )
    ax1.set_yticks( ylabs )

    print( "data plot")
    c_l = [ 'r', 'k', 'b' ]

    
 
    for i in range( len( EXP_l) ):
       ax1.plot( fts, ptend_l[i,:], c=c_l[i],
                 label=tit_l[i] )

    
    ax1.legend()

    tit = "Imbalance measured by the second time derivative of\nthe lowest-level pressure"

    ax1.text( 0.5, 1.01, tit,
              fontsize=12, transform=ax1.transAxes,
              ha='center',
              va='bottom',
            )


    opath = "pdf"
    #ofig = "1p_dpdt2.png"
    print(ofig)

    if not quick:
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       plt.show()



############

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"
EXP = "20201117/D4_500m_CTRL"

tit_l = [ "H1V1",
          "H4V4 (CTRL)",
          "H8V8",
        ]

EXP_l = [ 
          "20201117/D4_500m_H1V1",
          "20201117/D4_500m_CTRL",
          "20201117/D4_500m_H8V8",
#          "20201117/D4_500m_H4V1",
#          "20201117/D4_500m_H8V1",
        ]

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/dafcst".format( EXP )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

stime = datetime( 2019, 8, 24, 15, 0, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

#time0 = datetime( 2019, 8, 19, 13, 0, 0 )
#stime = datetime( 2019, 8, 19, 13, 0, 30 )
#etime = datetime( 2019, 8, 19, 14, 0, 0 )

fcst_zmax = 43

obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )

TMAX = 61
DT = 30.0
INFO = { "TOP": TOP,
         "EXP": EXP,
         "time0": time0,
         "FCST_DIR": FCST_DIR,
         "gz": fcst_zmax,
         "gy": lon2d.shape[0],
         "gx": lon2d.shape[1],
         "lon2d": lon2d,
         "lat2d": lat2d,
         "cz": cz,
         "DT": DT,
         "TMAX": TMAX,
       }




itime_l = []
time = stime
while time <= etime:
   itime_l.append( time )
   time += timedelta( seconds=30 )

main( INFO, itime_l=itime_l, EXP_l=EXP_l, tit_l=tit_l )

