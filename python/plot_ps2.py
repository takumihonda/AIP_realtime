import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta





def read_fcst_grads_ps( INFO, itime=datetime(2019,9,3,2,0,0)  ):

    #fn = os.path.join( "/data15/honda/SCALE-LETKF/AIP/verify/d4_500m_small",
    fn = os.path.join( INFO["TOP"], INFO["EXP"],
                       INFO["time0"].strftime('%Y%m%d%H0000'),  "dafcst/fcst_all_" +
                       itime.strftime('%Y%m%d-%H%M%S.grd') )

    try:
       infile = open(fn)
    except: 
       print("Failed to open")
       print( fn )
       sys.exit()
    
    gx = 256
    gy = 256
    gz = 1
    rec3d = gx*gy*gz*61 # 61 tlevs, 2D data

    nv = 1
 
#    # grads file starts from FT > 0s
#    gtlev = tlev - 1
#    rec = rec3d * nv * gtlev

    infile.seek( 0 )
    tmp3d = np.fromfile(infile, dtype=np.dtype('>f4'), count=rec3d)  # big endian   
    input3d = np.reshape( tmp3d, (61, gz,gy,gx) )

    return( input3d ) 


#############

def main( INFO, itime=datetime(2019,9,3,2,0), nlev=1 ):

  
    fcst4d = read_fcst_grads_ps( INFO, itime=itime  )
   
    dt = 30*nlev
    dpdt = np.mean( np.abs( np.diff( fcst4d, n=nlev, axis=0 ) / dt ), axis=(1,2,3) )
    #print( fcst4d.shape, dpdt.shape )

    return( dpdt )


###################

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"

EXP1 = "D4_500m_TEST_DEFAULT_0515_MEAN"

LAB1 = "Default" 

EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG02VG02_M2"
LAB2 = "H2V2"

EXP3 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG04_M2"
LAB3 = "H4V4"
EXP4 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG06VG06_M2"
LAB4 = "H6V6"
EXP5 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG08VG08_M2"
LAB5 = "H8V8"

stime = datetime( 2019, 8, 24, 15, 0, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )
#etime = datetime( 2019, 8, 24, 15, 1, 30 )

#etime = stime

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )

tmin = 1
tmax = 61 # max time dimension does not include FT=0

INFO = { "TOP": TOP,
         "EXP": EXP1,
         "time0": time0,
       }


tskip = 1


cnt = 0

nlev = 2 # 2nd derivative

time = stime
while (time <= etime):
  cnt += 1
  print( "Initial time:", time )


  INFO["EXP"] = EXP1
  dpdt1_ = main( INFO, itime=time, nlev=nlev )

  INFO["EXP"] = EXP2
  dpdt2_ = main( INFO, itime=time, nlev=nlev )

  INFO["EXP"] = EXP3
  dpdt3_ = main( INFO, itime=time, nlev=nlev )

  INFO["EXP"] = EXP4
  dpdt4_ = main( INFO, itime=time, nlev=nlev )

  INFO["EXP"] = EXP5
  dpdt5_ = main( INFO, itime=time, nlev=nlev )


  if cnt == 1:
     dpdt1 = np.copy( dpdt1_ )
     dpdt2 = np.copy( dpdt2_ )
     dpdt3 = np.copy( dpdt3_ )
     dpdt4 = np.copy( dpdt4_ )
     dpdt5 = np.copy( dpdt5_ )
  else: 
     dpdt1 += dpdt1_
     dpdt2 += dpdt2_
     dpdt3 += dpdt3_
     dpdt4 += dpdt4_
     dpdt5 += dpdt5_

 #     bs_l.append( bs )
 #     ftime_l.append( time + timedelta(seconds=tlev*30) )
  
#  odir = "ts_npz/" + INFO["EXP"]
#  os.makedirs( odir, exist_ok=True)
#  np.savez(os.path.join(odir,fn_ts), ts=np.array(ts_l), bs=np.array(bs_l), times=ftime_l )
#
  time += timedelta(seconds=30)

dpdt1 = dpdt1 / cnt
dpdt2 = dpdt2 / cnt
dpdt3 = dpdt3 / cnt
dpdt4 = dpdt4 / cnt
dpdt5 = dpdt5 / cnt

print( dpdt1.shape)
print( dpdt2.shape)

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm

if nlev == 2:
   t_l = np.arange( 30, 30*60, 30 )
   tit =  r'Domain averaged |$\partial^2 p_s/\partial$t$^2$| (Pa s$^{-2}$)'
elif nlev == 1:
   t_l = np.arange( 30, 30*61, 30 )
   tit =  r'Domain averaged |$\partial p_s/\partial$t| (Pa s$^{-1}$)'

print(t_l.shape)


fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, )

dat = dpdt1

dat_l = [ dpdt1, dpdt2, dpdt3, dpdt4, dpdt5,  ]
c_l = [ 'k', 'r', 'b', 'g', 'y']
lab_l = [ LAB1, LAB2, LAB3, LAB4, LAB5, ]

xmax = 1800
xmin = 0

ymin = 0
ymax = 3

for i, dat in enumerate( dat_l):
   ax1.plot(t_l, dat, lw=1.0, color=c_l[i], linestyle='solid',
             label=lab_l[i] )


ax1.set_xlim( xmin, xmax )
ax1.set_ylim( ymin, ymax )
ax1.set_xlabel( "Forecast time (s)", fontsize=14 )

dxsec = 60
xsecs = np.arange(dxsec, 1800, dxsec)
ax1.vlines( x=xsecs, ymin=0, ymax=ymax, ls='dashed', lw=0.5)

ax1.set_xticks( np.arange( 0, 2100, 300) )


ax1.legend( fontsize=14, loc='upper right')


ax1.text( 0.5, 1.08, tit, 
          fontsize=16, transform=ax1.transAxes,
          horizontalalignment='center',
          verticalalignment='top',
          )

plt.show() 

#print( ts_l )
sys.exit()

