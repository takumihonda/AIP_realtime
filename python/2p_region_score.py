import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

from tools_AIP import read_obs, read_mask, read_nc_lonlat, read_fcst_grads, read_obs_grads_latlon, read_fcst_grads_all

quick = False
quick = True

###########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"




EXP = "D4_500m_CTRL"

theight = 3000.0
#theight = 6000.0


#etime = stime
stime = datetime( 2019, 8, 24, 15, 0, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

stime = datetime( 2019, 8, 24, 15, 20, 0 )
etime = datetime( 2019, 8, 24, 15, 30, 0 )

#stime = datetime( 2019, 8, 19, 13, 0, 30 )
#etime = datetime( 2019, 8, 19, 14, 0, 0 )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )
#time0 = datetime( 2019, 8, 19, 13, 0, 0 )

tmin = 0
tmax = 61 # max time dimension does not include FT=0

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/20201117/{0:}/dafcst".format( EXP )
#OBS_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/{1:}/pawr_grads".format( EXP )

fcst_zmax = 43

obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )
mask = read_mask()


INFO = {
         "mask": mask,
         "EXP": EXP, 
#         "OBS_DIR": OBS_DIR,
         "FCST_DIR": FCST_DIR,
         "obsz": obsz,
         "olon2d": olon2d,
         "olat2d": olat2d,
         "lon2d": lon2d,
         "lat2d": lat2d,
         "hgt3d": hgt3d,
         "ohgt3d": ohgt3d,
         "cz": cz,
         "gz": fcst_zmax,
         "gy": lon2d.shape[0],
         "gx": lon2d.shape[1],
       }



dbz_thrs_l = [ 15.0, 30.0, ]

tskip = 1

itmax = int( ( etime - stime ).total_seconds() / 30 ) + 1

tlevs = np.arange( tmin, tmax, tskip, dtype=np.int32 )
tlev_max = np.shape( tlevs )[0]

ts_l = np.zeros( ( tlev_max, len( dbz_thrs_l) ) )
bs_l = np.zeros( ( tlev_max, len( dbz_thrs_l) ) )

ts_l[:] = np.nan
bs_l[:] = np.nan

odir = "ts_npz/" + INFO["EXP"]
os.makedirs( odir, exist_ok=True)




lats = 36.05
late = 36.1

lons = 139.7
lone = 139.8

vtime = datetime( 2019, 8, 24, 15, 30, 0 )

ts_l = []
wmax_l = []

w3d = read_fcst_grads_all( INFO, itime=stime, tlev=0 , FT0=True, nvar="w" ) 
lon3d = np.resize( lon2d, w3d.shape )
lat3d = np.resize( lat2d, w3d.shape )

times = []

time = stime
while (time <= etime):
  print( "Initial time:", time )

  tlev = int( ( vtime - time ).total_seconds() / 30 )
  print( tlev, time, vtime)

  dbz = 30.0
#  fn_ts = "TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}_lon{3:.2f}-{4:.2f}_lat{5:.2f}-{6:.2f}.npz".format( dbz, theight, time.strftime('%H%M%S_%Y%m%d'), lons, lone, lats, late )
  fn_ts = "TS_3d_thrs{0:.1f}dbz_i{1:}_lon{2:.2f}-{3:.2f}_lat{4:.2f}-{5:.2f}.npz".format( dbz, time.strftime('%H%M%S_%Y%m%d'), lons, lone, lats, late )

  fn = os.path.join( odir, fn_ts )
  try:
    dat = np.load( fn )[ 'ts' ]
    ts_l.append( dat[tlev] )
  except:
    print( 'Failed to read. Skip')
    sys.exit()
    time += timedelta(seconds=30)
    continue


  w3d = read_fcst_grads_all( INFO, itime=time, tlev=tlev , FT0=True, nvar="w" ) 
  w3d_ = w3d[ ( lon3d >= lons ) & ( lon3d <= lone ) &
              ( lat3d >= lats ) & ( lat3d <= late ) ]
  wmax_l.append( np.nanmax( w3d_ ) )

  times.append( time )

#  print( dat[tlev] )
  
  time += timedelta(seconds=30)


import matplotlib.pyplot as plt
import matplotlib.dates as mdates

fig, (ax1) = plt.subplots(1, 1, figsize=(7,5) )
fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, )


ax2 = ax1.twinx()

ax_l = [ ax1, ax2 ]
data_l = [ ts_l, wmax_l ]

lab_l = [ 'Threat score', 'Wmax' ]
color_l = [ 'k', 'magenta' ]


for i, ax in enumerate( ax_l ):
   ax.xaxis.set_major_locator(mdates.SecondLocator(interval=120) )
   ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

   ax.set_xlim( stime, etime )
   ax.plot( times, data_l[i], marker='o', color=color_l[i], 
            ms=4.0, label=lab_l[i] )

ax1.set_xlabel( 'Forecast initial time (UTC)', fontsize=12 )
ax1.set_ylim( 0.0, 0.99 )
ax2.set_ylim( 2, 14 )

ax1.vlines( x=datetime( 2019, 8, 24, 15, 25, 0) , ymin=0.0, ymax=0.99, colors='gray', ls='dashed',
        lw=0.5)

ax1.annotate('', xy=[ datetime( 2019, 8, 24, 15, 25, 0), 0.1], 
                 xytext=[ datetime( 2019, 8, 24, 15, 30, 0), 0.1], 
            arrowprops=dict( arrowstyle='<->',
                             #shrink=0, width=1, headwidth=8, 
                             #headlength=10, 
                             connectionstyle='arc3',
                             facecolor='gray', edgecolor='gray')
           )
ax1.text( datetime( 2019, 8, 24, 15, 27, 30), 0.11, 'JMA nowcast\n updated every 5 min',
         fontsize=12, #transform=ax.transAxes,
         ha='center',
         va='bottom' )

ax1.set_ylabel( 'Threat score', fontsize=12 )
ax2.set_ylabel( r'Wmax (m s$^{-1}$)', fontsize=12 )
ax2.yaxis.label.set_color( color_l[1] )

ax1.legend( loc='upper left', fontsize=12 )
ax2.legend( loc='upper left', fontsize=12, bbox_to_anchor=(0, 0.9), )
ax2.tick_params(axis='y', colors=color_l[1] )
fig.suptitle( 'Threat score and Wmax\nvalid at {0:}'.format( vtime.strftime('%H:%M:%S %m/%d/%Y') ), fontsize=13)

ofig = "2p_scores.png"
print(ofig)

if not quick:
   opath = "png"
   ofig = os.path.join(opath, ofig)
   plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
   print(ofig)
   plt.clf()
else:
   plt.show()



#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0,4) )
#fig.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9, )
#
#ax_l = [ ax1, ax2 ]
#
#tit_l = [ 'Threat score', 'Wmax']
#
#ax1.plot( times, ts_l, marker='o', color='k', ms=4.0 )
#ax2.plot( times, wmax_l, marker='o', color='k', ms=4.0 )
#
#for i, ax in enumerate( ax_l ):
#   ax.xaxis.set_major_locator(mdates.SecondLocator(interval=120) )
#   ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
#
#   ax.set_xlim( stime, etime )
#
#   ax.text( 0.5, 1.01, tit_l[i],
#            fontsize=13, transform=ax.transAxes,
#            ha='center',
#            va='bottom' )
#
#   if i == 0:
#      ax.text( 0.9, -0.1, 'Forecast initial time (UTC)',
#               fontsize=12, transform=ax.transAxes,
#               ha='left',
#               va='top' )

   

