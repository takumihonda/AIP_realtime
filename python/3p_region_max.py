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

theight = 2000.0


fstime = datetime( 2019, 8, 24, 15, 20, 0 )
fetime = datetime( 2019, 8, 24, 15, 30, 0 )
fitmax = int( ( fetime - fstime ).total_seconds() / 30 ) + 1

nstime = datetime( 2019, 8, 24, 15, 20, 0 )
netime = datetime( 2019, 8, 24, 15, 35, 0 )

stime = datetime( 2019, 8, 24, 15, 20, 0 )
etime = datetime( 2019, 8, 24, 15, 40, 0 )


stime_ = datetime( 2019, 8, 24, 15, 20, 0 )
etime_ = datetime( 2019, 8, 24, 15, 40, 0 )
#etime_ = datetime( 2019, 8, 24, 15, 50, 0 )






odir = "rainmax_npz/" + EXP
os.makedirs( odir, exist_ok=True )






lats = 36.05
late = 36.1
#late = 36.15

lons = 139.7
lone = 139.8


ms = 10.0 
mt = 'o'
ls_jma = 'dashed'
lw_jma = 2.0
cjma = 'gray'

lw_obs = 2.0
ls_obs = 'solid'
cobs = 'k'

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm


fig, (ax1, ax2, ax3 ) = plt.subplots(1, 3, figsize=(14,4.5) )
fig.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.92, )

ax_l = [ ax1, ax2, ax3 ]

# obs
#ax1.plot( otime_l, obs_l, color=cobs, linewidth=lw_obs, linestyle=ls_obs )

fac_l = [ 1.0, 1.e0, 1.0 ]

# fcst
time = fstime
it = 0
while (time <= fetime):
   fn_fcst = "fcst_max_lon{0:.2f}-{1:.2f}_lat{2:.2f}-{3:.2f}_i{4:}.npz".format( lons, lone, lats, late, time.strftime('%Y%m%d%H%M%S') )
   try: 
     data = np.load( os.path.join( odir, fn_fcst ), allow_pickle=True  )
     w_l = data['wmax']
     #qh_l = data['qhmax']
     water_l = data['awater']
     pw_l = data['apw']
     ftime_l = data['times']

     fcolor = cm.jet( it/fitmax ) 
     lw = 2.0

     data_l = [ w_l, water_l, pw_l ]
     for j, data in enumerate( data_l ):
         ax_l[j].plot( ftime_l, data*fac_l[j], color=fcolor, lw=lw )

         ax_l[j].plot( ftime_l[0], data[0]*fac_l[j], color=fcolor, lw=lw,
                   marker=mt, markeredgecolor='w', ms=ms )


   except:
     print( 'skip forecast', fn_fcst )

   it += 1
   time += timedelta( seconds=30 )

## nowcast
#time = nstime
#while (time <= netime):
#   fn_nowcast = "nowcast_rainmax_z{0:.1f}_lon{1:.2f}-{2:.2f}_lat{3:.2f}-{4:.2f}_i{5:}.npz".format( theight, lons, lone, lats, late, time.strftime('%Y%m%d%H%M%S') )
#   try: 
#     data = np.load( os.path.join( odir, fn_nowcast ), allow_pickle=True  )
#     jma_l = data['rmax']
#     jtime_l = data['times']
#     ax1.plot( jtime_l, jma_l, color=cjma, linewidth=lw_jma, 
#               linestyle=ls_jma )
#
#     ax1.plot( jtime_l[0], jma_l[0], color=cjma, lw=lw,
#               marker=mt, markeredgecolor='w', ms=ms )
#
##     ax1.vlines( x=jtime_l[0], ymin=ymax-20, ymax=ymax, color='gray', 
##                 lw=1.0, linestyle='solid',
##               )     
#
#   except:
#     print( 'skip nowcast', time )
#
#   time += timedelta( seconds=30 )


## legend fcst
#
#y_obs = ymax   - 10
#y_scale = ymax - 18
#y_jma = ymax   - 26
#
#xmin0_ = stime_ + timedelta( seconds=60 )
#dx_ = timedelta( seconds=5 )
#xmax0_ = xmin0_ + dx_*fitmax
#
#xmin_ = xmin0_
#for it in range( fitmax ):
#     fcolor = cm.jet( it/fitmax ) 
#     ax1.plot( [ xmin_, xmin_+dx_], [ y_scale, y_scale ], 
#         color=fcolor, lw=lw, )
#     xmin_ += dx_
#
#ax1.text( xmax0_+dx_*2, y_scale, 'SCALE-LETKF',
#          color='k',
#          fontsize=12, transform=ax1.transData,
#          ha='left',
#          va='center' )
#
## legend obs
#ax1.plot( [ xmin0_, xmax0_], [y_obs, y_obs], color=cobs, linewidth=lw_obs, linestyle=ls_obs )
#ax1.text( xmax0_+dx_*2, y_obs, 'MP-PAWR obs',
#          color='k',
#          fontsize=12, transform=ax1.transData,
#          ha='left',
#          va='center' )
#
## legend nowcast
#ax1.plot( [ xmin0_, xmax0_], [y_jma, y_jma], color=cjma, linewidth=lw_jma, 
#          linestyle=ls_jma )
#
#ax1.text( xmax0_+dx_*2, y_jma, 'JMA nowcast',
#          color='k',
#          fontsize=12, transform=ax1.transData,
#          ha='left',
#          va='center' )
#


#---

tit_l = [
        r'Maximum upward velocity (m s$^{-1}$)', 
        #r'Area-averaged total water path (kg m$^{-2}$)', 
        'Area-averaged total (liquid + ice)\n' + r'water path (kg m$^{-2}$)', 
        r'Area-averaged PW (mm m$^{-2}$)', 
        ]

#ax1.text( 1.0, 1.01, 'Z={0:.0f} km'.format( theight/1000.0 ),
#          color='k',
#          fontsize=11, transform=ax1.transAxes,
#          ha='right',
#          va='bottom' )


ymax = 120
ymin = 0

ymax_l = [ 16, 4, 54 ]
ymin_l = [ 0, 0, 48 ]

bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
         'edgecolor':'w' }

pnum_l = [ '(a)', '(b)', '(c)' ]

for j, ax in enumerate( ax_l ):
    ax.xaxis.set_major_locator( mdates.MinuteLocator(interval=5) )
    ax.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S') )
    #ax.xaxis.set_major_formatter( mdates.DateFormatter('%H%M') )
    
    ax.set_xlim( stime_, etime_ )
    ax.set_ylim( ymin_l[j], ymax_l[j] )

    dx = 0.0
    if j == 1:
      dx = 0.05
    ax.text( 0.5+dx, 0.99, tit_l[j],
             color='k',
             fontsize=12, transform=ax.transAxes,
             ha='center',
             va='top', 
             bbox=bbox )

    ax.text( 0.01, 0.99, pnum_l[j],
             color='k',
             fontsize=12, transform=ax.transAxes,
             ha='left',
             va='top', )



ax2.set_xlabel( 'Time (UTC)', fontsize=12 )
#ax1.set_ylabel( r'Maximum rainfall intensity (mm h$^{-1}$)', fontsize=12 )

tit = 'SCALE-LETKF forecast'
fig.suptitle( tit, fontsize=14)

ofig = "3p_max_lon{0:.2f}-{1:.2f}_lat{2:.2f}-{3:.2f}.png".format( lons, lone, lats, late, )
print(ofig)

if not quick:
   opath = "png"
   ofig = os.path.join(opath, ofig)
   plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
   print(ofig)
   plt.clf()
else:
   plt.show()



sys.exit()



