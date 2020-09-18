import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta


quick = True
quick = False

#theight = 3000.0
theight = 2000.0

# UTC
    
stime = datetime( 2020, 8, 24, 7, 0 )
etime = datetime( 2020, 9, 7,  0, 0 )

# Figure x range
stime_ = datetime( 2020, 8, 25, 0, 0 )
etime_ = datetime( 2020, 9, 7,  0, 0 )

tmin = 0
tmax = 61 # max time dimension includes FT=0
tskip = 20
    
itmax = int( ( etime - stime ).total_seconds() / 30 ) + 1

tlevs = np.arange( tmin, tmax, tskip, dtype=np.int32 )
tlev_max = np.shape( tlevs )[0]
       


dbz = 15.0
#dbz = 30.0

###########


def read_scores():

    odir = "ts_npz/realtime_score"
    
    fn_ts_t = os.path.join( odir, "TS_thrs{0:.1f}dbz_z{1:.1f}_tskip{2:0=3}_s{3:}_e{4:}.npz".format( dbz, theight, tskip, stime.strftime('%H%M%S_%Y%m%d'), etime.strftime('%H%M%S_%Y%m%d') ) )
    
    try:
       print( "Read ", fn_ts_t )
       data = np.load( fn_ts_t, allow_pickle=True )
       return( data['ts_l'], data['bs_l'], data['itimes_l'] )
    except:
       print( "Read all files")

       ts_l = np.zeros( ( itmax, tlev_max,  ) )
       bs_l = np.zeros( ( itmax, tlev_max,  ) )
       
       ts_l[:] = np.nan
       bs_l[:] = np.nan
       
       
       itimes_l = []
       
       time = stime
       cnt = 0
       while (time <= etime):
#          print( "Initial time:", time )
        
          itimes_l.append( time + timedelta(hours=9))
          
          fn_ts = "TS_thrs{0:.1f}dbz_z{1:.1f}_tskip{2:0=3}_i{3:}.npz".format( dbz, theight, tskip, time.strftime('%H%M%S_%Y%m%d') )
        
          try:
             data = np.load( os.path.join(odir,fn_ts) )
             ts_l[cnt,:] = data['ts_l']
             bs_l[cnt,:] = data['bs_l']
          except:
             print( "failed to read", time )
        
          time += timedelta(seconds=30)
          cnt += 1

       np.savez( fn_ts_t, ts_l=ts_l, bs_l=bs_l, itimes_l=itimes_l )
       return( ts_l, bs_l, itimes_l )


ts_l, bs_l, itimes_l =  read_scores()

# get lead times for nan
path = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825"
ofn = os.path.join( path, "result_leadtime/leadtime.npz" )
data = np.load( ofn, allow_pickle=True )
lt_l = data['lt_l']
ftime_l = data['ftime_l']

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ymin1 = 0.0
ymax1 = 1.0
ymin2 = 0.0
ymax2 = 5.0

print( tlevs.shape)
print( tlevs )

fig, (ax1,ax2) = plt.subplots( 2, 1, figsize=( 16,8 ) )
fig.subplots_adjust(left=0.05, bottom=0.07, right=0.98, top=0.95,
                    hspace=0.3 )

cc_l = [ "k", "r", "b", "g" ]
print( ts_l.shape )


lw = 0.3
alp = 0.5
for i, ft in enumerate( tlevs ):
    ft_ = int( ft * 30 / 60 ) # min

    itimes_l_ = itimes_l + timedelta( minutes=ft_ )
    

    ax1.plot( itimes_l_, ts_l[:,i], color=cc_l[i], label="FT={0:.0f}min".format(ft_), 
              lw=lw )
    ax2.plot( itimes_l_, bs_l[:,i], color=cc_l[i], lw=lw )

ax1.vlines( x=ftime_l[ np.isnan( lt_l ) ], ymin=ymin1, ymax=ymax1, color='gray', 
            ls='dashed', lw=0.01, alpha=alp )
ax2.vlines( x=ftime_l[ np.isnan( lt_l ) ], ymin=ymin2, ymax=ymax2, color='gray', 
            ls='dashed', lw=0.01, alpha=alp )


leg = ax1.legend( loc='lower left', fontsize=12, framealpha=1.0,
                  )

for line, text in zip(leg.get_lines(), leg.get_texts()):
    text.set_color(line.get_color())

#for i, text in enumerate( leg.get_texts() ):
#    ax1.setp( text, color=cc_l[i] )

bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':3 }

dmin = 120
dmin = 360
dmin = 720

tit_l = [ "Threat score", "Bias score" ]
ax_l = [ ax1, ax2 ]
for i, ax in enumerate( ax_l ):
   ax.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M\n%m/%d') )
   ax.xaxis.set_major_locator( mdates.MinuteLocator(interval=dmin) )
   ax.tick_params( axis='x', labelsize=8 )

   #ax.set_xlim( itimes_l[0], itimes_l[-1] )
   ax.set_xlim( stime_, etime_ )

   ax.text( 0.5, 1.01, tit_l[i],
          fontsize=14, transform=ax.transAxes,
          ha="center", va="bottom", )
          #bbox=bbox )
 
   ax.set_ylabel( tit_l[i], fontsize=12 )
    

ax2.hlines( y=1.0, xmin=stime, xmax=etime, ls='dashed',
            linewidth=0.5, color='gray' )  

ax2.text( 0.01, 0.05, "Threshold: {0:.0f}dBZ\nZ={1:.0f}km".format( dbz, theight/1000.0, ),
          fontsize=11, transform=ax2.transAxes, 
          ha="left", va="bottom", 
          bbox=bbox)

ax1.set_ylim( ymin1, ymax1 )
ax2.set_ylim( ymin2, ymax2 )

ax1.set_xlabel( "Valid time (JST)", fontsize=10 )
ax2.set_xlabel( "Valid time (JST)", fontsize=10 )


pnum_l = [ "(a)", "(b)" ]
ax_l = [ ax1, ax2 ] 
for i, ax in enumerate( ax_l ):
    ax.text( 0.01, 0.95, pnum_l[i],
              fontsize=11, transform=ax.transAxes,
              ha='left',
              va='top', 
              bbox=bbox )
    

ofig = "2p_realtime_score_thrs{0:.1f}dbz_z{1:.1f}.png".format( dbz, theight )

print( ofig )
if quick:
   plt.show()
else:
   plt.savefig( ofig,
                bbox_inches="tight", pad_inches = 0.1)
   plt.clf()
   plt.close('all')



