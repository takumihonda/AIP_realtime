import os
import sys
import numpy as np 
from datetime import datetime, timedelta

quick = True
quick = False

#fn_h = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/realtime_log_dafcst_nc20200807.txt"
fn_h = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/honda.txt"
fn_a = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/amemiya.txt"

#stime = datetime( 2020, 7,  31, 18, 0 )
#etime = datetime( 2020, 8,  7, 18, 0 )

stime = datetime( 2020, 7,  31, 11, 20 )
etime = datetime( 2020, 8,  7,  14, 0 )

stime_ = datetime( 2020, 7,  31, 12, 0 )
etime_ = datetime( 2020, 8,  7,  18, 0 )

dt = timedelta( seconds=30 )

FT_sec = 1800
ftime_delta = timedelta( seconds=FT_sec )

ftime_l = []
lt_l = []


time = stime
while time <= etime:

   ftime_l.append( time ) # initial time
   lt_l.append( np.nan )

   time += dt




def read_files( fn="", ftime_l=[], 
           stime=datetime( 2020, 7, 31, 12, ), 
           etime=datetime(2020, 8, 1, 12, 0), AMEMIYA=True ):

    j = 0
    
    f = open( fn )
    lines = f.readlines()
    for line in lines:
        print(line, end="")
        data = line.split(' ')
    
        cdate = data[0]
        ctime = data[1]
        cname = str( data[2] ).strip()
    
    
        # When 30-min forecast is available
        time_ = datetime( year=int(cdate[0:4]), 
                          month=int(cdate[5:7]),
                          day=int(cdate[9:11]),
                          hour=int(ctime[0:2]),
                          minute=int(ctime[3:5]),
                          second=int(ctime[6:8]),
                        )
       
        if j == 0:
           i = 9
        else:
           i = 5
    
        if AMEMIYA:
           i = 0
    
#        print( cname )
#        print( cname[i:i+4] )
    #    print( cname[i+4:i+6] )
    #    print( cname[i+6:i+8] )
    #    print( cname[i+9:i+11] )
    #    print( cname[i+11:i+13] )
    #    print( cname[i+13:i+15] )
    
    #    if j == 10:
    #       sys.exit()
    
        # 30-min forecast initial time
        ftime_ = datetime( year=int(cname[i:i+4]), 
                           month=int(cname[i+4:i+6]),
                           day=int(cname[i+6:i+8]), 
                           hour=int(cname[i+9:i+11]),
                           minute=int(cname[i+11:i+13]),
                           second=int(cname[i+13:i+15]),
                         )
    
    
        j += 1
    
        dt_ = ( ftime_ - stime ).total_seconds() 
        idx_ = int( dt_ / 30 )
        #print( dt_, ftime_, stime, idx_, ftime_l[idx_]  )
    
        lt_l[idx_] = ( ftime_ + ftime_delta - time_ ).total_seconds()
    
    return( lt_l )

##


lt_l = read_files( fn=fn_a, ftime_l=ftime_l, stime=stime, etime=etime, AMEMIYA=True )
lt_l = read_files( fn=fn_h, ftime_l=ftime_l, stime=stime, etime=etime, AMEMIYA=True )

ftime_l = np.array( ftime_l )
lt_l = np.array( lt_l )


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates

#plt.plot( time_l )
#plt.plot( ftime_l )

fig, ((ax)) = plt.subplots( 1, 1, figsize=( 12, 6.5 ) )
fig.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.95, )

# lead time

lt_l = lt_l / 60.0 # minute
ax.plot( ftime_l, lt_l, color='k', lw=1.0 )

#tlev_rm = 120*3
#kernel = np.ones(tlev_rm) / tlev_rm
#lt_rm = np.convolve( lt, kernel, mode='same') 
#ax.plot( ftime_l, lt_rm, color='r', lw=1.0 )

#ax.xaxis.set_major_locator(mdates.HourLocator(interval=1) )
ax.xaxis.set_major_locator(mdates.HourLocator(interval=12) )
ax.xaxis.set_major_formatter(mdates.DateFormatter('%HJST\n%m/%d'))

ymin_ = 0
ymax_ = 30
dy_ = 5
ax.set_ylim( ymin_, ymax_ )

ylab = "Forecast lead time (minute)"
ax.set_ylabel( ylab, fontsize=13 )

xlab = "Forecast initial time"
ax.set_xlabel( xlab, fontsize=13 )

ylevs = np.arange( ymin_, ymax_+dy_, dy_ )
ax.set_yticks( ylevs )

xlabs = []
time_ = stime
while time_ <= etime:
   xlabs.append( time_ )
   time_ += timedelta( hours=12 )

ax.hlines( y=ylevs, xmin=stime_, xmax=etime_, color='gray', ls='dashed', lw=0.5)
ax.vlines( x=xlabs, ymin=ymin_, ymax=ymax_, color='gray', ls='dashed', lw=0.5)

tit = "Lead time of 30-min forecasts"
ax.text( 0.5, 1.01, tit,
          fontsize=15, transform=ax.transAxes,
          horizontalalignment='center',
          verticalalignment='bottom' )

ax.set_xlim( stime_, etime_ )
print( ftime_l[0], ftime_l[-1])

ofig = "realtime0807_leadtime.png"

print( ofig )
if quick:
   plt.show()
else:
   odir = "png/realtime0807"
   os.makedirs( odir, exist_ok=True)
   plt.savefig( os.path.join(odir, ofig),
                bbox_inches="tight", pad_inches = 0.1)
   plt.clf()
   plt.close('all')


