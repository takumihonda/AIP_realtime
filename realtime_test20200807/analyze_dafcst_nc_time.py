import os
import sys
import numpy as np 
from datetime import datetime, timedelta

quick = True
#quick = False

#fn_h = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/realtime_log_dafcst_nc20200807.txt"
fn_h = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/honda.txt"
fn_a = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/amemiya.txt"

stime = datetime( 2020, 7,  31, 12, 0 )
etime = datetime( 2020, 8,  7, 14, 0 )

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

def get_JMA_rmax( rmin=30.0 ):


    stime_ = datetime( 2020, 7, 31, 0 )
    etime_ = datetime( 2020, 8,  7, 0 )


    of = "../python/JMA_rmax_rmin{0:.0f}_{1:}_{2:}.npz".format( rmin, 
                                            stime_.strftime( '%Y%m%d' ),    
                                            etime_.strftime( '%Y%m%d' ),    
                                          )

    data = np.load( of, allow_pickle=True )
    
    rmax_l = data["rmax_l"]
    rarea_l = data["rarea_l"]
    time_l = data["time_l"]
    
    # UTC2JST
    time_l += timedelta( hours=9 )

    return( rmax_l, rarea_l, time_l )


def read_files( fn="", ftime_l=[], 
           stime=datetime( 2020, 7, 31, 12, ), 
           etime=datetime(2020, 8, 1, 12, 0), AMEMIYA=True ):

    j = 0
    
    f = open( fn )
    lines = f.readlines()
    for line in lines:
#        print(line, end="")
        data = line.split(' ')
    
        cdate = data[0]
        ctime = data[1]
        cname = str( data[2] ).strip()
    
    
        # When 30-min forecast is available
        time_ = datetime( year=int(cdate[0:4]), 
                          month=int(cdate[5:7]),
                          day=int(cdate[8:10]),
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
#        if ftime_ < datetime( 2020, 8, 1, 0 ):
#           print( dt_, ftime_, time_, idx_, ftime_l[idx_]  )
    
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
fig.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, )

# lead time

lt_l = lt_l / 60.0 # minute
#lt_l[lt_l<0] = np.nan
#lt_l[lt_l>30] = np.nan
print( np.nanmax(lt_l), np.nanmin(lt_l) )
#sys.exit()
ax.plot( ftime_l, lt_l, color='k', lw=1.0 )

#tlev_rm = 120*3
#kernel = np.ones(tlev_rm) / tlev_rm
#lt_rm = np.convolve( lt, kernel, mode='same') 
#ax.plot( ftime_l, lt_rm, color='r', lw=1.0 )

#ax.xaxis.set_major_locator(mdates.HourLocator(interval=1) )

#ax.xaxis.set_major_locator(mdates.HourLocator(interval=6) )
#ax.xaxis.set_major_formatter(mdates.DateFormatter('%HJST\n%m/%d'))
ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6), tz=None) )
ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))

ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(12, 13, 1), tz=None) )
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H\n%m/%d'))

ymin_ = 0
ymax_ = 30
dy_ = 5
ax.set_ylim( ymin_, ymax_ )

ax.vlines( ftime_l[ np.isnan( lt_l ) ], ymin=ymin_, ymax=ymax_, color='gray', lw=0.01, alpha=0.5 )

ylab = "Forecast lead time (minute)"
ax.set_ylabel( ylab, fontsize=13 )

xlab = "Forecast initial time (JST)"
ax.set_xlabel( xlab, fontsize=13 )

ylevs = np.arange( ymin_, ymax_+dy_, dy_ )
ax.set_yticks( ylevs )

xlabs = []
time_ = stime
while time_ <= etime:
   xlabs.append( time_ )
   time_ += timedelta( hours=12 )

ax.hlines( y=ylevs, xmin=stime, xmax=etime, color='gray', ls='dashed', lw=0.5)
ax.vlines( x=xlabs, ymin=ymin_, ymax=ymax_, color='gray', ls='dashed', lw=0.5)

tit = "Lead time of 30-min forecasts"
ax.text( 0.5, 1.01, tit,
          fontsize=15, transform=ax.transAxes,
          horizontalalignment='center',
          verticalalignment='bottom' )

ax.set_xlim( stime, etime )
print( ftime_l[0], ftime_l[-1])
print( stime, etime )



ax2 = ax.twinx()

ymin2_ = 0.0
ymax2_ = 80.0

ylevs2 = [ 0, 8, 16, 24, 32, 40 ]

rmin = 30.0
rmin_l = [ 1.0, 20.0, ]
#rmin_l = [ 5.0, ]
cc_l = [ "cyan", "b", ]
for k, rmin in enumerate( rmin_l ):

    rmax_l, rarea_l, jtime_l = get_JMA_rmax( rmin=rmin )
    print( k, len( rmax_l ) )

    ax2.plot( jtime_l, rarea_l/100.0, color=cc_l[k], lw=0.5, 
         label='>={0:.0f}mm h$^{{-1}}$'.format( rmin ) )
    ax2.set_ylim( ymin2_, ymax2_ )
    ax2.set_yticks( ylevs2 )
    ax2.tick_params( axis='y', labelsize=8 ) 
    ax2.yaxis.label.set_color( 'b' )
#    ax2.set_xlim( stime_, etime_ )

ax2.legend( bbox_to_anchor=( 0.05, 0.1), 
             loc='lower left', fontsize=9, ).get_frame().set_alpha( 1.0 )

ylab2 = 'Precipitation area from JMA radar (x10$^2$km$^2$)\n(where >{0:.0f})mm h$^{{-1}}$)'.format( rmin )
ax2.set_ylabel( ylab2, fontsize=9 )

ax2.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(0, 24, 6), tz=None) )
ax2.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))

ax2.xaxis.set_major_locator(mdates.HourLocator(byhour=range(12, 13, 1), tz=None) )
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H\n%m/%d'))

ax2.tick_params(axis='y', colors='b')
ax2.yaxis.set_label_coords( 1.03, 0.15 )


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


