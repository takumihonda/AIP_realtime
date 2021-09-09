import os
import sys
import numpy as np 
from datetime import datetime, timedelta

quick = True
#quick = False

#fn_h = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/realtime_log_dafcst_nc20200807.txt"
fn_h = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/honda.txt"
fn_a = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/amemiya.txt"

# JST
stime = datetime( 2020, 7,  31, 11, 20 )
etime = datetime( 2020, 8,  7,  14, 0 )

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
#        print(line, end="")
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
print( len(lt_l) )
lt_l = read_files( fn=fn_h, ftime_l=ftime_l, stime=stime, etime=etime, AMEMIYA=True )
#print( len(lt_l) )
#sys.exit()

ftime_l = np.array( ftime_l )
lt_l = np.array( lt_l )


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates

#plt.plot( time_l )
#plt.plot( ftime_l )

fig, ((ax)) = plt.subplots( 1, 1, figsize=( 10, 6.4 ) )
fig.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.98, )

# lead time

lt_l = lt_l / 60.0 # minute

# When forecast becomes available after obs
ft_l = 30.0 - lt_l
#ft_l[ np.isnan(ft_l) ] = 31
ft_l[ ft_l > 30.0 ] = np.nan


total = len( ft_l[ ~np.isnan(ft_l) ] )
total_nan = len( ft_l[ np.isnan(ft_l) ] )


xmin = 0.0
xmax = 30.0
dx = 2.0

#xlabs = np.arange( xmin, 35.0, 5.0 )
xlabs = np.arange( xmin, xmax+dx, dx )

#xlabs_c = [ '{0:}'.format( x/2 ) for x in range(len(xlabs)) ]
#print(xlabs)
#sys.exit()



bins = int( ( xmax - xmin ) / 0.5 ) 

#print( np.nanmax(lt_l) )
#print( np.nanmin(lt_l) )
y, x, _ = ax.hist( ft_l[ ~np.isnan(ft_l) ], range=(xmin,xmax), 
             bins=bins, rwidth=0.8 )

ymax = y.max()
ymax = 8000

ylab = "Count"
ax.set_ylabel( ylab, fontsize=13 )

#xlab = "Time for 30-min forecast (min)"
xlab = "Process time for 30-min forecast (min)"
ax.set_xlabel( xlab, fontsize=13 )

ax.set_xlim( xmin, xmax )
ax.set_xticks( xlabs )
ax.set_yticks( np.arange( 0, 8000, 1000 ) )
#ax.set_xticklabels( xlabs_c )
ax.set_ylim( y.min(), ymax )

xobs = 10.0 / 60.0 
ax.vlines( x=xobs, ymin=0, ymax=ymax, color='r', ls='solid', lw=1.0)

arrow_dict = dict(arrowstyle = "->", color = "r", )
                  #connectionstyle = "angle, angleA=315, angleB=315")

ax.annotate("Obs\nreceived(10s)", xy=(xobs, 0), size=10, xytext=(xobs+0.5, -500),
            color="r", arrowprops=arrow_dict )

#ax.text( 31, -100, "No\nfcst",
#          fontsize=12, #transform=ax.transAxes,
#          ha='center', va='top' )

xlabs_ = np.arange( xmin, xmax, 1 )
ax.vlines( x=xlabs_, ymin=0, ymax=ymax, color='gray', ls='dotted', lw=0.3)


#tit = "When 30-min forecasts finished"
#ax.text( 0.5, 1.01, tit,
#          fontsize=15, transform=ax.transAxes,
#          horizontalalignment='center',
#          verticalalignment='bottom' )

rdtime = timedelta( seconds=30*total )
print( rdtime)
period = "Total: {0:} (active for {1:})\nFrom {2:}\nto {3:}".format( total, rdtime,
                                stime.strftime('%H:%M:%SJST %m/%d'),
                                etime.strftime('%H:%M:%SJST %m/%d'))
ax.text( 0.9, 0.9, period,
          fontsize=12, transform=ax.transAxes,
          ha='right', va="top" )


ofig = "hist_realtime0807_leadtime.png"

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




sys.exit()

