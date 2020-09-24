import numpy as np
from datetime import datetime, timedelta 
import matplotlib.dates as mdates

import matplotlib.pyplot as plt

quick = False
#quick = True

def prep_dat( stime=datetime( 2020, 8, 27, 0, 0 ), etime=datetime( 2020, 8, 28, 0, 0 ), 
              dtmin=1 ):
    
    times = []
    time = stime
    while time <= etime:
        times.append( time )
        time += timedelta( minutes=dtmin )
    
    dat1d = np.zeros( ( len(times) ) )
    
    return( times, dat1d )    


stime = datetime( 2020, 8, 27, 3, 0 )
etime = datetime( 2020, 8, 28, 3, 0 )

ymin = 0
ymax = 1200

ylevs = np.arange( 0, 1400, 200 )

times = []
time = stime
while time <= etime:
    times.append( time )
    time += timedelta( minutes=1 )

dat1d = np.zeros( ( len(times) ) )

# every 3 h
times3h = []
time = stime
while time <= etime:
    times3h.append( time )
    time += timedelta( minutes=180 )

lw = 0.5
ec = 'k'

fig, (ax1) = plt.subplots( 1, 1, figsize=( 10,5 ) )
fig.subplots_adjust(left=0.08, bottom=0.07, right=0.97, top=0.94, )

ax1.plot( times, dat1d )

ax1.set_ylim( ymin, ymax )
ax1.set_xlim( stime, etime )
ax1.set_yticks( ylevs )

ax1.xaxis.set_major_formatter( mdates.DateFormatter('%HJST') )
ax1.xaxis.set_major_locator( mdates.MinuteLocator( interval=180 ) )

ax1.hlines( xmin=stime, xmax=etime, y=np.arange( 1000, 2000, 1000), 
       ls='dashed', lw=0.5 )

ax1.vlines( ymin=ymin, ymax=ymax, x=times3h,
       ls='dashed', lw=0.5 )

ylab = "Node"
ax1.set_ylabel( ylab, fontsize=13 )

alp = 0.6

c4 = 'r'
c1 = 'gray'
c12 = 'b'
c3 = 'g'
c4i = 'orange'

# D4
d4_node = 992
labd4 = "D4 DA cycle & fcst\n({0:} nodes)".format( d4_node )
stime_ = stime
etime_ = stime + timedelta( hours=4 ) - timedelta( minutes=5 )
times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
ax1.fill_between( times, dat1d, dat1d+d4_node, facecolor=c4, alpha=alp,
                  label=labd4,
                  edgecolor=ec, lw=lw )

for i in range( 2 ):
    stime_ = etime_ + timedelta( minutes=5 )
    etime_ = stime_ + timedelta( hours=12 ) - timedelta( minutes=5 )
    times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
    ax1.fill_between( times, dat1d, dat1d+d4_node, facecolor=c4, alpha=alp,
                       label=None,
                       edgecolor=ec, lw=lw )

# D1 cycle
min_d1c = 20
d1c_node = 156

stime1 = stime + timedelta( minutes=40 )

# D1/D2 fcst
min_d12 = 45
d12_node = 182
d12it = 2

# D3 fcst
min_d3 = 45
d3_node = 208
d3it = 1

# D4 init
min_d4i = 10
d4i_node = 208
d4it = 4

for i in range( 4 ):

    if i == 0:
       labd4i = "D4 init\n({0:} nodes, {1:}min x {2:})".format( d4i_node, min_d4i, d4it)
       labd3 = "D3 fcst\n({0:} nodes, {1:}min)".format( d3_node, min_d3 )
       labd12 = "D1&D2 fcst\n({0:} nodes, {1:}min x {2})".format( d12_node, min_d12, d12it )
       labd1 = "D1 DA cycle\n({0:} nodes, {1:}min)".format( d1c_node, min_d1c )
    else:
       labd4i = None
       labd3 = None
       labd12 = None
       labd1 = None

    stime_ = stime1 + timedelta( hours=6*i )
    etime_ = stime_ + timedelta( minutes=min_d1c )
    
    times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
    ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d1c_node, 
                       facecolor=c1, alpha=alp, label=labd1,
                       edgecolor=ec, lw=lw )
    
    # D1/D2 fcst
    for j in range( d12it ):
        if j > 0:
           labd12 = None
        stime_ = etime_
        etime_ = stime_ + timedelta( minutes=min_d12 )
        times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
        ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d12_node, 
                           facecolor=c12, alpha=alp, label=labd12,
                           edgecolor=ec, lw=lw )
    
    
    # D3 fcst
    stime_ = etime_
    etime_ = stime_ + timedelta( minutes=min_d3 )
    
    times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
    ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d3_node, 
                       facecolor=c3, alpha=alp, label=labd3,
                       edgecolor=ec, lw=lw )
    
    
    # D4 init
    for j in range( d4it ):
        if j > 0:
           labd4i = None
        stime_ = etime_
        etime_ = stime_ + timedelta( minutes=min_d4i )
        times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
        ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d4i_node, 
                           facecolor=c4i, alpha=alp, label=labd4i,
                           edgecolor=ec, lw=lw )
    
ax1.legend( loc="lower left", fontsize=11, framealpha=0.9 )

tit = "Job schedule"
ax1.text( 0.5, 1.01, tit,
          fontsize=15, transform=ax1.transAxes,
          ha='center',
          va='bottom' )

ofig = "1p_job_schedule.png"
    
print( ofig )
if quick:
   plt.show()
else:
   plt.savefig( ofig,
                bbox_inches="tight", pad_inches = 0.1)
   plt.clf()
   plt.close('all')

