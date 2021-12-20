import numpy as np
from datetime import datetime, timedelta 
import matplotlib.dates as mdates

import matplotlib.pyplot as plt
import matplotlib

quick = False
#quick = True

nmpi = 64

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
ymax3 = 1200
ymax = 1200 * nmpi

#ylevs = np.arange( 0, 1400*nmpi, 200*nmpi )
ylevs = np.arange( 0, 80000, 10000 )
ylevs = np.append( ylevs, 1200*nmpi)
ylevs3 = np.arange( 0, 1400, 200 )

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
fig.subplots_adjust(left=0.09, bottom=0.1, right=0.91, top=0.94, )

ax1.plot( times, dat1d )

dat1d_ = np.copy( dat1d )
dat1d_[:] = np.nan

ax3 = ax1.twinx()
ax3.plot( times, dat1d_ )
ax3.set_yticks( ylevs3 )
ax3.set_ylim( ymin, ymax3 )

ax3.get_yaxis().set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda y, p: format(int(y), ',')))

ax1.set_ylim( ymin, ymax )
ax1.set_xlim( stime, etime )
ax1.set_yticks( ylevs )

ax1.xaxis.set_major_formatter( mdates.DateFormatter('%HJST') )
ax1.xaxis.set_major_locator( mdates.MinuteLocator( interval=180 ) )

ax1.get_yaxis().set_major_formatter(
    matplotlib.ticker.FuncFormatter(lambda y, p: format(int(y), ',')))


ax2 = ax1.twiny()
#ax2.spines["right"].set_position(("axes", 1.2))
#ax2.spines["right"].set_visible(True)
ax2.plot( times, dat1d_ )
ustime = stime - timedelta(hours=9)
uetime = etime - timedelta(hours=9)
ax2.set_xlim( ustime, uetime )
ax2.xaxis.set_major_formatter( mdates.DateFormatter('(%HUTC)') )
ax2.xaxis.set_major_locator( mdates.MinuteLocator( interval=180 ) )
ax2.spines["top"].set_position(("axes", -0.12))
ax2.spines["top"].set_visible( False )


ax1.tick_params( axis='x', labelsize=10 )
ax2.tick_params( axis='x', labelsize=9 )

ax1.hlines( xmin=stime, xmax=etime, y=np.arange( 1000, 2000, 1000), 
       ls='dotted', lw=0.4 )

ax1.vlines( ymin=ymin, ymax=ymax, x=times3h,
       ls='dotted', lw=0.4 )

ylab = "Node"
ylab = "Number of MPI processes"
ax1.set_ylabel( ylab, fontsize=13 )

ylab3 = "Number of Computation nodes on OFP"
ax3.set_ylabel( ylab3, fontsize=13 )

alp = 0.6

c4 = 'r'
c1 = 'gray'
c12 = 'b'
c3 = 'g'
c4i = 'orange'

# D4
d4_node = 992
d4_prc = d4_node*nmpi
labd4 = "D4 DA cycle & fcst\n({0:,} MPI prc)".format( d4_prc )
#stime_ = stime
stime_ = datetime( 2020, 8, 27, 4, 0 ) - timedelta( hours=12 )
etime_ = stime_ - timedelta( minutes=5 )
times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
#ax1.fill_between( times, dat1d, dat1d+d4_node, facecolor=c4, alpha=alp,
ax1.fill_between( times, dat1d, dat1d+d4_prc, facecolor=c4, alpha=alp,
                  label=labd4,
                  edgecolor=ec, lw=lw )

for i in range( 3 ):
    stime_ = etime_ + timedelta( minutes=5 )
    etime_ = stime_ + timedelta( hours=12 ) - timedelta( minutes=5 )
    times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
    #ax1.fill_between( times, dat1d, dat1d+d4_node, facecolor=c4, alpha=alp,
    ax1.fill_between( times, dat1d, dat1d+d4_prc, facecolor=c4, alpha=alp,
                       label=None,
                       edgecolor=ec, lw=lw )

# D1 cycle
min_d1c = 20
d1c_node = 156
d1c_prc = d1c_node * nmpi

stime1 = stime + timedelta( minutes=40 )

# D1/D2 fcst
min_d12 = 30
d12_node = 182
d12_prc = d12_node * nmpi
d12it = 2

# D3 fcst
#min_d3 = 45
min_d3 = 105
d3_node = 208
d3_prc = d3_node * nmpi
d3it = 1

# D4 init
min_d4i = 15
d4i_node = 208
d4i_prc = d4i_node * nmpi
d4it = 4

for i in range( 4 ):

    if i == 0:
       labd4i = "D4 boundary\n({0:,} MPI prc, {1:} min x {2:})".format( d4i_prc, min_d4i, d4it)
       labd3 = "D3 fcst\n({0:,} MPI prc, {1:} min)".format( d3_prc, min_d3 )
       labd12 = "D1&D2 fcst\n({0:,} MPI prc, {1:} min x {2})".format( d12_prc, min_d12, d12it )
       labd1 = "D1 DA cycle\n({0:,} MPI prc, {1:} min)".format( d1c_prc, min_d1c )
    else:
       labd4i = None
       labd3 = None
       labd12 = None
       labd1 = None

    stime_ = stime1 + timedelta( hours=6*i )
    etime_ = stime_ + timedelta( minutes=min_d1c )
    
    times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
    #ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d1c_node, 
    ax1.fill_between( times, dat1d+d4_prc, dat1d+d4_prc+d1c_prc, 
                       facecolor=c1, alpha=alp, label=labd1,
                       edgecolor=ec, lw=lw )
    
    # D1/D2 fcst
    for j in range( d12it ):
        if j > 0:
           labd12 = None
        stime_ = etime_
        etime_ = stime_ + timedelta( minutes=min_d12 )
        times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
        #ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d12_node, 
        ax1.fill_between( times, dat1d+d4_prc, dat1d+d4_prc+d12_prc, 
                           facecolor=c12, alpha=alp, label=labd12,
                           edgecolor=ec, lw=lw )
    
    
    # D3 fcst
    stime_ = etime_
    etime_ = stime_ + timedelta( minutes=min_d3 )
    
    times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
    #ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d3_node, 
    ax1.fill_between( times, dat1d+d4_prc, dat1d+d4_prc+d3_prc, 
                       facecolor=c3, alpha=alp, label=labd3,
                       edgecolor=ec, lw=lw )
    
    
    # D4 init
    for j in range( d4it ):
        if j > 0:
           labd4i = None
        stime_ = etime_
        etime_ = stime_ + timedelta( minutes=min_d4i )
        times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
        #ax1.fill_between( times, dat1d+d4_node, dat1d+d4_node+d4i_node, 
        ax1.fill_between( times, dat1d+d4_prc, dat1d+d4_prc+d4i_prc, 
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

