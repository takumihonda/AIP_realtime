import numpy as np
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#quick = False
quick = True

# color
c1 = 'k'
c3 = 'b'
c4 = 'r'

fig, (ax) = plt.subplots( 1, 1, figsize=( 10,5 ) )
fig.subplots_adjust(left=0.08, bottom=0.1, right=0.97, top=0.94, )

#plt.show()

xmin_ = -24
xmax_ = 54

xmin = 0
xmax = 24

dx24 = 24
dy24 = 8 

# D1 DA cycle

y1 = 0
dx = 6

boxdic1da = {
    "facecolor" : "black",
    "edgecolor" : "black",
    "boxstyle" : "Round",
    "linewidth" : 2
}

for x in range( xmin_, xmax_, 6 ):

   ax.annotate(s='',xy=(x+dx,y1),xytext=(x,y1),xycoords='data',
               annotation_clip=False,
               arrowprops=dict(facecolor=c1, width=20.0,
                               headwidth=40.0,headlength=7.0,
                               shrink=0.01)
              )

   if x <= 18:
       ax.text( x, y1+1, '{0:0=4} UTC'.format( x*100 ),
                 fontsize=10, #transform=ax1.transAxes,
                 ha='center',
                 va='bottom' )

   if x == 18:
      ax.text( x+3, y1+1.2, 'D1 DA cycling',
               fontsize=12,
               color='w',
               bbox=boxdic1da,
               ha='center',
               va='center' )



# D1 & D2 fcst

y2 = y1 - 0.4
dy = -dy24
dx = dx24

for x in range( xmin_, xmax_, 6 ):
   ax.annotate(s='',xy=(x+dx,y2+dy),xytext=(x,y2),xycoords='data',
               annotation_clip=False,
               arrowprops=dict(facecolor='k', width=5.0,
                               headwidth=10.0,headlength=7.0,
                               ) )

   if x == 0:
    
      x_ = x + dx24/24*6 + 2
      y_ = y2 - dy24/24*6 #-4
      angle = np.rad2deg( np.arctan( dy/dx ) )
      ax.text( x_, y_, 'D1 & D2 1-day extended forecast',
               rotation=angle,
               fontsize=12,
               ha='center',
               va='center' )


# D3 fcst

y3 = y2 - 0.5
dx = dx24 / 24 * 7
dy = -dy24 / 24 * 7

xs = dx24 / 24 * 17
ys = -dy24 / 24 * 17

y3 += ys

for x in range( xmin_, xmax_, 6 ):
   x_ = x + xs
   ax.annotate(s='',xy=(x_+dx,y3+dy),xytext=(x_,y3),xycoords='data',
               annotation_clip=False,
               arrowprops=dict(facecolor=c3, width=5.0,
                               headwidth=10.0,headlength=7.0,
                               shrink=0.00)
              )

   if x == 6:
    
      x_ = x_ - dx24/24*9
      y_ = y3 + dy24/24*3
      angle = np.rad2deg( np.arctan( dy/dx ) )
      ax.text( x_, y_, 'D3 7-h extended forecast',
               rotation=angle,
               fontsize=12,
               color=c3,
               ha='center',
               va='center' )

# D4 data
y3 -= 1.0
dx = dx24 / 24 * 6
dy = -dy24 / 24 * 6

xs = dx24 / 24 * 18
ys = -dy24 / 24 * 18

for x in range( xmin_, xmax_, 6 ):
   x_ = x + xs
   ax.annotate(s='',xy=(x_+dx,y3+dy),xytext=(x_,y3),xycoords='data',
               annotation_clip=False,
               arrowprops=dict( color=c4, linewidth=2.0, arrowstyle="<->",
                                mutation_scale=25, 
                                #linestyle="--",
                              ),
               )
                               #headwidth=10.0, headlength=7.0,
                               #shrink=0.01, ))#arrowstyle="<->"),)
   if x != -1:
    
      x_ = x_ + 2 #- dx24/24*9
      y_ = y3 - 1.5 #+ dy24/24*3
      sh = ( x+18 )*100
      sh2 = ( x+24 )*100
      ax.text( x_, y_, 'D4 boundary data\n({0:0=4}-{1:0=4})'.format( sh, sh2 ),
               rotation=angle,
               fontsize=12,
               color=c4,
               ha='center',
               va='center' )



ax.set_xlim( xmin, xmax )
#ax.set_xlim( -1, 60 )
ax.set_ylim( -10, 2 )
plt.show()





ax.set_xlim( xmin, xmax )
#ax.set_xlim( -1, 60 )
ax.set_ylim( -10, 1 )
plt.show()

sys.exit()


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

ax1.plot( times, dat1d )

ax1.set_ylim( ymin, ymax )
ax1.set_xlim( stime, etime )
ax1.set_yticks( ylevs )

ax1.xaxis.set_major_formatter( mdates.DateFormatter('%HJST') )
ax1.xaxis.set_major_locator( mdates.MinuteLocator( interval=180 ) )


dat1d_ = np.copy( dat1d )
dat1d_[:] = np.nan

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
#stime_ = stime
stime_ = datetime( 2020, 8, 27, 4, 0 ) - timedelta( hours=12 )
etime_ = stime_ - timedelta( minutes=5 )
times, dat1d = prep_dat( stime=stime_, etime=etime_, dtmin=1 )
ax1.fill_between( times, dat1d, dat1d+d4_node, facecolor=c4, alpha=alp,
                  label=labd4,
                  edgecolor=ec, lw=lw )

for i in range( 3 ):
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
       labd4i = "D4 boundary\n({0:} nodes, {1:}min x {2:})".format( d4i_node, min_d4i, d4it)
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

