import numpy as np
import sys

import matplotlib.pyplot as plt
from datetime import datetime, timedelta

quick = False
quick = True

# color
c1 = 'firebrick'
c2 = c1

c_mp = 'goldenrod'

xsize = 5
ysize = 5

fig, (ax) = plt.subplots( 1, 1, figsize=( xsize,ysize ) )
fig.subplots_adjust(left=0.08, bottom=0.1, right=0.97, top=0.94, )

ax.axis( 'off' )

#plt.show()

xmin_ = -24
xmax_ = 48

xmin = 0
xmax = 19

# dx=1 corresponds 10 s
dx30m = 1*6*30
dy30m = 120 

# D1 DA cycle

y1 = 0
dx = 3

boxdic1da = {
    "facecolor" : c1,
    "edgecolor" : c1,
    "boxstyle" : "Round",
    "linewidth" : 2
}

xmp = 12
ymp = y1+4

time0 = datetime( 2000, 1, 1, 0, 0, 0 )
for x in range( xmin_, xmax_, 3 ):
   time = time0 + timedelta( seconds=10*x )
   ctime = time.strftime( '%H:%M:%S' )    

   ax.annotate(s='',xy=(x+dx,y1),xytext=(x,y1),xycoords='data',
               annotation_clip=False,
               arrowprops=dict(facecolor=c1, width=20.0,
                               edgecolor=c1,
                               headwidth=40.0,headlength=7.0,
                               shrink=0.01)
              )

   if x < 18:
   
       ax.text( x-1, y1+1, '{0:}'.format( ctime ),
                 fontsize=10, #transform=ax1.transAxes,
                 ha='left',
                 va='bottom' )

   # MP-PAWR
   if x >= -3 and x <= 18:
      ax.annotate(s='',xytext=(xmp,ymp),xy=(x,y1+1.5),
                  xycoords='data',
                  annotation_clip=False,
                  arrowprops=dict(facecolor=c_mp, width=1.0,
                                  edgecolor=c_mp,
                                  headwidth=4.0, headlength=2.0,
                                  shrink=0.01)
                 )

   if x == 12:
      ax.text( x+3.0, y1-1.3, 'D4 DA cycling',
               zorder=10,
               fontsize=12,
               color='w',
               bbox=boxdic1da,
               ha='center',
               va='center' )


ax.set_xlim( xmin, xmax )
#ax.set_xlim( -1, 60 )
ax.set_ylim( -10, 4 )


# D4 fcst

y2 = y1 - 0.4
dy = -dy30m
dx = dx30m

for x in range( xmin_, xmax_, 3 ):
   ax.annotate(s='',xy=(x+dx,y2+dy),xytext=(x,y2),xycoords='data',
               annotation_clip=False,
               arrowprops=dict(facecolor=c2, width=1.0,
                               edgecolor=c2,
                               headwidth=10.0,headlength=7.0,
                               ) )

   for j in range( 1, 3, 1 ):
       ax.annotate(s='',xy=(x+dx/30*j,y2+dy/30*j),xytext=(x,y2),xycoords='data',
                   annotation_clip=False,
                   arrowprops=dict(facecolor=c2, width=1.0,
                                   edgecolor=c2,
                                   headwidth=10.0,headlength=7.0,
                                   ) )

   if x == 0:
    
      x_ = x + dx30m/20 + 0.5
      y_ = y2 - dy30m/20 + 0.5
      dx_ = dx/30
      dy_ = dy/30 
      angle = np.rad2deg( np.arctan( dy_/dx_ * xsize/ysize  ) ) - 6
      ax.text( x_, y_, 'D4 extended forecast (initiated every 30 s)',
               color=c2,
               rotation=angle,
               fontsize=12,
               ha='center',
               va='center' )

plt.show()

sys.exit()




#ofig = "pdf/D1-D3_schematic.pdf"
ofig = "png/D1-D3_schematic.png"

print( ofig )
if quick:
   plt.show()
else:
   plt.savefig( ofig,
                bbox_inches="tight", )#pad_inches = 0.1)
   plt.clf()
   plt.close('all')



sys.exit()


ax1.legend( loc="lower left", fontsize=11, framealpha=0.9 )

tit = "Job schedule"
ax1.text( 0.5, 1.01, tit,
          fontsize=15, transform=ax1.transAxes,
          ha='center',
          va='bottom' )

    
