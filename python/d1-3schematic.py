import numpy as np
import sys

import matplotlib.pyplot as plt

quick = False
quick = True

# color
c1 = 'gray'
c2 = 'k'
c3 = 'royalblue'
c4 = 'firebrick'

xsize = 7
ysize = 5
fig, (ax) = plt.subplots( 1, 1, figsize=( xsize,ysize ) )
fig.subplots_adjust(left=0.08, bottom=0.1, right=0.97, top=0.94, )

ax.axis( 'off' )

#plt.show()

xmin_ = -24
xmax_ = 24

xmin = 0
xmax = 24

dx24 = 24
dy24 = 8 

# D1 DA cycle

y1 = 0
dx = 6

boxdic1da = {
    "facecolor" : c1,
    "edgecolor" : c1,
    "boxstyle" : "Round",
    "linewidth" : 2
}

for x in range( xmin_, xmax_, 6 ):

   ax.annotate(s='',xy=(x+dx,y1),xytext=(x,y1),xycoords='data',
               annotation_clip=False,
               arrowprops=dict(facecolor=c1, width=20.0,
                               edgecolor=c1,
                               headwidth=40.0,headlength=7.0,
                               shrink=0.01)
              )

   if x <= 18:
       ax.text( x-1, y1+1, 'NCEP PREPBUFR/GFS\n{0:0=4} UTC'.format( x*100 ),
                 fontsize=10, #transform=ax1.transAxes,
                 ha='left',
                 va='bottom' )

   if x == 18:
      ax.text( x+3.2, y1-1.2, 'D1 DA cycling',
               fontsize=12,
               zorder=10,
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
               arrowprops=dict(facecolor=c2, width=5.0,
                               edgecolor=c2,
                               headwidth=10.0,headlength=7.0,
                               ) )

   if x == 0:
    
      x_ = x + dx24/24*6 + 2
      y_ = y2 - dy24/24*6 #-4

      #x0_ax, y0_ax = ax.transData.transform( (x, y2 ))
      #x1_ax, y1_ax = ax.transData.transform( (x+dx, y2+dy ))


#      axis_to_data = ax.transAxes + ax.transData.inverted()
#      data_to_axis = axis_to_data.inverted()
#
##      inv = ax.transData.inverted()
#      x0_ax, y0_ax = axis_to_data.transform( (x,y2) )
#      x1_ax, y1_ax = axis_to_data.transform( (x+dx,y2+dy) )
#      dx_ = x1_ax - x0_ax
#      dy_ = y1_ax - y0_ax
#      print( x0_ax, x1_ax )
#      print( y0_ax, y1_ax )
      angle = np.rad2deg( np.arctan( dy/ dx * xsize/ysize ) ) 
      print( angle )
      ax.text( x_, y_, 'D1 & D2 1-day extended forecast',
               color=c2,
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
                               edgecolor=c3,
                               shrink=0.00)
              )

   if x == 6:
    
      x_ = x_ - dx24/24*9 -1.5
      y_ = y3 + dy24/24*3 +0.5
      angle = np.rad2deg( np.arctan( dy/dx * xsize/ysize ) )
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
   if x >= -18:
    
      x_ = x_ + 3 #- dx24/24*9
      y_ = y3 - 2.0 #+ dy24/24*3
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

    
