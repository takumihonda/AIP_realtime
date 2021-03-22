import numpy as np
from matplotlib import pyplot as plt

fig, ax = plt.subplots( 1, 1, figsize=(13,4) )
#ax.plot(x,np.sin(x))
#ax.set_xlabel('x')
#ax.set_ylabel('y')

xmin = -1
xmax = 19

xticks = np.arange( xmin, xmax+1, 1 )

ymin = -0.1
ymax = 1

ax.arrow( x=0, y=0, dx=xmax, dy=0,
          width=0.02, head_width=0.05, head_length=0.3, 
          length_includes_head=True, color='k' )

ax.set_xlim( xmin, xmax )
ax.set_ylim( ymin, ymax )
ax.set_xticks( xticks )
#ax.set_xlimit( xmin, xmax )

ax.set_frame_on( False )
ax.axes.get_yaxis().set_visible( False )
plt.show()
