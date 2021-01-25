import numpy as np
import sys

range_ = 10.e3
height_ = 11.e3
elv = np.rad2deg( np.arctan2( height_, range_ ) )
print( elv )

print( "beam interval" )
da = 360 / 300 # delta azimith # degree
dist = np.arange( 20.0, 70.0, 10.0 ) # km
print( dist )
print( dist * np.sin( np.deg2rad( da )) )


sys.exit()

dx = 0.5
lats = 30
late = 40
dlat = 0.25
 
lons = 130
lone = 140
dlon = 0.25

x1d = np.arange( lons, lone+dlon, dlon )
y1d = np.arange( lats, late+dlat, dlat )

x2d, y2d = np.meshgrid( x1d, y1d )

xc = np.mean( x1d )
yc = np.mean( y1d )

dist = np.sqrt( np.square( x2d - xc ) +  np.square( y2d - yc ) ) 

dlon = x2d - xc
dlat = y2d - yc
azm = np.rad2deg( np.arctan2( dlon * np.cos( np.deg2rad(yc) ), dlat ) )
azm[ azm < 0.0] += 360.0

u2d = np.zeros( azm.shape ) 
v2d = np.zeros( azm.shape ) + 20.0
# assume elevation angle = 0.0
elv = 10.0
elv = 0.0
vr2d = u2d * np.cos( np.deg2rad(elv) ) * np.sin( np.deg2rad(azm) ) + \
       v2d * np.cos( np.deg2rad(elv) )  * np.cos( np.deg2rad(azm) )


# plot
import matplotlib.pyplot as plt

fig, (ax1 ) = plt.subplots(1, 1, figsize=(7,5))
#ax1.set_aspect('equal')

#ax1.contourf(x2d, y2d, dist )
cmap = plt.cm.get_cmap( "hot" )
levs = np.arange(0, 360, 1 )
var2d = azm

cmap = plt.cm.get_cmap( "RdBu_r" )
levs = np.arange( -20, 21, 1 )

var2d = vr2d

SHADE = ax1.contourf(x2d, y2d, var2d, cmap=cmap, levels=levs )

ax1.set_xlabel( "Lon" )
ax1.set_ylabel( "Lat" )


plt.colorbar( SHADE )



plt.show()


