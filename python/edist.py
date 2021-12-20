import numpy as np


# define 3d grid
xmin = -10
xmax = 10
dx = 1
x1d = np.arange( xmin, xmax+dx, dx ) * 500.0

ymin = -10
ymax = 10
dy = 1
y1d = np.arange( ymin, ymax+dy, dy ) * 500.0

zmin = -10
zmax = 10
dz = 1
z1d = np.arange( zmin, zmax+dz, dz ) * 500.0


x3d, y3d, z3d = np.meshgrid( x1d, y1d, z1d )

# thinning interval
h = 4
v = 4
dist_l = np.arange( 1000.0, 4000.0, 25.0 )

h = 8
v = 8
dist_l = np.arange( 1000.0, 8000.0, 25.0 )

h = 1
v = 1
dist_l = np.arange( 1000.0, 8000.0, 25.0 )

# grid location ( should be positive for all 3 components )
cx, cy, cz = 0.1, 0.1, 0.1

dist3d = np.sqrt( ( x3d - cx )**2 + ( y3d - cy )**2 + ( z3d - cz )**2 )

flag3d = np.zeros( dist3d.shape )

# nearest point
ix, iy, iz = np.unravel_index( np.argmin( dist3d ), dist3d.shape ) 

# chose nearest eight points
#print( dist3d[ix,iy,iz], dist3d[ix+1,iy,iz], dist3d[ix-1,iy,iz])
flag3d[ix:ix+2,iy:iy+2,iz:iz+2] = 1.0
print( np.sum( flag3d ) )

# thinning

sx = ix % h
sy = iy % h
sz = iz % v

flag3d[ix+1+h::h,:,:] = 1.0
flag3d[sx:ix:h,:,:] = 1.0

flag3d[:,iy+1+h::h,:] = 1.0
flag3d[:,sy:iy:h,:] = 1.0

flag3d[:,:,iz+1+v::v] = 1.0
flag3d[:,:,sz:iz:v] = 1.0

print( flag3d[:,iy,iz])
print( dist3d.shape )
print( np.sum( flag3d ))

for dist in dist_l:
    num = np.sum( np.where( ( dist3d < dist ) & ( flag3d > 0 ), 1, 0  ) )
    #num = np.sum( np.where( ( dist3d < dist ), 1, 0  ) )
#    num = np.where( ( dist3d < dist )&( flag3d > 0 ), 1, 0  )
    print( dist, num  )
#    print( np.sum( num ), num.shape )

    if num > 100:
       break

