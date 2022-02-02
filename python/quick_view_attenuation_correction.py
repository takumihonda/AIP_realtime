
import numpy as np

fn = '/data_ballantine02/miyoshi-t/otsuka/nowcast_pawr/saitama/obs/500mKdpR/2019/08/25/00/00/14/corrected_zh_polar_00000.bin'

#800         301         114

na = 301
nr = 800
ne = 114

infile = open( fn )
tmp = np.fromfile( infile, dtype=np.dtype('<f4'), count=na*nr*ne )
z3d = np.reshape( tmp, (nr,na,ne) )

import matplotlib.pyplot as plt
print( z3d.shape )
plt.plot( z3d[:,10,20] )
plt.show()

