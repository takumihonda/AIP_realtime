
import numpy as np


def read( fn  ):

   try:
      f = open(fn)
   except:
      print("Failed to open")
      sys.exit()
      
   lines = f.read().split('\n')

   data = []

   for l in range(len(lines)):
      if lines[l] is '':
         break
      print(lines[l])
      data.append( float(lines[l]) )

   f.close()

   return( data )

# NOMELT
fn = "/data15/honda/SCALE-LETKF/AIP_SAFE/DEBUG256p/D4_500M_mercator_50mem_NOMELT/exp/2833233_cycle_20190610080000/RMSE.txt"
rmse_nomelt = read( fn )

fn = "/data15/honda/SCALE-LETKF/AIP_SAFE/DEBUG256p/D4_500M_mercator_50mem_NOMELT/exp/2833233_cycle_20190610080000/BIAS.txt"
bias_nomelt = read( fn )

# MELT

fn = "/data15/honda/SCALE-LETKF/AIP_SAFE/DEBUG256p/D4_500M_mercator_50mem_MELT/exp/2833286_cycle_20190610080000/RMSE.txt"
rmse_melt = read( fn )

fn = "/data15/honda/SCALE-LETKF/AIP_SAFE/DEBUG256p/D4_500M_mercator_50mem_MELT/exp/2833286_cycle_20190610080000/BIAS.txt"
bias_melt = read( fn )

print( len(rmse_nomelt) )

import matplotlib.pyplot as plt


x = np.arange(1, len(rmse_melt)+1)
x[0::2] = np.arange(len(rmse_melt)/2) + 1
x[1::2] = np.arange(len(rmse_melt)/2) + 1
print( x )

fig, (ax) = plt.subplots(1, 1, figsize=(5.0, 5.5))
fig.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.92, wspace=0.2, hspace=0.0)

ax.plot( x, rmse_melt, color='k', label='RMSE (DEFAULT)' )
ax.plot( x, rmse_nomelt, color='r', label='RMSE (NOMELT)' )
ax.plot( x, bias_melt, color='k', label='BIAS (DEFAULT)', ls='dashed' )
ax.plot( x, bias_nomelt, color='r', label='BIAS (NOMELT)', ls='dashed' )

ax.set_ylim( -10, 15 )
ax.set_xlim( 0, int( len(rmse_melt)/2 ) )
ax.set_xticks( np.arange( 0, int( len(rmse_melt)/2 ), 2 ) )
ax.set_yticks( np.arange( -10, 20, 2 ) )

ax.set_xlabel("# of DA cycles")

ax.legend()
#plt.legend8)

fig.suptitle( "RMSE/BIAS in reflectivity (dBZ)", fontsize=12)

plt.show()


