import numpy as np
import sys
import os

from netCDF4 import Dataset

zlev = 14
HALO = 2

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/DEBUG20200827/20200827090000"

quick = False
#quick = True

def read_cz( fn="" ):
    nc = Dataset( fn, "r", format="NETCDF4" )
    var = nc.variables["CZ"][HALO+zlev]
    nc.close()

    return( var )

def read_nc3d( fn="", nvar="QR" ):
    nc = Dataset( fn, "r", format="NETCDF4" )
    var = nc.variables[nvar][:]
    #var = nc.variables[nvar][:] #[HALO:-HALO:2,HALO:-HALO:2,zmax]
    nc.close()

    return( var[HALO:-HALO:2,HALO:-HALO:2,zlev] )
   

def get_dbz( m=1 ):
    mm = str(m).zfill(4)
    if m == 0:
       mm = "mean"
    fn = os.path.join( TOP, "anal_sno_np00001/{0:}/anal.pe000000.nc".format( mm ) )
    print( fn )
    cz = read_cz( fn=fn )
    qr = read_nc3d( fn=fn, nvar="QR" )
    qs = read_nc3d( fn=fn, nvar="QS" )
    qg = read_nc3d( fn=fn, nvar="QG" )
    rho = read_nc3d( fn=fn, nvar="DENS" )
  
    dbz =  q2dbz( qr, qs, qg, rho )
    return( dbz )

def q2dbz( qr, qs, qg, rho ):

    qrp = qr
    qsp = qs
    qgp = qg

    zr = 2.53e4 * np.power( rho*qrp*1.e3, 1.84 )
    zs = 3.48e3 * np.power( rho*qsp*1.e3, 1.66 )
    zg = 5.54e3 * np.power( rho*qgp*1.e3, 1.7 )

    ref = zr + zs + zg

    return( np.where( ref > 0.0, 10*np.log10( zr + zs + zg ), 0.0 ) )

def main():

    fn = os.path.join( TOP, "anal_sno_np00001/{0:}/anal.pe000000.nc".format( "mean" ) )
    cz = read_cz( fn=fn )
    dbz = get_dbz( m=0 )
    dbz[:] = 0.0 

    thrs = 10.0
 
    levs = np.arange( 1, 6, 1 )
 
    for m in range( 1, 51 ):
       dbz_ = get_dbz( m=m )
       dbz += np.where( dbz_ > thrs, 1, 0 )

    import matplotlib.pyplot as plt

    fig, (ax1) = plt.subplots(1, 1, figsize=(6,6) )
    fig.subplots_adjust(left=0.07, bottom=0.07, right=0.93, top=0.93, )

    ax1.set_aspect('equal')

    cmap = plt.cm.get_cmap("jet")
    SHADE = ax1.contourf( dbz, levels=levs, cmap=cmap,
                          extend='max' )



    pos = ax1.get_position()
    cb_h = pos.height
    ax_cb = fig.add_axes( [pos.x1+0.005, pos.y0, 0.01, cb_h] )
    cb = plt.colorbar( SHADE, cax=ax_cb, )

    tit = "Rainy member, Z={0:.1f}km, {1:}dBZ".format( cz/1000.0, thrs )
    fig.suptitle( tit, fontsize=14 ) 

    ofig = "member_cnt"

    if quick:
       plt.show()
    else:
       odir = "png/2d_map_all"
       os.makedirs( odir, exist_ok=True)
       plt.savefig( os.path.join(odir, ofig),
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')

main() 
