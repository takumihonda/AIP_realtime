from netCDF4 import Dataset
import matplotlib.pyplot as plt

from tools_AIP import prep_proj_multi

quick = True

def read_nc( prj="MER" ):

    fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/DEBUG_{0:}/const/topo_sno_np00001/topo.pe000000.nc".format( prj )

    nc = Dataset(fn, "r", format="NETCDF4")
    lon2d = nc.variables["lon"][:]
    lat2d = nc.variables["lat"][:]
    topo2d = nc.variables["topo"][:]
    nc.close()

    return( lon2d, lat2d, topo2d )

def main():
    mmlon2d, mmlat2d, mmtopo2d = read_nc( prj="MER_MLAT" )
    mlon2d, mlat2d, mtopo2d = read_nc( prj="MER" )
    llon2d, llat2d, ltopo2d = read_nc( prj="LAM" )

    fig, ( ax1,ax2,ax3 ) = plt.subplots( 1, 3, figsize=( 11, 4.0 ) )
    fig.subplots_adjust( left=0.04, bottom=0.03, right=0.97, top=0.97,
                         wspace=0.15, hspace=0.05)

    ax_l = [ ax1, ax2, ax3 ]

    if quick:
       res = "l"
    else:
       res = "f"

    lons = 138.7
    lone = 140.5
    lats = 35.2
    late = 36.5
    method = "merc"
    lon_r = 139.609
    lat_r = 35.861
    contc = "palegreen"
    contc = "burlywood"
    oc = "w"
    lon_0 = lon_r
    lat_0 = lat_r
    if quick:
       res = 'l'
    else:
       res = 'f'

    m_l = prep_proj_multi( method, ax_l, fs=6, res=res, lw=0.0, 
                           ll_lon=lons, ur_lon=lone, ll_lat=lats, ur_lat=late, 
                           pdlon=0.2, pdlat=0.2, blon=lon_r, blat=lat_0,
                           contc=contc, oc=oc, )

    lon2d_l = [ mmlon2d, mlon2d, llon2d ]
    lat2d_l = [ mmlat2d, mlat2d, llat2d ]

    var2d_l = [ mmtopo2d, mtopo2d, ltopo2d ]
    tit_l = [ "MER(35deg)", "MER(0deg)", "LAMBERT",]

    for i, ax in enumerate( ax_l ):
        x2d, y2d = m_l[0]( lon2d_l[i], lat2d_l[i] )

        SHADE = ax.pcolormesh( x2d, y2d, var2d_l[i], )

        ax.text( 0.5, 1.01, tit_l[i],
                va='bottom',
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=10, )
        
    plt.show()

main()
