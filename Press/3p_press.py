import numpy as np
from datetime import datetime, timedelta
import sys
import os

from netCDF4 import Dataset

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

quick = True
quick = False

#    print("proj")

def dbz2rain( dbz ):
    return( np.power( np.power( 10, dbz*0.1) / 200.0, 1.0/1.6 ) ) 

def prep_proj_multi( METHOD, axs, res="c", ll_lon=120, ur_lon=150, ll_lat=20, ur_lat=50,
                     fs=10,zorder=2,contc='burlywood',cont_alp=0.2, lw=0.1 ):


    lat2 = 40.0

    blon = 139.609
    blat = 35.861

    m_l = []
    cc = "k"
    #contc = 'lightgrey'
    contc = contc

    pdlat = 0.1
    pdlon = 0.1

    lc = 'k'
    fs = fs
    lw = lw

    for ax in axs:
      m = Basemap(projection=METHOD,resolution = res,
              llcrnrlon = ll_lon,llcrnrlat = ll_lat,
              urcrnrlon = ur_lon,urcrnrlat = ur_lat,
              lat_0 = blat, lat_1 = lat2,
              lat_2 = lat2, lon_0 = blon,
              ax = ax)
      m_l.append(m)

      m.drawcoastlines(linewidth = 0.5, color = cc, zorder=zorder)
      m.fillcontinents(color=contc,lake_color='w', zorder=0, alpha =cont_alp)
      m.drawparallels(np.arange(0,70,pdlat),labels=[1,0,0,0],fontsize=fs,color=lc,linewidth=lw)
      m.drawmeridians(np.arange(0,180,pdlon),labels=[0,0,0,1],fontsize=fs,color=lc,linewidth=lw)


    return( m_l )

def read_nowcast_hires( stime=datetime(2019,9,10,9), ft=timedelta(minutes=5) ):

    top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/JMA/data/jma_nowcast_highres"

    time = stime + ft

    fn = os.path.join( top, time.strftime('%Y%m%d%H%M%S.nc') )

    nc = Dataset( fn, "r", format="NETCDF4" )

    rain2d = nc.variables["rain"][:,:]
    lon1d = nc.variables["longitude"][:]
    lat1d = nc.variables["latitude"][:]

    nc.close()

    lon2d, lat2d = np.meshgrid( lon1d, lat1d )

    return( rain2d, lon2d, lat2d )

def read_nowcast( stime=datetime(2019,9,10,9), ft=timedelta(minutes=5) ):
  
    dts = 300 # sec
    fts = ft.total_seconds()
    ft0 = -1
 
    top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/JMA/data"
    fn = os.path.join( top, stime.strftime('jma_nowcast_%Y%m%d%H%M%S.dat') )

    if fts == 0.0:
       # FT0 is composite radar 
       # http://www.jmbsc.or.jp/jp/online/file/f-online30200.html  
       fn = os.path.join( top, stime.strftime('jma_radar_%Y%m%d%H%M%S.dat') )
       ft0 = 0

    print( "Nowcast: ", fn )
    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       sys.exit()


  
    # Nowcast does not include FT0
    tlev = int( fts / dts ) - ft0

    gx = 2560
    gy = 3360
    gt = 12
    rec3d = gx*gy*gt
    rec2d = gx*gy

    nv = 1
    rec = 0

    infile.seek( rec2d*4*tlev )
    tmp3d = np.fromfile( infile, dtype=np.dtype('>f4'), count=rec2d )  # big endian   
    input3d = np.reshape( tmp3d, (gy,gx) )

    lons = 118.006250
    lats = 20.004167
    dlon = 0.012500
    dlat = 0.008333
    lon1d = np.arange( lons, lons+dlon*gx, dlon )
    lat1d = np.arange( lats, lats+dlat*gy, dlat )
    
    dt = 5*60
    t1d = np.arange( dt, dt*(12+1), dt )

    lon2d, lat2d = np.meshgrid( lon1d, lat1d )

    return( input3d, lon2d, lat2d, t1d  )


def read_fcst( stime=datetime(2019, 8, 24, 15, 30, 0 ), ft=timedelta(minutes=5) ):

    top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/Press"

    fn = os.path.join( top, stime.strftime('%Y%m%d-%H%M%S.nc') )

    print( "Forecast: ", fn )

    dts = 300 # sec
    fts = ft.total_seconds()
  
    tlev = int( fts / dts )

    nc = Dataset( fn, "r", format="NETCDF4" )
    ref3d = nc.variables["Reflectivity"][tlev,:,:,:]
    lon1d = nc.variables["Longitude"][:]
    lat1d = nc.variables["Latitude"][:]
    z1d = nc.variables["Height"][:]
    fts = nc.variables["time"][:]

    nc.close()

    return( np.transpose( ref3d, axes=(2,0,1) ), lon1d, lat1d, z1d, fts )

def read_obs_grads_latlon( kmax=1, ):
    
    # tlev starts from "0" (not from "1")
    gx = 161
    gy = 161

    count = gx*gy
    
    lonlat = {"lon":None, "lat":None}
    
    for tvar in ["lon", "lat"]: 
        fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m/" + tvar + ".bin"

        try:
           infile = open(fn)
        except: 
           print("Failed to open")
           print( fn )
           sys.exit()
    
        tmp2d = np.fromfile(infile, dtype=np.dtype('<f8'), count=count)  # big endian   
        input2d = np.reshape( tmp2d, (gy,gx) )

        lonlat[tvar] = input2d
    
    dz = 500.0
    oz = np.arange( 0.0, dz*kmax, dz )
           
    return( lonlat["lon"], lonlat["lat"], oz )

def read_obs( utime=datetime(2019,9,3,2,0,0), ):

    jtime = utime + timedelta(hours=9)

    OBS_EXIST = False
    for sec in range( -15, 16, 1 ):
        jtime2 = jtime + timedelta( seconds=sec )
        fn = os.path.join("/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m/",
                          jtime2.strftime('%Y/%m/%d/%H/%M/%S'),
                          "rain_cart_0002.nc")

        if os.path.isfile( fn ):
           OBS_EXIST = True
           break

    if not OBS_EXIST:
       print( "Not found OBS ", fn)
       sys.exit()

    print( "Obs: ", fn )


    try:
       nc = Dataset( fn, "r", format="NETCDF4" )
    except:
       print("Failed to open")
       print( fn )
       sys.exit()

    obs = nc.variables['rain'][0,:,:,:]

    lon2d, lat2d, z1d = read_obs_grads_latlon( kmax=obs.shape[0] )

    return( obs, lon2d, lat2d, z1d )


def main( stime=datetime( 2019, 8, 24, 15, 30), ft=timedelta(minutes=5) ):

    height = 3000.0

    vtime = stime + ft

    ref3d, lon1d, lat1d, z1d, fts, = read_fcst( stime=stime, ft=ft )
    lon2d, lat2d = np.meshgrid( lon1d, lat1d )

    zidx = np.argmin( np.abs( z1d - height ) )

    obs3d, olon2d, olat2d, oz1d = read_obs( utime=vtime )
    ozidx = np.argmin( np.abs( oz1d - height ) )

    #now2d, nlon2d, nlat2d, nt1d = read_nowcast( stime=stime, ft=ft )
    now2d, nlon2d, nlat2d = read_nowcast_hires( stime=stime, ft=ft )


    fig, (( ax1,ax2,ax3 )) = plt.subplots( 1, 3, figsize=( 13, 4.5 ) )
    fig.subplots_adjust( left=0.05, bottom=0.01, right=0.96, top=0.95,
                         wspace=0.15, hspace=0.02)

    lons = np.min( olon2d )
    lone = np.max( olon2d )
    lats = np.min( olat2d )
    late = np.max( olat2d )


    ax_l = [ ax1, ax2, ax3, ]

    if quick:
       res = "c"
    else:
       res = "h"

    m_l = prep_proj_multi('merc', ax_l, fs=7, res=res,
                           ll_lon=lons, ur_lon=lone, ll_lat=lats, ur_lat=late, lw=0.0 )

    lon2d_l = [ 
                nlon2d, 
                lon2d, 
                olon2d,
                ]
    lat2d_l = [ 
                nlat2d, 
                lat2d, 
                olat2d,
                ]

    print( "obs", obs3d.shape )


    print( now2d.shape, nlon2d.shape, nlat2d.shape  )

    var_l = [ 
              now2d[:,:],
              dbz2rain( ref3d[zidx,:,:] ), 
              dbz2rain( obs3d[ozidx,:,:] ),
            ]

    fts = ft.total_seconds()

    tit_l = [ 
              "JMA Nowcast (FT={0:0=2}min)".format( int( fts/60 ) ),
              "SCALE-LETKF (FT={0:0=2}min)".format( int( fts/60 ) ),
              "PAWR Obs",
              ]


    cmap = plt.cm.get_cmap("jet")
    extend = 'max'

    levs = np.array( [ 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40 ] )

    for i, ax in enumerate( ax_l ):
        print( i )
        x2d, y2d = m_l[0](lon2d_l[i], lat2d_l[i])

        SHADE = ax.contourf( x2d, y2d, var_l[i],
                                levels=levs, 
                                cmap=cmap, 
                                extend=extend )

        ax.text( 0.5, 1.01, tit_l[i],
                fontsize=13, transform=ax.transAxes,
                horizontalalignment='center',
                verticalalignment='bottom', )
       

    pos = ax3.get_position()
    cb_width = 0.01
    cb_height = pos.height*1.0
    ax_cb = fig.add_axes( [pos.x1, pos.y0+0.01, cb_width, cb_height] )
    cb = plt.colorbar( SHADE, cax=ax_cb, orientation = 'vertical', ticks=levs[::1] )
    cb.ax.tick_params( labelsize=6 )

    ax3.text( 0.95, 1.01, '(mm/h)',
            fontsize=9, transform=ax.transAxes,
            horizontalalignment='left',
            verticalalignment='bottom', )
       

    tit = 'Valid: {0:}'.format( vtime.strftime('%H:%M:%S %m/%d/%Y') )

    fig.suptitle( tit, fontsize=16 )

    ofig = "report_3p_s{0:}_FT{1:0=2}m_hires".format( stime.strftime('%Y%m%d%H%M%S'), int( fts / 60 ) )
    print( ofig )

    if not quick:

       opath = "png/3p_press_hires"
       os.makedirs(opath, exist_ok=True)
 
       ofig = os.path.join(opath, ofig + ".pdf")
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       print(ofig)
       plt.clf()
    else:
       print(ofig)
       plt.show()



##############


stime = datetime( 2019, 8, 24, 15, 30, 0 ) 

ft = 10
while ft <= 10:
   ft_ = timedelta( minutes=ft )
   main( stime=stime, ft=ft_ )

   ft += 5

