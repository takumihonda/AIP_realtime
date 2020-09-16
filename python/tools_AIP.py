import numpy as np

import os
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from datetime import timedelta

def read_CZ( ):
    fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/topo/init.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" )
    cz = nc.variables['CZ'][HALO:-HALO]
    nc.close()

    return( cz )

def read_nc_lonlat( fcst_zmax=43, obsz=np.arange(2) ):

    #fn = "/work/jh200062/share/honda/SCALE-LETKF/monitor/topo/topo.d4.nc"
    fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/topo/topo.d4.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" )
    lon = nc.variables['lon'][HALO:-HALO,HALO:-HALO]
    lat = nc.variables['lat'][HALO:-HALO,HALO:-HALO]
    topo = nc.variables['TOPO'][HALO:-HALO,HALO:-HALO]
    nc.close()

    cz = read_CZ( )

    hgt3d = np.zeros( ( fcst_zmax, topo.shape[0],topo.shape[1]) )
    for k in range( fcst_zmax ):
        hgt3d[k,:,:] = topo[:,:] + cz[k]

    ohgt3d = np.zeros( ( obsz.shape[0], topo.shape[0],topo.shape[1]) )
    for k in range( obsz.shape[0] ):
        ohgt3d[k,:,:] = topo[:,:] + cz[k]

    return( lon, lat, hgt3d, cz, ohgt3d )

def get_JMA_lonlat():
    JMAr_gx = 2560
    JMAr_gy = 3360

    dlat = 0.008333
    dlon = 0.012500

    latmin = 20.004167
    lonmin = 118.006250

    latmax =latmin + dlat * (JMAr_gy-1)
    jlat = np.arange( latmin, latmax, dlat )
    lonmax = lonmin + dlon * (JMAr_gx-1)
    jlon = np.arange( lonmin, lonmax, dlon )

    jlon2d, jlat2d = np.meshgrid( jlon, jlat )
    return( jlon2d, jlat2d )

def get_grads_JMA_snap( itime, ):
    print("get_JMA",itime)

    JMAr_gx = 2560
    JMAr_gy = 3360
    JMAr_top = "/data11/honda/SCALE-LETKF/JMA/radar/grads"

    latmax = 20.004167 + 0.008333 * (JMAr_gy-1)
    jlat = np.arange(20.004167, latmax,0.008333)
    lonmax = 118.006250 + 0.012500 * (JMAr_gx-1)
    jlon = np.arange(118.006250, lonmax,0.012500)

    fn = os.path.join( JMAr_top, itime.strftime('%Y/%m/%d'),
              "jmaradar_" + itime.strftime('%Y%m%d%H%M%S') + ".grd")

    try:
       infile = open(fn)
    except:
      print("Failed to open",fn)
      sys.exit()

    tmp2d = np.fromfile(infile, dtype=np.dtype('<f4'), count=JMAr_gx*JMAr_gy)  # little endian
    input2d = np.reshape(tmp2d, (JMAr_gy,JMAr_gx))
    input2d[input2d < 0.0] = 0.0 #np.nan # undef is negative in JMA radar
    infile.close

    return( input2d )


def get_Him8_FLDK(itime,band):
    print("get_Him8_FLDK")

    BB = "B" + str(band).zfill(2)

    otop = "/data9/honda/himawari_real/HIMAWARI-8/HISD/Hsfd"
    fn = os.path.join( otop, itime.strftime('%Y%m'), itime.strftime('%d'),
                 itime.strftime('%Y%m%d%H00'), itime.strftime('%M'), BB, 
                 "DLON0.02_DLAT0.02_HS_H08_" + itime.strftime('%Y%m%d_%H%M_')  + BB + "_FLDK.nc")

    obsnc = Dataset(fn,'r',format='NETCDF4')
    olon = obsnc.variables['longitude'][:]
    olat = obsnc.variables['latitude'][:]
    otbb = obsnc.variables['tbb'][:,:]

    obsnc.close()

    lon2d, lat2d = np.meshgrid( olon, olat )

    return( otbb, lon2d, lat2d )

def get_grads_JMA( itime, FT=1, ACUM=True ):
    print("get_JMA",itime)

    JMAr_gx = 2560
    JMAr_gy = 3360
    JMAr_top = "/data11/honda/SCALE-LETKF/JMA/radar/grads"

    jetime = itime + timedelta(seconds=3600) * FT

    jstime = itime + timedelta(seconds=600)
# comment out 9/14/2020
#    if not ACUM:
#      jstime += timedelta(seconds=3600) * (FT - 1)

    print(jstime,jetime)

    JMA_data_npz = os.path.join( JMAr_top, itime.strftime('%Y/%m/%d'),
             "jmaradar_" + jstime.strftime('%Y%m%d%H%M%S') + "-" + jetime.strftime('%Y%m%d%H%M%S') +  ".npz" )
    JMA_data_nc = os.path.join( JMAr_top, itime.strftime('%Y/%m/%d'),
            "jmaradar_" + jstime.strftime('%Y%m%d%H%M%S') + "-" + jetime.strftime('%Y%m%d%H%M%S') +  ".nc" )

    latmax = 20.004167 + 0.008333 * (JMAr_gy-1)
    JMAr_lat = np.arange(20.004167, latmax,0.008333)
    lonmax = 118.006250 + 0.012500 * (JMAr_gx-1)
    JMAr_lon = np.arange(118.006250, lonmax,0.012500)

    # npz is available?
    ISNPZ = os.path.isfile(JMA_data_npz)

    # nc is available?
    ISNC = os.path.isfile(JMA_data_nc)

    if ISNC:
      JMArain2d,JMAr_lon,JMAr_lat =  read_JMA_nc(JMA_data_nc)
      lon2d, lat2d = np.meshgrid( JMAr_lon, JMAr_lat )
      return( JMArain2d, lon2d, lat2d )

    if ISNPZ and not ISNC:
      JMArain2d = np.load(JMA_data_npz)['arr_0']
      print("read JMA data from: ", JMA_data_npz)

      write_JMA_nc(JMA_data_nc,JMAr_lon,JMAr_lat,jstime,JMArain2d)

    elif not ISNPZ and not ISNC:

      JMArain2d = np.zeros((JMAr_gy,JMAr_gx))

      jtime2 = jstime
      while (jtime2 <= jetime):

        jtime1 = jtime2 - timedelta(seconds=600)
        fn1 = os.path.join( JMAr_top, jtime1.strftime('%Y/%m/%d'),
                  "jmaradar_" + jtime1.strftime('%Y%m%d%H%M%S') + ".grd")
        fn2 = os.path.join( JMAr_top, jtime2.strftime('%Y/%m/%d'),
                  "jmaradar_" + jtime2.strftime('%Y%m%d%H%M%S') + ".grd")

        for fn in fn1,fn2:
           try:
             infile = open(fn)
           except:
             print("Failed to open",fn)
             sys.exit()

           tmp2d = np.fromfile(infile, dtype=np.dtype('<f4'), count=JMAr_gx*JMAr_gy)  # little endian
           input2d = np.reshape(tmp2d, (JMAr_gy,JMAr_gx))
           input2d[input2d < 0.0] = 0.0 #np.nan # undef is negative in JMA radar
           infile.close
           if fn == fn1:
              rain2d = input2d * 0.5
           else:
              rain2d += input2d * 0.5

        JMArain2d += rain2d / 6 # mm/h -> mm/10min
        jtime2 += timedelta(seconds=600)

      write_JMA_nc( JMA_data_nc, JMAr_lon, JMAr_lat, jstime, JMArain2d )
      print("New file stored: ", JMA_data_nc)
      #sys.exit()
    lon2d, lat2d = np.meshgrid( JMAr_lon, JMAr_lat )
    return( JMArain2d, lon2d, lat2d )

def write_JMA_nc(JMA_data_nc,JMA_lon,JMA_lat,jstime,JMArain2d):
      
    nc = Dataset(JMA_data_nc, "w", format="NETCDF4")
    
    nc.createDimension("latitude", len(JMA_lat))
    nc.createDimension("longitude", len(JMA_lon))
    nc.createDimension("level", 1)
    nc.createDimension("time", 1)
    
    XX = nc.createVariable("longitude","f4",("longitude",))
    XX.units = "degrees_east"
    
    YY = nc.createVariable("latitude","f4",("latitude",))
    YY.units = "degrees_north"
    
    ZZ = nc.createVariable("level","i4",("level",))
    ZZ.units = "km"
    
    times = nc.createVariable("time","f4",("time",))
    nc.description = "Superobs of Himawari-8 IR"
    
    times.units = "hours since " + str(jstime)
    times.calendar = "gregorian"
    times[0] = 0
    
    XVAR = nc.createVariable("RAIN","f4",("time","level","latitude","longitude"))
    XVAR.units = "mm"
    
    XVAR[:,:,:,:] = JMArain2d
    XX[:] = JMA_lon
    YY[:] = JMA_lat
    ZZ[:] = 0

def read_JMA_nc(JMA_data_nc):

    print(JMA_data_nc)
    nc = Dataset(JMA_data_nc, "r", format="NETCDF4")
    JMArain2d = nc.variables['RAIN'][0,0,:,:]
    JMA_lon = nc.variables['longitude'][:]
    JMA_lat = nc.variables['latitude'][:]

    nc.close()

    return( JMArain2d, JMA_lon, JMA_lat )

def get_GFS_grads(itime,var,zdim):
    print("get_GFS_grads",itime)

    GFS_top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/GFS/grads/20190824120000/mean"

    GFS_hPa = [ 1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
                 600, 550, 500, 450, 400, 350, 300, 250, 200, 150,
                 100,  70,  50,  30,  20,  10,   7,   5,  3,   2,
                   1]


    # surface
    sfc_fn = os.path.join(GFS_top, "sfc_" + itime.strftime('%Y%m%d%H%M%S') + ".grd")

    # atmos
    atm_fn = os.path.join(GFS_top, "atm_" + itime.strftime('%Y%m%d%H%M%S') + ".grd")

    GFS_gx = 361
    GFS_gy = 241
    GFS_lon = np.arange(  90.0, 180.25, 0.25 )
    GFS_lat = np.arange( 10.0, 70.25, 0.25 )
    GFS_lon2d, GFS_lat2d = np.meshgrid(GFS_lon,GFS_lat)


    count = GFS_gx * GFS_gy
    zmax = 31

    # variable dependent setting
    if zdim >= 1:
       fn = atm_fn
    else:
       fn = sfc_fn
       zdim = 1

    if var == "MSLETmsl":
       nv = 1
    if var == "UGRDprs":
       nv = 2
    if var == "VGRDprs":
       nv = 3
    if var == "TMPprs":
       nv = 4
    if var == "RHprs":
       nv = 5

    if fn == sfc_fn and (var == "UGRDprs" or var == "VGRDprs"):
       nv += 1

    rec = count * (zdim - 1) + count * zmax * (nv - 1)


    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       sys.exit()

    infile.seek(rec*4)
    tmp2d = np.fromfile(infile, dtype=np.dtype('>f4'), count=count)  # big endian
    input2d = np.reshape(tmp2d, (GFS_gy,GFS_gx))
    input2d[input2d > 9.999e19] = np.nan #np.nan
    infile.close

    return( input2d, GFS_lon2d, GFS_lat2d )

def read_nc( dom=1 ):

    top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/domains"

    fn = os.path.join( top, "topo.d" + str(dom) + ".nc")

    nc = Dataset(fn, "r", format="NETCDF4")
    lon2d = nc.variables["lon"][:]
    lat2d = nc.variables["lat"][:]
    topo2d = nc.variables["TOPO"][:]
    nc.close()

    return( lon2d, lat2d, topo2d )



def prep_proj_multi( METHOD, axs, res="c", ll_lon=120, ur_lon=150, ll_lat=20, ur_lat=50,
                     fs=10,zorder=2,contc='burlywood',cont_alp=0.2, lw=0.1, pdlon=0.1, pdlat=0.1 ):


    lat2 = 40.0

    blon = 139.609
    blat = 35.861

    m_l = []
    cc = "k"
    #contc = 'lightgrey'
    contc = contc


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


def read_dat( fn ):
    print( fn ) 

    recl = 11

    head = ("head",">i")
    tail = ("tail",">i")

    dtyp = np.dtype( [head, 
                      ("elm",">f"), 
                      ("lon",">f"), 
                      ("lat",">f"), 
                      ("lev",">f"), 
                      ("dat",">f"), 
                      ("err",">f"), 
                      ("typ",">f"), 
                      ("dif",">f"), 
                      ("qc",">f"), 
                      ("ob",">f"), 
                      ("oa",">f"), 
                      tail] )

    chunk = np.fromfile( fn, dtype=dtyp, count=-1 )

    print(chunk.shape)

    return( chunk[:]["elm"], chunk[:]["lon"], chunk[:]["lat"], chunk[:]["lev"], chunk[:]["dat"], chunk[:]["ob"], chunk[:]["qc"], chunk[:]["oa"] )


def get_nbin( vmin, vmax, dat ):
   h = 3.5*np.std( dat, ddof=1 ) / np.power( len(dat), 1.0/3.0 )
   nbin = int( (vmax - vmin) / h )

   return( nbin )

def desroz_diag_R( dat_a, dat_b ):
    dat = ( dat_a - np.mean(dat_a) ) * ( dat_b - np.mean(dat_b) )

    return( np.sqrt( np.mean(dat) ) )
