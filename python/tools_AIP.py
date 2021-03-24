import numpy as np

import os
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


from datetime import timedelta, datetime

def read_CZ( ):
    fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/topo/init.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" )
    cz = nc.variables['CZ'][HALO:-HALO]
    nc.close()

    return( cz )

def read_nc_lsmask( NEW=False ):

    if NEW:
       fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/TEST_INPUT/D4_500m_20201117/const/topo_sno_np00001/topo.pe000000.nc"
    else:
       fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/topo/topo.d4.nc"

    HALO = 2
    nc = Dataset( fn, "r", format="NETCDF4" )
    lsmask = nc.variables['lsmask'][HALO:-HALO,HALO:-HALO]
    nc.close()

    return( lsmask )

def read_nc_lonlat( fcst_zmax=43, obsz=np.arange(2), NEW=False ):

    #fn = "/work/jh200062/share/honda/SCALE-LETKF/monitor/topo/topo.d4.nc"
    if NEW:
       fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/TEST_INPUT/D4_500m_20201117/const/topo_sno_np00001/topo.pe000000.nc"
    else:
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

def get_GFS_grads( itime, var, zdim ):
    print("get_GFS_grads",itime)

    GFS_top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/GFS/grads/{0:}/mean".format( itime.strftime('%Y%m%d%H%M%S')  )

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

    top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/domains_20201117"

    fn = os.path.join( top, "topo.d" + str(dom) + ".nc")

    nc = Dataset(fn, "r", format="NETCDF4")
    lon2d = nc.variables["lon"][:]
    lat2d = nc.variables["lat"][:]
    topo2d = nc.variables["TOPO"][:]
    nc.close()

    return( lon2d, lat2d, topo2d )



def prep_proj_multi( METHOD, axs, res="c", ll_lon=120, ur_lon=150, ll_lat=20, ur_lat=50,
                     fs=10,zorder=2,contc='burlywood',cont_alp=0.2, lw=0.1, pdlon=0.1, pdlat=0.1, 
                     blon=139.609, blat=35.861, oc="aqua" ):


    lat2 = 40.0

    m_l = []
    cc = "k"


    lc = 'k'
    fs = fs
    lw = lw

    for ax in axs:
      m = Basemap(projection=METHOD,resolution=res,
              llcrnrlon = ll_lon,llcrnrlat = ll_lat,
              urcrnrlon = ur_lon,urcrnrlat = ur_lat,
              lat_0 = blat, lat_1 = lat2,
              lat_2 = lat2, lon_0 = blon,
              ax = ax)
      m_l.append(m)

      m.drawcoastlines( linewidth = 0.5, color=cc, zorder=zorder)
#      m.fillcontinents(color=contc,lake_color='w', zorder=0, alpha=cont_alp)
      m.drawlsmask( land_color=contc, ocean_color=oc, lakes=True, 
                    alpha=cont_alp,
                    resolution=res )
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

def read_obs( utime=datetime(2019,9,3,2,0,0), mask=np.array([]) ):

    jtime = utime + timedelta(hours=9)
    jsec = jtime.second

    if jsec == 30:
       itime = datetime( jtime.year, jtime.month, jtime.day, jtime.hour, jtime.minute, 15 )
    elif jsec == 0:
       itime = datetime( jtime.year, jtime.month, jtime.day, jtime.hour, jtime.minute, 45 ) - timedelta( minutes=1 )


    OBS_EXIST = False
    for sec in range( 0, 30, 1 ):
        jtime2 = itime + timedelta( seconds=sec )
        fn = os.path.join("/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_new/",
                          jtime2.strftime('%Y/%m/%d/%H/%M/%S'),
                          "rain_cart_0002.nc")

        if os.path.isfile( fn ):
           OBS_EXIST = True
           break

    if not OBS_EXIST:
       print( "Not found OBS ", fn, jtime)
       return( None, False )

    try:
       nc = Dataset( fn, "r", format="NETCDF4" )
    except:
       print("Failed to open")
       print( fn )
       sys.exit()

    obs = nc.variables['rain'][0,:,:,:]
    nc.close()

    obs = np.where( mask > 0.0, np.nan, obs )

    return( obs, True )

def read_obs_grads( INFO, itime=datetime(2019,9,10,9), ores='500m' ):
    fn = os.path.join( INFO["TOP"], INFO["EXP"],
                       INFO["time0"].strftime('%Y%m%d%H0000'),  "pawr_grads/pawr_ref3d_" +
                       itime.strftime('%Y%m%d-%H%M%S.grd') )
    gz = 22
    dz = 500.0
    gx = 241
    gy = 241
    dlon = 0.00554812
    dlat = 0.00449640

    if ores == "100m":
       fn = os.path.join( INFO["TOP"], "20201117/D4_500m_NODA_O100M",
                          INFO["time0"].strftime('%Y%m%d%H0000'),  "pawr_grads/pawr_ref3d_" +
                          itime.strftime('%Y%m%d-%H%M%S.grd') )
       gz = 110
       dz = 100.0
       gx = 1201
       gy = 1201
       dlon = 0.00110962
       dlat = 0.00089928


    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       sys.exit()

    rec3d = gx*gy*gz

    nv = 1
    rec = 0

    infile.seek(rec*4)
    tmp3d = np.fromfile(infile, dtype=np.dtype('>f4'), count=rec3d)  # big endian   
    input3d = np.reshape( tmp3d, (gz,gy,gx) )

    lons = 138.94319715
    lats = 35.32193244
    lon1d = np.arange( lons, lons+dlon*gx, dlon )
    lat1d = np.arange( lats, lats+dlat*gx, dlat )
    z1d = np.arange( 0.0, gz*dz, dz )

    lon2d, lat2d = np.meshgrid( lon1d, lat1d )

#    gx = 241
#    gy = 241
##    gz = 241
#    for nvar in "lon", "lat":
#        fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_full/shadow/{0:}.bin".format( nvar )   
#        try:
#           infile = open(fn)
#        except:
#           print("Failed to open")
#           sys.exit()
#    
#        rec = gx * gy
#        tmp = np.fromfile( infile, dtype=np.dtype('<f8'), count=rec )
#        
#        if nvar == "lon":
#           lon2d = np.reshape( tmp, (gy,gx) )
#        elif nvar == "lat":
#           lat2d = np.reshape( tmp, (gy,gx) )
#



    return( input3d, lon2d, lat2d, z1d  )


def read_nc_topo( dom=1 ):

    top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/domains_20201117"

    fn = os.path.join( top, "topo.d" + str(dom) + ".nc") 

    nc = Dataset(fn, "r", format="NETCDF4")
    lon2d = nc.variables["lon"][:]
    lat2d = nc.variables["lat"][:]
    topo2d = nc.variables["TOPO"][:]
    nc.close()

    return( lon2d, lat2d, topo2d )

def read_mask_grads():
    fn = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/git/AIP_realtime/python/dat/saitama-shadow-mask.dat"

    print( fn )
    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       return( None, False )
       sys.exit()

    tmp = np.fromfile( infile, dtype=np.dtype('<i4'), count=2 ) 
    dims = np.reshape( tmp, (2) )
    gz = dims[0]
    ga = dims[1]
    print( dims )

    infile.seek(2*4)
    rec = gz*ga
    tmp2d = np.fromfile( infile, dtype=np.dtype('<i2'), count=rec )
    input2d = np.reshape( tmp2d, (gz,ga) )

    return( input2d )

def read_mask_full():
    fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_full/shadow/mask.nc"   
    nc = Dataset( fn, "r", format="NETCDF4" )
    mask = nc.variables['mask'][0,:,:,:]
       
    nc.close()


    gx = 241
    gy = 241
#    gz = 241
    for nvar in "lon", "lat":
        fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_full/shadow/{0:}.bin".format( nvar )   
        try:
           infile = open(fn)
        except:
           print("Failed to open")
           sys.exit()
    
        rec = gx * gy
        tmp = np.fromfile( infile, dtype=np.dtype('<f8'), count=rec )
        
        if nvar == "lon":
           lon2d = np.reshape( tmp, (gy,gx) )
        elif nvar == "lat":
           lat2d = np.reshape( tmp, (gy,gx) )


    return( mask, lon2d, lat2d )



def read_mask():
    fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_new/shadow/mask.nc"   
    nc = Dataset( fn, "r", format="NETCDF4" )
    mask = nc.variables['mask'][0,:,:,:]
       
    nc.close()

    return( mask )

def read_fcst_grads( INFO, itime=datetime(2019,9,3,2,0,0), tlev=0 , FT0=True, ):

    fn = os.path.join( INFO["FCST_DIR"],
                       itime.strftime('%Y%m%d-%H%M%S.grd') )

    print( fn )
    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       return( None, False )
       sys.exit()

    # tlev starts from "0" (not from "1")
    #gz = 43 # only below 15km is stored in fcst
    gz = INFO["gz"]
    gx = INFO["gx"]
    gy = INFO["gy"]

    rec3d = gx*gy*gz

    nv = 1

    if FT0:
      rec = rec3d * nv * tlev
    else:
      rec = rec3d * nv * ( tlev - 1 )

    try:
       infile.seek(rec*4)
       tmp3d = np.fromfile(infile, dtype=np.dtype('>f4'), count=rec3d)  # big endian   
       input3d = np.reshape( tmp3d, (gz,gy,gx) )
       fstat = True
    except:
       input3d = None
       fstat = False

    return( input3d, fstat )

def read_obs_grads_latlon( ores='500m' ):

    # tlev starts from "0" (not from "1")
    gx = 161
    gy = 161
    gz = 29

    count = gx*gy

    lonlat = {"lon":None, "lat":None}

    for tvar in ["lon", "lat"]:
        fn = "/lfs01/otsuka/_OLD_DATA12_/nowcast_pawr/saitama/obs/500m_new/" + tvar + ".bin"

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
    oz = np.arange( 0.0, dz*gz, dz )

    return( oz, lonlat["lon"], lonlat["lat"] )

def draw_rec( ax, m, lon2d, lat2d, lc='k', lw=1.0 ):

    x, y = m( lon2d[0,:], lat2d[0,:] )
    ax.plot( x, y, color=lc, lw=lw )

    x, y = m( lon2d[-1,:], lat2d[-1,:] )
    ax.plot( x, y, color=lc, lw=lw )

    x, y = m( lon2d[:,0], lat2d[:,0] )
    ax.plot( x, y, color=lc, lw=lw )

    x, y = m( lon2d[:,-1], lat2d[:,-1] )
    ax.plot( x, y, color=lc, lw=lw )

def dist( rlon=130.0, rlat=30.0, lons=np.zeros(1), lats=np.zeros(1) ):
    lon1 = rlon * np.pi / 180.0
    lat1 = rlat * np.pi / 180.0

    lon2 = lons * np.pi / 180.0
    lat2 = lats * np.pi / 180.0

    cosd = np.sin( lat1 ) * np.sin( lat2 ) + np.cos( lat1 ) * np.cos( lat2 ) * np.cos( lon2 - lon1 )

    cosd[ cosd > 1.0 ] = 1.0
    cosd[ cosd < -1.0 ] = -1.0

    dist = np.arccos( cosd ) * 6371.3e3
    return( dist )

def read_oerr_npz( range_inl=[1.0], elv_inl=[5.0], dr=1.0, de=1.0, mode="az",
                   azm_inl=[10.0], da=1.0, DR=500.0, exp="20201117/D4_500m_H1V1_Z", otyp=4002 ):

    amax = 10

    if mode == "az":
       azm_inl = [ 999 ]
       x1d = np.arange( -amax, amax+1 ) * da

    elif mode == "ra":
       range_inl = [ 999 ]
       x1d = np.arange( -amax, amax+1 ) * dr * DR

    elif mode == "el":
       elv_inl = [ 999 ]
       x1d = np.arange( -amax, amax+1 ) * de

    cnt1d = np.zeros( amax*2+1)
    cor1d = np.zeros( amax*2+1)
    oerr = 0.0

    cnt = 0
    for elv_in in elv_inl:
        for azm_in in azm_inl:
            for range_in in range_inl:
 
                dir_in = "dat/{0:}/{1:05}".format( exp, otyp )
                ofile = os.path.join( dir_in, "oerr_{0:}_e{1:03}_{2:03}_r{3:03}_{4:03}_a{5:03}_{6:03}_{7:03}.npz".format( mode, elv_in, de, range_in, dr, azm_in, da, int( DR ) ) )
 
                try:
                   #print( "Found npz file ", ofile )
                   data = np.load( ofile )
                   cor1d_ = data["cor1d"]
                   cnt1d_ = data["cnt1d"]
                   oerr_ = data["oerr"]
                except:
                   print( "No npz file ", ofile)
                   print( ofile )
                   sys.exit()

                if not np.isnan( oerr_ ) and not np.isnan( cor1d_[amax] ) :
                 
                   cor1d_[ ( cnt1d_ < 1 ) | np.isnan( cnt1d_ ) ] = 0.0
                   cnt1d_[ ( cnt1d_ < 1 ) | np.isnan( cnt1d_ ) ] = 0.0

                   cor1d += cor1d_ 
                   cnt1d += cnt1d_ 
                   oerr += oerr_
                   cnt += 1

    cor1d = cor1d / cnt1d
    cor1d = cor1d / cor1d[amax]
    if cnt > 0:
       oerr = oerr / cnt
   
    return( cor1d, cnt1d, oerr, x1d )

def read_fcst_grads_all( INFO, itime=datetime(2019,9,3,2,0,0), tlev=0 , FT0=True, nvar="p", ):

    fn = os.path.join( INFO["FCST_DIR"],
                       itime.strftime('fcst_all_%Y%m%d-%H%M%S.grd') )

    print( fn )
    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       return( None, False )
       sys.exit()

    # tlev starts from "0" (not from "1")
    gz = 60
    gx = INFO["gx"]
    gy = INFO["gy"]

    rec3d = gx*gy*gz
    rec2d = gx*gy

    # u, v, w, t, p, qv, qc, qr, qi, qs, qg
    nv3d = 11
    # PW PRCP
    nv2d = 2

    if FT0:
      rec = ( rec3d * nv3d + rec2d * nv2d ) * tlev
    else:
      rec = ( rec3d * nv3d + rec2d * nv2d ) * ( tlev - 1 )

    nv2d_ = 0
    if nvar == "u":
       nv3d_ = 0
    elif nvar == "v":
       nv3d_ = 1
    elif nvar == "w":
       nv3d_ = 2
    elif nvar == "t":
       nv3d_ = 3
    elif nvar == "p":
       nv3d_ = 4
    elif nvar == "qv":
       nv3d_ = 5
    elif nvar == "qc":
       nv3d_ = 6
    elif nvar == "qr":
       nv3d_ = 7
    elif nvar == "qi":
       nv3d_ = 8
    elif nvar == "qs":
       nv3d_ = 9
    elif nvar == "qg":
       nv3d_ = 10

    rec = rec + ( rec3d * nv3d_ + rec2d * nv2d_ )

    try:
       infile.seek( rec*4 )
       tmp3d = np.fromfile( infile, dtype=np.dtype('>f4'), count=rec3d )  # big endian   
       input3d = np.reshape( tmp3d, (gz,gy,gx) )
    except:
       input3d = None

    return( input3d )

def read_ga_grads_all( INFO, itime=datetime(2019,9,3,2,0,0), nvar="p", typ="g" ):
    # read guess/analsis ensemble mean  

    if typ == "g":
       fn_ = itime.strftime( 'gues%Y%m%d-%H%M%S.grd')
    elif typ == "a":
       fn_ = itime.strftime( 'anal%Y%m%d-%H%M%S.grd')

    fn = os.path.join( INFO["TOP"], INFO["EXP"],
                       INFO["time0"].strftime('%Y%m%d%H%M%S'),
                       "mean_grads", fn_ )

    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       return( None, False )
       sys.exit()

    gz = 60
    gx = INFO["gx"]
    gy = INFO["gy"]

    rec3d = gx*gy*gz
    rec2d = gx*gy

    # u, v, w, t, p, qv, qc, qr, qi, qs, qg
    nv3d = 11
    # PW PRCP
    nv2d = 0

    nv2d_ = 0
    nv3d_ = 0
    if nvar == "u":
       nv3d_ = 0
    elif nvar == "v":
       nv3d_ = 1
    elif nvar == "w":
       nv3d_ = 2
    elif nvar == "t":
       nv3d_ = 3
    elif nvar == "p":
       nv3d_ = 4
    elif nvar == "qv":
       nv3d_ = 5
    elif nvar == "qc":
       nv3d_ = 6
    elif nvar == "qr":
       nv3d_ = 7
    elif nvar == "qi":
       nv3d_ = 8
    elif nvar == "qs":
       nv3d_ = 9
    elif nvar == "qg":
       nv3d_ = 10

    rec = rec3d * nv3d_ + rec2d * nv2d_

    try:
       infile.seek( rec*4 )
       tmp3d = np.fromfile( infile, dtype=np.dtype('>f4'), count=rec3d )  # big endian   
       input3d = np.reshape( tmp3d, (gz,gy,gx) )
    except:
       input3d = None

    return( input3d )

def get_cfeature( typ='land', res='10m' ):

    if typ == 'land': 
       name = 'land'
       fcolor = cfeature.COLORS['land']
       ecolor = 'face'
    elif typ == 'coastline': 
       name = 'coastline'
       fcolor = 'none'
       ecolor = 'k'

    feature = cfeature.NaturalEarthFeature( 'physical', name, res, 
                                            edgecolor=ecolor,
                                            facecolor=fcolor )

    return( feature )

def setup_grids_cartopy( ax, xticks=np.array([]), yticks=np.array([]), lw=0.5, 
                         fs=10, fc='k', xfs=-1, yfs=-1 ):
       gl = ax.gridlines( crs=ccrs.PlateCarree(), linewidth=lw, 
                          draw_labels=True  )
       gl.xlabels_top = False
       gl.ylabels_right = False 
       gl.xlocator = mticker.FixedLocator( xticks ) 
       gl.ylocator = mticker.FixedLocator( yticks ) 
       if lw == 0.0:
          gl.xlines = False
          gl.ylines = False
       else:
          gl.xlines = True
          gl.ylines = True
       gl.xformatter = LONGITUDE_FORMATTER
       gl.yformatter = LATITUDE_FORMATTER

       if xfs < 0:
          xfs = fs
       if yfs < 0:
          yfs = fs

       gl.xlabel_style = {'size': xfs, 'color': fc, }
       gl.ylabel_style = {'size': yfs, 'color': fc, }

def prep_proj_multi_cartopy( fig, xfig=1, yfig=1, proj='none', latitude_true_scale=35.0 ):

    if proj == 'none':
       projection = ccrs.PlateCarree()

    elif proj == 'merc':
       projection = ccrs.Mercator( latitude_true_scale=latitude_true_scale, ) 

    ax_l = []
    for i in range( 1, xfig*yfig+1 ):
       ax_l.append( fig.add_subplot( yfig,xfig,i, projection=projection ) )

    return( ax_l )

def read_nowcast_hires( stime=datetime(2019,9,10,9), ft=timedelta(minutes=5) ):

    top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/JMA/data/jma_nowcast_highres"

    time = stime + ft

#    fn = os.path.join( top, time.strftime('%Y%m%d%H%M%S.nc') )
    fn = os.path.join( top, stime.strftime('%Y%m%d%H%M%S'),
              time.strftime('%Y%m%d%H%M%S.nc') )
    print( fn )
    print( "" )

    nc = Dataset( fn, "r", format="NETCDF4" )

    rain2d = nc.variables["rain"][:,:]
    lon1d = nc.variables["longitude"][:]
    lat1d = nc.variables["latitude"][:]

    nc.close()

    lon2d, lat2d = np.meshgrid( lon1d, lat1d )

    return( rain2d, lon2d, lat2d )

def dbz2rain( dbz ):
    # Suezawa
    B_ = 115.7
    beta_ = 1.681
    const_ = np.log10( B_ )
    return( np.power( 10, ( dbz*0.1 - const_ ) / beta_ ) )
    #return( np.power( np.power( 10, dbz*0.1) / 200.0, 1.0/1.6 ) ) 

def read_ga_grads_dbz( INFO, itime=datetime(2019,9,3,2,0,0), typ="g" ):
    # read guess/analsis ensemble mean  

    if typ == "g":
       fn_ = itime.strftime( 'gues_pawr_%Y%m%d-%H%M%S.grd')
    elif typ == "a":
       fn_ = itime.strftime( 'anal_pawr_%Y%m%d-%H%M%S.grd')

    fn = os.path.join( INFO["TOP"], INFO["EXP"],
                       INFO["time0"].strftime('%Y%m%d%H%M%S'),
                       "mean_grads", fn_ )

    try:
       infile = open(fn)
    except:
       print("Failed to open")
       print( fn )
       return( None, False )
       sys.exit()

    print( fn )
    gz = 60
    gx = INFO["gx"]
    gy = INFO["gy"]

    rec3d = gx*gy*gz
    rec2d = gx*gy

    # u, v, w, t, p, qv, qc, qr, qi, qs, qg
    nv3d = 11
    # PW PRCP
    nv2d = 0

    nv2d_ = 0
    nv3d_ = 0

    rec = rec3d * nv3d_ + rec2d * nv2d_

    try:
       infile.seek( rec*4 )
       tmp3d = np.fromfile( infile, dtype=np.dtype('>f4'), count=rec3d )  # big endian   
       input3d = np.reshape( tmp3d, (gz,gy,gx) )
    except:
       input3d = None

    return( input3d )

def draw_rec_4p( ax, lon_l=[], lat_l=[], lc='k', lw=1.0, transform=None ):

    ax.plot( [ lon_l[0], lon_l[0] ], [ lat_l[0], lat_l[1] ],  color=lc, lw=lw, 
             transform=transform )

    ax.plot( [ lon_l[1], lon_l[1] ], [ lat_l[0], lat_l[1] ],  color=lc, lw=lw, 
             transform=transform )
    ax.plot( [ lon_l[0], lon_l[1] ], [ lat_l[0], lat_l[0] ],  color=lc, lw=lw, 
             transform=transform )

    ax.plot( [ lon_l[0], lon_l[1] ], [ lat_l[1], lat_l[1] ],  color=lc, lw=lw, 
             transform=transform )

#    x, y = m( lon2d[-1,:], lat2d[-1,:] )
#    ax.plot( x, y, color=lc, lw=lw )
#
#    x, y = m( lon2d[:,0], lat2d[:,0] )
#    ax.plot( x, y, color=lc, lw=lw )
#
#    x, y = m( lon2d[:,-1], lat2d[:,-1] )
#    ax.plot( x, y, color=lc, lw=lw )

#############################
if __name__ == "__main__":
    dist() 
#    read_mask_full()
