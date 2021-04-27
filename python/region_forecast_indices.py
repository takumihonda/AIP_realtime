import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

from tools_AIP import read_obs, read_mask, read_nc_lonlat, read_obs_grads_latlon, read_fcst_grads_all, read_fcst_qh_grads_all

quick = True
#quick = False




def vint_fcst( hgt3d, fcst3d, theight=3000.0 ):

    # get ratio for averaging
    tmp3d = np.where(hgt3d - theight >= 0.0, 0.0, hgt3d )
    hgt_below = np.where(tmp3d == np.max(tmp3d, axis=0), tmp3d, 0.0)

    tmp3d = np.where(hgt3d - theight < 0.0, 999999.0, hgt3d)
    hgt_above = np.where(tmp3d == np.min(tmp3d, axis=0), tmp3d, 0.0)

    rat2d = ( np.max(hgt_above, axis=0) - theight) / \
            ( np.max(hgt_above, axis=0) - np.max(hgt_below, axis=0) )

    fcst_below = np.max( np.where( hgt_below > 0.0, fcst3d, 0.0 ), axis=0 )
    fcst_above = np.max( np.where( hgt_above > 0.0, fcst3d, 0.0 ), axis=0 )

    return( fcst_below*rat2d + fcst_above*(np.ones(rat2d.shape)-rat2d) )

#############

def main( INFO, itime=datetime(2019,9,3,2,0), tlev=0, theight=3000, 
          lons=0.0, lone=0.0, lats=0.0, late=0.0 ):

    lon3d = INFO["lon3d"]
    lat3d = INFO["lat3d"]


    tlev = int( tlev )
    ftime = itime + timedelta( seconds=int(tlev)*30 )

    Rd = 287.04 # Markowski and Richardson 2010

    p3d  = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar='p' ) 
    tk3d = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar='t' ) 
    qv3d = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar='qv' ) 
    qv3d = qv3d / ( 1.0 - qv3d ) # specific humidity => mixing ratio
    tv3d = ( 1.0 + 0.61*qv3d ) * tk3d 
    rho3d = p3d / tv3d / Rd
    qv3d = qv3d * rho3d 

    w3d  = read_fcst_grads_all( INFO, itime=itime, tlev=tlev , FT0=True, nvar='w' ) 
    qh3d = read_fcst_qh_grads_all( INFO, itime=itime, tlev=tlev , FT0=True ) * rho3d # (kg/kg) => (kg/m3)
    nvar_l = [ w3d, ]
    max_l = []
    for nvar3d in nvar_l:
       nvar3d_  =  nvar3d[ ( lon3d >= lons ) & ( lon3d <= lone ) &
                           ( lat3d >= lats ) & ( lat3d <= late ) ]
       max_l.append( np.nanmax( nvar3d_ ) )

    # Area-averaged total (ice + liquid) water path & area-averaged PW 

    qv3d_ = np.where( ( lon3d >= lons ) & ( lon3d <= lone ) &
                      ( lat3d >= lats ) & ( lat3d <= late ),
                      qv3d, 
                      np.nan )

    qh3d_ = np.where( ( lon3d >= lons ) & ( lon3d <= lone ) &
                      ( lat3d >= lats ) & ( lat3d <= late ),
                      qh3d, 
                      np.nan )

    qv1d_ = np.nanmean( qv3d_, axis=(1,2) )
    qh1d_ = np.nanmean( qh3d_, axis=(1,2) )
    pw = qv1d_[0] * INFO["cz"][0]
    water = qh1d_[0] * INFO["cz"][0]
    for k in range( len( INFO["cz"] ) - 1  ):
        pw += ( qv1d_[k+1] + qv1d_[k] ) * 0.5 * ( INFO["cz"][k+1] - INFO["cz"][k] )
        water += ( qh1d_[k+1] + qh1d_[k] ) * 0.5 * ( INFO["cz"][k+1] - INFO["cz"][k] )

    max_l.append( water )
    max_l.append( pw )

    return( max_l )




###########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"




EXP = "D4_500m_H4V1"
EXP = "D4_500m_H8V8"

EXP = "D4_500m_CTRL_NOCLRZ"
EXP = "D4_500m_CTRL"
#EXP = "D4_500m_CTRL_MELT"
#EXP = "D4_500m_CTRL_NOVR"

theight = 2000.0


#etime = stime
stime = datetime( 2019, 8, 24, 15, 0, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )

stime = datetime( 2019, 8, 24, 15, 20, 0 )
etime = datetime( 2019, 8, 24, 15, 40, 0 )
#etime = datetime( 2019, 8, 24, 15, 50, 0 )
#etime = stime

#stime = datetime( 2019, 8, 19, 13, 0, 30 )
#etime = datetime( 2019, 8, 19, 14, 0, 0 )

# data should be stored in EXP/[time0]/dafcst
time0 = datetime( 2019, 8, 24, 15, 0, 0 )
#time0 = datetime( 2019, 8, 19, 13, 0, 0 )

tmin = 0
tmax = 61 # max time dimension does not include FT=0

FCST_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/20201117/{0:}/dafcst".format( EXP )
#OBS_DIR = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/{0:}/{1:}/pawr_grads".format( EXP )

fcst_zmax = 43

obsz, olon2d, olat2d = read_obs_grads_latlon()
lon2d, lat2d, hgt3d, cz, ohgt3d = read_nc_lonlat( fcst_zmax=fcst_zmax, obsz=obsz, NEW=True )
mask = read_mask()


INFO = {
         "mask": mask,
         "EXP": EXP, 
#         "OBS_DIR": OBS_DIR,
         "FCST_DIR": FCST_DIR,
         "obsz": obsz,
         "olon2d": olon2d,
         "olat2d": olat2d,
         "lon2d": lon2d,
         "lat2d": lat2d,
         "hgt3d": hgt3d,
         "ohgt3d": ohgt3d,
         "cz": cz,
         "gz": fcst_zmax,
         "gy": lon2d.shape[0],
         "gx": lon2d.shape[1],
       }


w3d_ = read_fcst_grads_all( INFO, itime=stime, tlev=0, FT0=True, nvar='w' ) 
lon3d = np.resize( lon2d, w3d_.shape )
lat3d = np.resize( lat2d, w3d_.shape )
INFO["lon3d"] = lon3d
INFO["lat3d"] = lat3d


tskip = 1
#tskip = 2
#tskip = 10
#tskip = 60

itmax = int( ( etime - stime ).total_seconds() / 30 ) + 1

tlevs = np.arange( tmin, tmax, tskip, dtype=np.int32 )
tlev_max = np.shape( tlevs )[0]



odir = "rainmax_npz/" + INFO["EXP"]
os.makedirs( odir, exist_ok=True )




lats = 36.05
late = 36.1
#late = 36.15

lons = 139.7
lone = 139.8



FCST = True

time = stime
while (time <= etime):
  print( "Initial time:", time )


  fn_fcst = "fcst_max_lon{0:.2f}-{1:.2f}_lat{2:.2f}-{3:.2f}_i{4:}.npz".format( lons, lone, lats, late, time.strftime('%Y%m%d%H%M%S') )

  ftime_l = []
  w_l = []
  qh_l = []
  pw_l = []
  for i, tlev in enumerate( tlevs ):
      w_, qh_, pw_ = main( INFO, itime=time, tlev=tlev, theight=theight, 
                            lons=lons, lone=lone, lats=lats, late=late )
      print( "###", tlev, w_, qh_, pw_ )


      if FCST:
         w_l.append( w_ )
         qh_l.append( qh_ )
         pw_l.append( pw_ )
         ftime_l.append( time + timedelta(seconds=int( tlev*30 )) )


  # fcst
  if FCST:
     np.savez( os.path.join( odir, fn_fcst ), wmax=np.array( w_l[:] ), 
#                qhmax=np.array( qh_l[:] ), 
                awater=np.array( qh_l[:] ), 
                apw=np.array( pw_l[:] ), 
                times=np.array( ftime_l[:] ) ) 

  time += timedelta(seconds=30)


  

