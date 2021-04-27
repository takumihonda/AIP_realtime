import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

from tools_AIP import read_obs, read_mask, read_nc_lonlat, read_fcst_grads, read_obs_grads_latlon, dbz2rain, read_nowcast_hires

quick = True
quick = False




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

    tlev = int( tlev )
    ftime = itime + timedelta( seconds=int(tlev)*30 )

    obs3d, ostat = read_obs( utime=ftime, mask=INFO["mask"] )

    fcst3d, fstat,= read_fcst_grads( INFO, itime=itime, tlev=tlev, )

    ifcst2d = vint_fcst( INFO["hgt3d"], fcst3d, theight=theight )

    olon2d = INFO["olon2d"]
    olat2d = INFO["olat2d"]

    from scipy.interpolate import griddata
    # undef values
    if ostat: 
       obs3d[ obs3d < 10.0 ] = 5.0
       ok = np.argmin( np.abs( INFO["obsz"] - theight ) )

       iobs2d = obs3d[ok,:,:]
       obs_  =  iobs2d[ ( olon2d >= lons ) & ( olon2d <= lone ) &
                        ( olat2d >= lats ) & ( olat2d <= late ) ]
    else:
       obs_ = np.nan


    ifcst2d = griddata( ( INFO["lon2d"].ravel(), INFO["lat2d"].ravel()), ifcst2d.ravel(),
                       (INFO["olon2d"], INFO["olat2d"]),
                       method='cubic',
                      )



    fcst_ = ifcst2d[ ( olon2d >= lons ) & ( olon2d <= lone ) &
                     ( olat2d >= lats ) & ( olat2d <= late ) ]


    try:
       jma2d, jlon2d, jlat2d = read_nowcast_hires( stime=itime, ft=timedelta(seconds=tlev*30) )
       ijma2d = griddata( ( jlon2d.ravel(), jlat2d.ravel()), jma2d.ravel(),
                          (INFO["olon2d"], INFO["olat2d"]),
                          method='cubic',
                         )
       jma_ = ijma2d[ ( olon2d >= lons ) & ( olon2d <= lone ) &
                      ( olat2d >= lats ) & ( olat2d <= late ) ]
       jmax = np.nanmax( jma_ )
    except:
       print( 'Failed to read JMA nowcast' )
       jmax = np.nan

    return( np.nanmax( dbz2rain( fcst_ ) ), np.nanmax( dbz2rain( obs_ ) ), jmax )




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
#etime = datetime( 2019, 8, 24, 15, 40, 0 )
etime = datetime( 2019, 8, 24, 15, 50, 0 )
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

fn_obs = "obs_rainmax_z{0:.1f}_lon{1:.2f}-{2:.2f}_lat{3:.2f}-{4:.2f}.npz".format( theight, lons, lone, lats, late )

obs_l = []
otime_l = []

OBS = True
FCST = False
NOWCAST = False

time = stime
while (time <= etime):
  print( "Initial time:", time )


  fn_fcst = "fcst_rainmax_z{0:.1f}_lon{1:.2f}-{2:.2f}_lat{3:.2f}-{4:.2f}_i{5:}.npz".format( theight, lons, lone, lats, late, time.strftime('%Y%m%d%H%M%S') )
  fn_nowcast = "nowcast_rainmax_z{0:.1f}_lon{1:.2f}-{2:.2f}_lat{3:.2f}-{4:.2f}_i{5:}.npz".format( theight, lons, lone, lats, late, time.strftime('%Y%m%d%H%M%S') )

  ftime_l = []
  fcst_l = []
  jtime_l = []
  jma_l = []
  for i, tlev in enumerate( tlevs ):
      fcst_, obs_, jma_ = main( INFO, itime=time, tlev=tlev, theight=theight, 
                             lons=lons, lone=lone, lats=lats, late=late )
      print( tlev, fcst_, obs_, jma_ )

      if i == 0:
         obs_l.append( obs_ )
         otime_l.append( time )

      if FCST:
         fcst_l.append( fcst_ )
         ftime_l.append( time + timedelta(seconds=int( tlev*30 )) )

      if NOWCAST:
         if not np.isnan( jma_ ): 
            jma_l.append( jma_ )
            jtime_l.append( time + timedelta(seconds=int( tlev*30 )) )


      if not FCST and not NOWCAST:
         break

#  print( obs_l )
#  print( otime_l )
#  print( fcst_l )
#  print( ftime_l )
#  print( jma_l )
#  print( jtime_l )
#  sys.exit()

  # fcst
  if FCST:
     np.savez( os.path.join( odir, fn_fcst ), rmax=np.array( fcst_l[:] ), times=np.array( ftime_l[:] ) ) 

  # nowcast
  if NOWCAST:
     if len( jma_l ) > 0: 
        np.savez( os.path.join( odir, fn_nowcast ), rmax=np.array( jma_l[:] ), times=np.array( jtime_l[:] ) ) 

  time += timedelta(seconds=30)

# obs
if OBS:
   print( obs_l )
   np.savez( os.path.join( odir, fn_obs ), rmax=np.array( obs_l[:] ), times=np.array( otime_l[:] ) ) 

  

