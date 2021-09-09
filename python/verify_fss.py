import numpy as np
import sys
import os
from netCDF4 import Dataset

from datetime import datetime, timedelta

from tools_AIP import read_obs, read_mask, read_nc_lonlat, read_fcst_grads, read_obs_grads_latlon, calc_fss

PLOT = True
quick = True
PLOT = False
#quick = False

REGION = True
REGION = False

def get_ts_bs( fcst2d, obs2d, thrs=10.0 ):

    fo = np.where( (fcst2d >= thrs) & (obs2d >= thrs), 1, 0 ) 
    fx = np.where( (fcst2d >= thrs) & (obs2d <  thrs), 1, 0 ) 
    xo = np.where( (fcst2d <  thrs) & (obs2d >= thrs), 1, 0 ) 

    return( np.sum(fo) / ( np.sum(fo)+np.sum(fx)+np.sum(xo) ), 
            ( np.sum(fo) + np.sum(fx) ) / ( np.sum(fo) + np.sum(xo) ) )

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

def main( INFO, itime=datetime(2019,9,3,2,0), tlev=0, theight=3000, dbz_thrs_l=[], 
          lons=0.0, lone=0.0, lats=0.0, late=0.0, ng_l=[] ):

    tlev = int( tlev )
    ftime = itime + timedelta( seconds=int(tlev)*30 )

    obs3d, ostat = read_obs( utime=ftime, mask=INFO["mask"] )

    if not ostat:
       return( None, None, ostat)

    fcst3d, fstat,= read_fcst_grads( INFO, itime=itime, tlev=tlev, )

    if not fstat:
       print( "No fcst" )
       return( None, None, fstat)

    ifcst2d = vint_fcst( INFO["hgt3d"], fcst3d, theight=theight )


    from scipy.interpolate import griddata
    # undef values
    obs3d[ obs3d < 10.0 ] = 5.0
    ok = np.argmin( np.abs( INFO["obsz"] - theight ) )


    iobs2d = obs3d[ok,:,:]
    ifcst2d = griddata( ( INFO["lon2d"].ravel(), INFO["lat2d"].ravel()), ifcst2d.ravel(),
                       (INFO["olon2d"], INFO["olat2d"]),
                       method='cubic',
                      )


    olon2d = INFO["olon2d"]
    olat2d = INFO["olat2d"]

    dbz_thrs = dbz_thrs_l[0]

    fss_l = np.zeros( ( len( ng_l), len( dbz_thrs_l) ) )
    fss_l[:,:] = np.nan

    for i, ng in enumerate( ng_l ):
       for j, dbz in enumerate( dbz_thrs_l ):
          fss_l[i,j] = calc_fss( ng=ng, thrs=dbz, fcst2d=ifcst2d, obs2d=iobs2d )

    return( fss_l )

###########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"




EXP = "D4_500m_H4V1"
EXP = "D4_500m_H8V8"

EXP = "D4_500m_CTRL_NOCLRZ"
EXP = "D4_500m_CTRL"
#EXP = "D4_500m_CTRL_MELT"
#EXP = "D4_500m_CTRL_NOVR"

theight = 3000.0
#theight = 6000.0


#etime = stime
stime = datetime( 2019, 8, 24, 15, 0, 30 )
etime = datetime( 2019, 8, 24, 16, 0, 0 )
stime = etime

#stime = datetime( 2019, 8, 24, 15, 20, 0 )
#etime = datetime( 2019, 8, 24, 15, 40, 0 )
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



dbz_thrs_l = [ 15.0, 30.0, ]

tskip = 1
tskip = 10 # debug

itmax = int( ( etime - stime ).total_seconds() / 30 ) + 1

tlevs = np.arange( tmin, tmax, tskip, dtype=np.int32 )
tlev_max = np.shape( tlevs )[0]

ng_l = [ 0, 1, 2, 3, 4, 6, ]

fss_l = np.zeros( ( tlev_max, len( ng_l), len( dbz_thrs_l ) ) )

fss_l[:] = np.nan

odir = "ts_npz/" + INFO["EXP"]
os.makedirs( odir, exist_ok=True)




lats = 36.05
late = 36.1

lons = 139.7
lone = 139.8

time = stime
while (time <= etime):
  print( "Initial time:", time )


  ftime_l = []
  for i, tlev in enumerate( tlevs ):
      fss_l_ = main( INFO, itime=time, tlev=tlev, theight=theight, 
                    dbz_thrs_l=dbz_thrs_l, lons=lons, lone=lone, lats=lats, late=late,
                    ng_l=ng_l,
                   )


#      if stat:
      fss_l[i,:,:] = fss_l_
      ftime_l.append( time + timedelta(seconds=int( tlev*30 )) )
  

  print( ng_l )
  print( fss_l[:,:,0] )
  print( fss_l[:,:,1] )

  sys.exit()
  for i, dbz in enumerate( dbz_thrs_l ):
      if REGION:
         fn_ts = "TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}_lon{3:.2f}-{4:.2f}_lat{5:.2f}-{6:.2f}.npz".format( dbz, theight, time.strftime('%H%M%S_%Y%m%d'), lons, lone, lats, late )
      else:
         fn_ts = "TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
      np.savez( os.path.join(odir,fn_ts), ts=np.array(ts_l[:,i]), bs=np.array(bs_l[:,i]), 
             times=ftime_l )
  
  time += timedelta(seconds=30)

