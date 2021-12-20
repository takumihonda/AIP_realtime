import numpy as np
from datetime import datetime, timedelta

import os
import sys

quick = True
#quick = False

USE_ARCH_DAT = True
USE_ARCH_DAT = False


def read_ts( fn_ts, time=datetime(2019,6,10,8,10) ):
    data = np.load( fn_ts, allow_pickle=True )

    ts_l = data['ts'] 
    bs_l = data['bs'] 
    time_l = data['times'] 


    return( ts_l, bs_l, time_l )
 

def main( INFO, stime_l=[], etime_l=[],
           theight=3000, thrs_dbz=15.0, lab_l=[]):

    data_path = "../../dat4figs_JAMES/Fig16_20211209"
    ofig = "Fig16_20211209.pdf"
    os.makedirs( data_path, exist_ok=True )
    fn = '{0:}/data.npz'.format( data_path, )

    itmax = int( ( etime_l[0] - stime_l[0] ).total_seconds()/30.0 + 1 )

    ts_l = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )
    bs_l = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )

    t_l = np.arange( 0, 30*INFO["TMAX"], 30 )

    for i in range( len( stime_l ) ):

        fn = '{0:}/data{1:}.npz'.format( data_path, i )

        if not USE_ARCH_DAT:

           it = 0
           time = stime_l[i]
           while time <= etime_l[i]:
       
               for n in range( INFO["NEXP"]) :
                   odir = "ts_npz/" + INFO["EXP" + str(n+1)]
                   #fn_ts = odir + "/TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
                   fn_ts = odir + "/20211202TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
                   ts_l[i,n,it,:], bs_l[i,n,it,:], _ = read_ts( fn_ts )
                   print( fn_ts, ts_l[i,n,it,0] )
       
               it += 1
               time += timedelta( seconds=30 ) 

#           np.savez( fn, ts_l=ts_l[i,:,:,:],  bs_l=bs_l[i,:,:,:] )
        else:

           ts_l[i,:,:,:] = np.load( fn )["ts_l"]
           bs_l[i,:,:,:] = np.load( fn )["bs_l"]


    print( ts_l.shape )

    e1 = 2
    e2 = 1


    tmin = 0
    tmax = 61
    t_l = np.zeros( tmax )
    t_l[:] = np.nan

    p_thrs = 0.05 # 95%

    for t in range( tmin, tmax ):

        data1 = ts_l[:,e1,:,t].flatten()
        data2 = ts_l[:,e2,:,t].flatten()
    
    
        dif = data1 - data2
    
        data1 = data1[~np.isnan( dif ) ]
        data2 = data2[~np.isnan( dif ) ]
    
        dif = data1 - data2
        print( dif.shape )
    
        import scipy.stats
        result= scipy.stats.ttest_ind( dif,
                                       np.zeros( len( dif ) ),
                                       equal_var=False,
                                     )
        print( result )

        p_value = result[1]
        if p_value < p_thrs:
           t_l[t] = t
 
    print( lab_l[e1], lab_l[e2])       
    print( t_l*30/60 )
    #    sys.exit()
    
        
    
#        import matplotlib.pyplot as plt
#        fig, ax = plt.subplots( 1, 1, figsize=( 8, 4 ) )
#    
#    
#        data = dif
#    
#        xmin = -0.3
#        xmax = 0.3
#        sigma = np.std( data, ddof=1 )
#        mu = np.mean( data, )
#        h = 3.5*sigma / np.power( data.size, 1.0/3.0 )
#        print( data )
#        print( "check ", sigma, h )
#        nbin = int( ( xmax - xmin ) / h )
#        ax.hist( data, bins=nbin, range=( xmin, xmax ), density=True )
#        ax.set_xlim( xmin, xmax )
#    
#        import scipy.stats as stats
#        x_l = np.arange( xmin, xmax+h, h )
#        ax.plot( x_l, stats.norm.pdf(x_l, mu, sigma))
#    
#        ax.text( 0.5, 1.02, '{0:}-{1:}'.format( lab_l[e1], lab_l[e2] ),
#                 fontsize=12, transform=ax.transAxes,
#                 ha='center',
#                 va='bottom' )
#    
#        plt.show()
#    
    sys.exit()

    bootstrap( A=ts_l[0,0,:,0].tolist(), B=ts_l[0,1,:,0].tolist() )
    sys.exit()

def bootstrap( A=[], B=[], nmax=100 ):
    org_dif = np.mean( A ) - np.mean( B )

    import random
    import statistics

    dif_l = []

    G = A + B 

    for n in range( nmax ):
       random.shuffle( G )
       new_A = G[:len( A )] 
       new_B = G[len( A ):] 

       dif_ = statistics.mean( new_A ) - statistics.mean( new_B )
       dif_l.append( dif_ )

    print( dif_l )

    import matplotlib.pyplot as plt
    plt.hist( dif_l )
    plt.axvline( org_dif, color='k' )
    print( np.quantile( dif_l, [0.95,0.05] ) )

    plt.show()

##################

otyp = 4002 # vr
#otyp = 4004 # zero Z





########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"



tmax = 61

EXP1 = None
LAB1 = None
EXP2 = None
LAB2 = None
EXP3 = None
LAB3 = None
EXP4 = None
LAB4 = None
EXP5 = None
LAB5 = None
EXP6 = None
LAB6 = None
EXP7 = None
LAB7 = None
EXP8 = None
LAB8 = None



lab_l = [
        "H1V1",
        "H4V4 (CTRL)",
        "H8V8",
        ]

nexp = 3

EXP1 = "D4_500m_H1V1"
LAB1 = "D4_500m_H1V1"

EXP2 = "D4_500m_CTRL"
LAB2 = "D4_500m_CTRL"

EXP3 = "D4_500m_H8V8"
LAB3 = "D4_500m_H8V8"

#EXP2 = EXP1
#LAB2 = LAB1
#EXP3 = EXP1
#LAB3 = LAB1
#EXP4 = EXP1
#LAB4 = LAB1
#EXP5 = EXP1
#LAB5 = LAB1

INFO = { "TOP": TOP,
         "NEXP": nexp,
         "EXP1": EXP1,
         "EXP2": EXP2,
         "EXP3": EXP3,
         "EXP4": EXP4,
         "EXP5": EXP5,
         "EXP6": EXP6,
         "EXP7": EXP7,
         "EXP8": EXP8,
         "LAB1": LAB1,
         "LAB2": LAB2,
         "LAB3": LAB3,
         "LAB4": LAB4,
         "LAB5": LAB5,
         "LAB6": LAB6,
         "LAB7": LAB7,
         "LAB8": LAB8,
         "TMAX": tmax,
       }

theight = 3000.0
#theight = 6000.0
thrs_dbz = 15.0
thrs_dbz = 30.0


stime1 = datetime( 2019, 8, 24, 15, 0, 30)
etime1 = datetime( 2019, 8, 24, 16, 0, 0)

stime2 = datetime( 2019, 8, 19, 13, 0, 30)
etime2 = datetime( 2019, 8, 19, 14, 0, 0)

stime_l = [ stime1, stime2 ]
etime_l = [ etime1, etime2 ]


main( INFO, stime_l=stime_l, etime_l=etime_l, theight=theight, thrs_dbz=thrs_dbz, lab_l=lab_l )

