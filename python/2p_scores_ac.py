import numpy as np
from datetime import datetime, timedelta

import os
import sys

quick = True
#quick = False

#USE_ARCH_DAT = True
USE_ARCH_DAT = False


def read_ts( fn_ts, time=datetime(2019,6,10,8,10) ):
    data = np.load( fn_ts, allow_pickle=True )

    ts_l = data['ts'] 
    bs_l = data['bs'] 
    time_l = data['times'] 
    fss_l = data['fss'] 
    ng_l = data['ng_l'] 


    return( ts_l, bs_l, time_l, fss_l, ng_l )
 

def main( INFO, stime_l=[], etime_l=[],
           theight=3000, thrs_dbz=15.0, lab_l=[], AC=False,
           exp_l=[], nglev=1 ):

#    data_path = "../../dat4figs_JAMES/Fig16_20211209"
#    ofig = "Fig16_20211209.pdf"
#    os.makedirs( data_path, exist_ok=True )
#    fn = '{0:}/data.npz'.format( data_path, )

    itmax = int( ( etime_l[0] - stime_l[0] ).total_seconds()/30.0 + 1 )

    ts_l = np.zeros( ( len( stime_l), len( exp_l ), itmax, INFO["TMAX"] )  )
    bs_l = np.zeros( ( len( stime_l), len( exp_l ), itmax, INFO["TMAX"] )  )

    t_l = np.arange( 0, 30*INFO["TMAX"], 30 )

    ng_l = [ 0, 1, 2, 3, 4, 6, ]
    fss_l = np.zeros( ( len( stime_l), len( exp_l ), itmax, INFO["TMAX"], len( ng_l ) ) )
    fss_l[:] = np.nan

    for i in range( len( stime_l ) ):

#        fn = '{0:}/data{1:}.npz'.format( data_path, i )

        if not USE_ARCH_DAT:

           it = 0
           time = stime_l[i]
           while time <= etime_l[i]:
       
               for n in range( len( exp_l ) ) :
#                   odir = "ts_npz/" + INFO["EXP" + str(n+1)]

                   if AC:
                      odir = "ts_npz/{0:}_acObs".format( exp_l[n] )
                   else:
                      odir = "ts_npz/{0:}_noacObs".format( exp_l[n] )

                   fn_ts = odir + "/20220216_ac_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
                   ts_l[i,n,it,:], bs_l[i,n,it,:], _,  fss_l[i,n,it,:,:], ng_l = read_ts( fn_ts )
                   print( fn_ts, ts_l[i,n,it,0] )
       
               it += 1
               time += timedelta( seconds=30 ) 

#           np.savez( fn, ts_l=ts_l[i,:,:,:],  bs_l=bs_l[i,:,:,:] )
        else:

           ts_l[i,:,:,:] = np.load( fn )["ts_l"]
           bs_l[i,:,:,:] = np.load( fn )["bs_l"]


    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.cm as cm


#    fig, ( ( ax1, ax2 ), ( ax3, ax4 ) ) = plt.subplots( 2, 2, figsize=( 11, 8 ) )
#    fig.subplots_adjust( left=0.05, bottom=0.08, right=0.98, top=0.95, 
#                         wspace=0.1, hspace=0.3 )


    #fig, ( ax1, ax2, ax3 ) = plt.subplots( 1, 3, figsize=( 11, 4 ) )
    fig, ( ax1, ax2, ) = plt.subplots( 1, 2, figsize=( 11, 4 ) )
    fig.subplots_adjust( left=0.05, bottom=0.15, right=0.98, top=0.9, 
                         wspace=0.1 )

    print( ts_l.shape )
    # count sample
    for i in range( len( exp_l ) ):
        tmp = ts_l[:,i,:,0]
        nnan = len( tmp[ np.isnan( tmp )] ) 
        print( i, "Nan: ", nnan )
#    sys.exit()
    # average
    ts_l = np.nanmean( ts_l, axis=2 )
    bs_l = np.nanmean( bs_l, axis=2 )
    nglev = 1
    fss_l = np.nanmean( fss_l[:,:,:,:,nglev], axis=2 )

    lw = 2.0
    c_l = [ 'k', 'r' ]
    ax_l = [ ax1, ax2, ] #ax3 ]
#    lab_l = [ "August 24", "August 19" ]
    ls_l = [ "solid", "solid", "solid" ]
#    ls_l = [ "solid", "dashed" ]
 
    t_l = t_l / 60 #sec 

#    for j in range( len(stime_l) ):
#        for n in range( INFO["NEXP"]) :
#            ax1.plot( t_l, ts_l[j,n,:], lw=lw, color=c_l[n], label=lab_l[n], ls=ls_l[j] ) 
#            ax2.plot( t_l, bs_l[j,n,:], lw=lw, color=c_l[n], label=lab_l[n], ls=ls_l[j] ) 

    ymin1 = 0
    ymax1 = 0.8
    ymin2 = 0.0
    ymax2 = 3.0
    ymin3 = 0.0
    ymax3 = 1.0

    ymin_l = [ ymin1, ymin2, ymin3 ]
    ymax_l = [ ymax1, ymax2, ymax3 ]
    tit_l = [ "Threat score", "Bias score", "FSS" ]
    note = "Z={:.1f} km\n{:.1f} dBZ".format(theight/1000, thrs_dbz )
    pnum_l = [ "(a)", "(b)", "(c)", "(d)" ] 

    xlab = "Forecast time (min)"
    print( "chk", ts_l.shape )

    for i in range( len(ax_l) ):
        if i == 0:
           ax = ax1
           ii = 0 # TS
        elif i == 1:
           ax = ax2
           ii = 1 # BS
        elif i == 2:
           ax = ax3
           ii = 2 # FSS
        elif i == 3:
           ax = ax4
           ii = 1 # BS

        if i<= 2:
           st = 0
        else:
           st = 1

        if ii == 0:
           data = ts_l[st,:,:]
           lloc = 'upper right'
        elif ii == 1:
           data = bs_l[st,:,:]
           lloc = 'lower left'
        elif ii == 2:
           data = fss_l[st,:,:]
           lloc = 'lower left'

        for n in range( INFO["NEXP"]) :
            ax.plot( t_l, data[n,:], lw=lw, color=c_l[n], 
                     label=lab_l[n], ls=ls_l[0] ) 

        ax.legend( fontsize=12, loc=lloc )
        ax.grid( ls='dashed', lw=1.0 )
        ax.set_xlim(0, 30)
        ax.set_ylim( ymin_l[ii], ymax_l[ii] )

        ax.text( 0.5, 1.02, tit_l[ii],
                 fontsize=15, transform=ax.transAxes,
                 ha='center',
                 va='bottom' )
 
        ax.text( 0.0, 1.02, pnum_l[i],
                 fontsize=14, transform=ax.transAxes,
                 ha='left',
                 va='bottom' )

        ax.text( 1.0, 1.01, note,
                 fontsize=10, transform=ax.transAxes,
                 ha='right',
                 va='bottom' )
  
        ax.set_xlabel( xlab, fontsize=12 )


    ofig = "TESTTT2p_scores_thinning.png"
 
    print( ofig )
    if quick:
       plt.show()
    else:
       odir = "pdf/"
       os.makedirs( odir, exist_ok=True)
       plt.savefig( os.path.join(odir, ofig), 
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')




##################

otyp = 4002 # vr
#otyp = 4004 # zero Z





########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"



tmax = 61


nexp = 2

INFO = { "TOP": TOP,
         "NEXP": nexp,
         "TMAX": tmax,
       }

theight = 3000.0
#theight = 6000.0
thrs_dbz = 15.0
thrs_dbz = 30.0


stime1 = datetime( 2019, 8, 24, 15, 0, 30)
etime1 = datetime( 2019, 8, 24, 16, 0, 0)
#etime1 = stime1

stime2 = datetime( 2019, 8, 19, 13, 0, 30)
etime2 = datetime( 2019, 8, 19, 14, 0, 0)

stime_l = [ stime1, stime2 ]
etime_l = [ etime1, etime2 ]


stime_l = [ stime1, stime2 ]
AC = True
exp_l = [ 'd4', 'd4_500m_ac', ]
#exp_l = [ 'd4', 'd4_500m_ac_vr', ]

lab_l = [
        "CTRL",
        "TEST",
        ]

nglev = 2

main( INFO, stime_l=stime_l, etime_l=etime_l, theight=theight, thrs_dbz=thrs_dbz, lab_l=lab_l, AC=AC, exp_l=exp_l, nglev=nglev )

