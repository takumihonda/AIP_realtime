import numpy as np
from datetime import datetime, timedelta

import os
import sys

quick = True
#quick = False

USE_ARCH_DAT = True
#USE_ARCH_DAT = False

def read_ts( fn_ts, time=datetime(2019,6,10,8,10) ):
    data = np.load( fn_ts, allow_pickle=True )
    print( fn_ts )

    ts_l = data['ts'] 
    bs_l = data['bs'] 
    fa_l = data['fa'] 
    time_l = data['times'] 


    return( ts_l, bs_l, fa_l, time_l )
 

def main( INFO, stime_l=[], etime_l=[],
           theight=3000, thrs_dbz=15.0 ):

    data_path = "../../dat4figs_JAMES/Fig15"
    ofig = "Fig15.pdf"
    os.makedirs( data_path, exist_ok=True )
    fn = '{0:}/data.npz'.format( data_path, )


    itmax = int( ( etime_l[0] - stime_l[0] ).total_seconds()/30.0 + 1 )

    ts_l = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )
    bs_l = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )
    fa_l = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )

    ts_l_persist = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )

#    ts_l_new = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )
#    bs_l_new = np.zeros( ( len( stime_l), INFO["NEXP"], itmax, INFO["TMAX"] )  )

    t_l = np.arange( 0, 30*INFO["TMAX"], 30 )

    for i in range( len( stime_l ) ):
        fn = '{0:}/data{1:}.npz'.format( data_path, i )

        if not USE_ARCH_DAT:

           it = 0
           time = stime_l[i]
           while time <= etime_l[i]:
       
               for n in range( INFO["NEXP"]) :
                   odir = "ts_npz/" + INFO["EXP" + str(n+1)]
                   fn_ts = odir + "/20211202TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
                   ts_l[i,n,it,:], bs_l[i,n,it,:], fa_l[i,n,it,:], _ = read_ts( fn_ts )
       
                   fn_ts_persist = odir + "/persistent_20211202TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
                   ts_l_persist[i,n,it,:], _, _,  _ = read_ts( fn_ts_persist )

                   #fn_ts_new = odir + "/20211202TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
                   #ts_l_new[i,n,it,:], bs_l_new[i,n,it,:], _ = read_ts( fn_ts )

               it += 1
               time += timedelta( seconds=30 ) 
   
           np.savez( fn, ts_l=ts_l[i,:,:,:],  bs_l=bs_l[i,:,:,:], fa_l=fa_l[i,:,:,:] )
        else:

           ts_l[i,:,:,:] = np.load( fn )["ts_l"]
           bs_l[i,:,:,:] = np.load( fn )["bs_l"]
           fa_l[i,:,:,:] = np.load( fn )["fa_l"]

    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.cm as cm


    fig, ( ax1, ax2, ax3 ) = plt.subplots( 1, 3, figsize=( 14, 3.8 ) )
    fig.subplots_adjust( left=0.04, bottom=0.15, right=0.99, top=0.9, 
                         wspace=0.15 )

    # average
    ts_l = np.nanmean( ts_l, axis=2 )
    bs_l = np.nanmean( bs_l, axis=2 )
    fa_l = np.nanmean( fa_l, axis=2 )

    ts_l_persist = np.nanmean( ts_l_persist, axis=2 )

#    ts_l_new = np.nanmean( ts_l_new, axis=2 )
#    bs_l_new = np.nanmean( bs_l_new, axis=2 )

    lw = 2.0
    c_l = [ 'k', 'k', 'k' ]
    c_l_old = [ 'b', 'b' ]
    ax_l = [ ax1, ax2, ax3 ]
    lab_l = [ "August 24", "August 19" ]
    ls_l = [ "solid", "dashed" ]
 
    t_l = t_l / 60 #sec 

    for j in range( len(stime_l) ):
        ax1.plot( t_l, ts_l[j,0,:], lw=lw, color=c_l[j], label=lab_l[j], ls=ls_l[j] ) 
        ax2.plot( t_l, bs_l[j,0,:], lw=lw, color=c_l[j], label=lab_l[j], ls=ls_l[j] ) 
        ax3.plot( t_l, fa_l[j,0,:], lw=lw, color=c_l[j], label=lab_l[j], ls=ls_l[j] ) 

        ax1.plot( t_l, ts_l_persist[j,0,:], lw=lw, color='gray', label=lab_l[j]+"\n(persistent)", ls=ls_l[j] ) 

#        ax1.plot( t_l, ts_l_new[j,0,:], lw=lw, color=c_l_old[j], label=lab_l[j]+ " new", ls=ls_l[j], zorder=0 ) 
#        ax2.plot( t_l, bs_l_new[j,0,:], lw=lw, color=c_l_old[j], label=lab_l[j]+" new", ls=ls_l[j], zorder=0 ) 

    ymin1 = 0
    ymax1 = 1
    ymin2 = 0.0
    ymax2 = 3.0
    ymin3 = 0.0
    ymax3 = 1.0
    ymin_l = [ ymin1, ymin2, ymin3 ]
    ymax_l = [ ymax1, ymax2, ymax3 ]
    tit_l = [ "Threat score", 
              "Bias score",
              "False alarm ratio" ]
    note = "Z={:.1f} km\n{:.1f} dBZ".format(theight/1000, thrs_dbz )
    pnum_l = [ "(a)", "(b)", "(c)" ] 

    xlab = "Forecast time (min)"

    lloc_l = [ 'upper right', 
               'upper right',
               'lower right']

    for i, ax in enumerate( ax_l ):
        #ax.legend( fontsize=12, loc='lower left' )
        ax.legend( fontsize=12, loc=lloc_l[i] )
        ax.set_xlim(0, 30)
        ax.set_ylim( ymin_l[i], ymax_l[i] )

        ax.text( 0.5, 1.02, tit_l[i],
                 fontsize=15, transform=ax.transAxes,
                 ha='center',
                 va='bottom' )
 
        ax.text( 0.0, 1.02, pnum_l[i],
                 fontsize=14, transform=ax.transAxes,
                 ha='left',
                 va='bottom' )

        if i == 2:
           ax.text( 1.0, 1.01, note,
                    fontsize=10, transform=ax.transAxes,
                    ha='right',
                    va='bottom' )
  
        ax.set_xlabel( xlab, fontsize=12 )


#    ofig = "2p_scores.png"
 
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




nexp = 1

EXP1 = "D4_500m_CTRL"
LAB1 = "D4_500m_CTRL"

EXP2 = "D4_500m_H4V1"
LAB2 = "D4_500m_H4V1"

EXP3 = "D4_500m_H8V1"
LAB3 = "D4_500m_H8V1"

EXP4 = "D4_500m_H1V1"
LAB4 = "D4_500m_H1V1"

EXP5 = "D4_500m_H8V8"
LAB5 = "D4_500m_H8V8"

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

main( INFO, stime_l=stime_l, etime_l=etime_l, theight=theight, thrs_dbz=thrs_dbz )

