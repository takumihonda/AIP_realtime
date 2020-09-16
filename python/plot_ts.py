import numpy as np
from datetime import datetime, timedelta

import sys

AVE = True
AVE = False

FT0 = True
BS = True


def read_ts( fn_ts, time=datetime(2019,6,10,8,10) ):
    data = np.load( fn_ts, allow_pickle=True )

    ts_l = data['ts'] 
    bs_l = data['bs'] 
    time_l = data['times'] 

    return( ts_l, bs_l, time_l )
 

def main( INFO, stime=datetime(2019,6,10,8,10,0), etime=datetime(2019,8,24,15,0),
           theight=3000, thrs_dbz=15.0 ):

    odir = "ts_npz/" + INFO["EXP"]

    itmax = int( ( etime - stime ).total_seconds()/30.0 + 1 )


    it = 0
    time = stime
    while time <= etime:
        it += 1
#        fn_ts = odir + "/TS_" + "thrs{:.1f}dbz_".format(thrs_dbz) + \
#                itime.strftime('i%H%M%S_%Y%m%d.npz')
        fn_ts = odir + "/TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
    
        ts_l, bs_l, time_l = read_ts( fn_ts )

        if it > 1:
           ave_ts_l += ts_l
           ave_bs_l += bs_l
        elif it == 1:
           ave_ts_l = ts_l
           ave_bs_l = bs_l

           if FT0:
              t_l = np.arange( 0, 30*len(ts_l), 30 )
           else:
              t_l = np.arange( 30, 30*(len(ts_l)+1), 30 )
       
           tss_l = np.zeros( ( len(t_l), itmax ) )
           bss_l = np.zeros( ( len(t_l), itmax ) )

        tss_l[:,it-1] = ts_l
        bss_l[:,it-1] = bs_l

        time += timedelta( seconds=30 ) 

    ave_ts_l = ave_ts_l / it


    if AVE:
       data_l = ave_ts_l
    else:
       data_l = ts_l

    data_l = np.mean( tss_l, axis=1 )  
    if BS:
       data_l = np.mean( bss_l, axis=1 )  

    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.cm as cm


    fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, )

#    t_l = np.arange(30, (len(time_l)+1)*30, 30)

    if AVE:
       ax1.plot(t_l, data_l, lw=4.0, color='k', linestyle='solid')
       ax1.set_xlim(0, 1800)
    else:
       for it in range(itmax):
          if BS:
             ax1.plot( t_l, bss_l[:,it], lw=0.5, linestyle='solid', alpha=0.5,
                       color=cm.jet(it/itmax))
          else:
             ax1.plot( t_l, tss_l[:,it], lw=0.5, linestyle='solid', alpha=0.5,
                       color=cm.jet(it/itmax))

#       it = 0
#       time = stime
#       while time <= etime:
       ax1.set_xlim(0, 1800)
          
    print( "itmax:", itmax)
    ymax = 1.0
    if BS:
       ymax = 5.0

    ax1.set_ylim(0, ymax )

    dxsec = 300
    xsecs = np.arange(dxsec, 1800, dxsec)
    ax1.vlines( x=xsecs, ymin=0, ymax=ymax, ls='dashed', lw=0.5)

    ax1.set_xticks( xsecs )

    ylab = "Threat score"
    tit = "Threat score"
    if BS:
       ylab = "Bias score"
       tit = "Bias score"
       ax1.hlines( xmin=0, xmax=1800, y=1.0, ls='dashed', lw=0.5)

    ax1.set_xlabel("Forecast time (s)", fontsize=12)
    ax1.set_ylabel( ylab, fontsize=12)
#    seclocator = mdates.SecondLocator(bysecond=[20, 40]) 
#    minlocator = mdates.MinuteLocator(byminute=range(600))  # range(60) is the default
#    
##    seclocator.MAXTICKS  = 40000
##    minlocator.MAXTICKS  = 40000
#    
##    majorFmt = mdates.DateFormatter('%Y-%m-%d, %H:%M:%S')  
#    minorFmt = mdates.DateFormatter('%H:%M:%S')  
#    
#    ax1.xaxis.set_major_locator(minlocator)
#    ax1.xaxis.set_major_formatter(majorFmt)
 
    ax1.text( 0.5, 1.02, tit,
              fontsize=15, transform=ax1.transAxes,
              horizontalalignment='center',
              verticalalignment='bottom' )
   
    note = "Z={:.1f}km\n{:.1f}dBZ".format(theight/1000, thrs_dbz)
    ax1.text( 1.0, 1.01, note,
              fontsize=12, transform=ax1.transAxes,
              horizontalalignment='right',
              verticalalignment='bottom' )

    plt.show()

########

TOP = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY"

EXP = "D4_500m_TEST_DEFAULT"
EXP = "D4_500m_TEST_DEFAULT_MEAN"

INFO = { "TOP": TOP,
         "EXP": EXP,
       }

theight = 3000.0
#theight = 6000.0
thrs_dbz = 25.0
thrs_dbz = 15.0
thrs_dbz = 30.0


#time = datetime(2019,9,3,2,50,0)
#time = datetime(2019, 6, 10, 8, 20, 0)
stime = datetime( 2019, 8, 24, 15, 0, 30)
etime = datetime( 2019, 8, 24, 16, 0, 0)

main( INFO, stime=stime, etime=etime, theight=theight, thrs_dbz=thrs_dbz )

