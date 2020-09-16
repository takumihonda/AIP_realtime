import numpy as np
from datetime import datetime, timedelta

import sys

AVE = True
#AVE = False

FT0 = True

BS = True
#BS = False

def read_ts( fn_ts, time=datetime(2019,6,10,8,10) ):
    data = np.load( fn_ts, allow_pickle=True )

    ts_l = data['ts'] 
    bs_l = data['bs'] 
    time_l = data['times'] 

    return( ts_l, bs_l, time_l )
 

def main( INFO, stime=datetime(2019,6,10,8,10,0), etime=datetime(2019,8,24,15,0),
           theight=3000, thrs_dbz=15.0 ):


    itmax = int( ( etime - stime ).total_seconds()/30.0 + 1 )


    it = 0
    time = stime
    while time <= etime:
        it += 1
#        fn_ts = odir + "/TS_" + "thrs{:.1f}dbz_".format(thrs_dbz) + \
#                itime.strftime('i%H%M%S_%Y%m%d.npz')
        odir = "ts_npz/" + INFO["EXP1"]
        fn_ts = odir + "/TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
    
        ts_l, bs_l, time_l = read_ts( fn_ts )

        odir = "ts_npz/" + INFO["EXP2"]
        fn_ts2 = odir + "/TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
        ts2_l, bs2_l, time2_l = read_ts( fn_ts2 )

        if it > 1:
           ave_ts_l += ts_l
           ave_ts2_l += ts2_l
        elif it == 1:
           ave_ts_l = ts_l
           ave_ts2_l = ts2_l

           if FT0:
              #t_l = np.arange( 0, 30*len(ts_l), 30 )
              t_l = np.arange( 0, 600*len(ts_l), 600 )
           else:
              t_l = np.arange( 30, 30*(len(ts_l)+1), 30 )
       
           tss_l = np.zeros( ( len(t_l), itmax ) )
           tss2_l = np.zeros( ( len(t_l), itmax ) )
           bss_l = np.zeros( ( len(t_l), itmax ) )
           bss2_l = np.zeros( ( len(t_l), itmax ) )

        tss_l[:,it-1] = ts_l
        tss2_l[:,it-1] = ts2_l
        bss_l[:,it-1] = bs_l
        bss2_l[:,it-1] = bs2_l

        time += timedelta( seconds=30 ) 

    ave_ts_l = ave_ts_l / it
    ave_ts2_l = ave_ts2_l / it


    if AVE:
       data_l = ave_ts_l
    else:
       data_l = ts_l


    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.cm as cm


    fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, )

    print( EXP1 )
    print( tss_l )

    print( EXP2 )
    print( tss2_l )
    print( "" )

    print( EXP1 )
    print( bss_l )

    print( EXP2 )
    print( bss2_l )

#    t_l = np.arange(30, (len(time_l)+1)*30, 30)

    if AVE:
       v1 = np.mean( tss_l, axis=1 )
       v2 = np.mean( tss2_l, axis=1 )
       if BS:
          v1 = np.mean( bss_l, axis=1 )
          v2 = np.mean( bss2_l, axis=1 )

       ax1.plot(t_l, v1, lw=1.0, color='k', linestyle='solid', label=INFO["LAB1"])
       ax1.plot(t_l, v2, lw=1.0, color='r', linestyle='solid', label=INFO["LAB2"])
       ax1.set_xlim(0, 1800)
       plt.legend()
    else:
       for it in range(itmax):
          ax1.plot( t_l, tss_l[:,it], lw=0.5, linestyle='solid', alpha=0.5,
                    color=cm.jet(it/itmax))
#       it = 0
#       time = stime
#       while time <= etime:
       ax1.set_xlim(0, 1800)
          
    print( "itmax:", itmax)

    ylab = "Threat score"
    tit = "Threat score"
    ymax = 1.0
    if BS:
       ymax = 2.0
       ylab = "Bias score"
       tit = "Bias score"
       ax1.hlines( xmin=0, xmax=1800, y=1.0, ls='dashed', lw=0.5)

    ax1.set_ylim(0, ymax )

    dxsec = 300
    xsecs = np.arange(dxsec, 1800, dxsec)
    ax1.vlines( x=xsecs, ymin=0, ymax=ymax, ls='dashed', lw=0.5)

    ax1.set_xticks( xsecs )

    ax1.set_xlabel("Forecast time (s)", fontsize=12 )
    ax1.set_ylabel( ylab , fontsize=12 )
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
 
    tit = tit
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

#EXP2 = "D4_500m_TEST_DEFAULT_MEAN"
#EXP = "D4_500m_TEST_DEFAULT_MEAN"

EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_NOMAX400"

LAB1 = "NOBS050"
LAB2 = "NOBS400"

EXP1 = "D4_500m_TEST_DEFAULT_0515_MEAN"
LAB1 = "Default"

EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_GRID2"
#EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG02VG02"
#EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG02VG01"

EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG02VG02_M2"
EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG04_M2"
EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG06VG06_M2"
#EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG08VG08_M2"


EXP1 = "D4_500m_TEST_DEFAULT_0604_MEAN"
LAB1 = "Default"

EXP2 = "D4_500m_TEST_DEFAULT_0708_NOBS100_NEAR_HV2"
EXP2 = "D4_500m_TEST_DEFAULT_0708_NOBS100"

EXP2 = "D4_500m_TEST_DEFAULT_0708_NOBS100_NEAR_HV4"

EXP2 = "D4_500m_M4_NEAR2_HT4"
LAB2 = "New"

EXP1 = "TEST_DEFAULT"
EXP2 = "TEST_DEFAULT_MEM01"

INFO = { "TOP": TOP,
         "EXP1": EXP1,
         "EXP2": EXP2,
         "LAB1": LAB1,
         "LAB2": LAB2,
       }

theight = 3000.0
#theight = 6000.0
thrs_dbz = 15.0
#thrs_dbz = 30.0


#time = datetime(2019,9,3,2,50,0)
#time = datetime(2019, 6, 10, 8, 20, 0)
etime = datetime( 2019, 8, 24, 15, 30, 0)

stime = etime

#stime = datetime( 2019, 8, 24, 15, 0, 30)
#etime = datetime( 2019, 8, 24, 15, 30, 0)

main( INFO, stime=stime, etime=etime, theight=theight, thrs_dbz=thrs_dbz )

