import numpy as np
from datetime import datetime, timedelta

import sys

AVE = True
#AVE = False


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

    ts_l = np.zeros( ( INFO["NEXP"], itmax, INFO["TMAX"] )  )
    bs_l = np.zeros( ( INFO["NEXP"], itmax, INFO["TMAX"] )  )

    t_l = np.arange( 0, 30*INFO["TMAX"], 30 )

    it = 0
    time = stime
    while time <= etime:

        for n in range( INFO["NEXP"]) :
            odir = "ts_npz/" + INFO["EXP" + str(n+1)]
            fn_ts = odir + "/TS_thrs{0:.1f}dbz_z{1:.1f}_i{2:}.npz".format( thrs_dbz, theight, time.strftime('%H%M%S_%Y%m%d') )
            ts_l[n,it,:], bs_l[n,it,:], _ = read_ts( fn_ts )

        it += 1
        time += timedelta( seconds=30 ) 


    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.cm as cm


    fig, (ax1) = plt.subplots(1, 1, figsize=(7,5))
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, )

#    t_l = np.arange(30, (len(time_l)+1)*30, 30)


    dat = ts_l
    ylab = "Threat score"
    tit = "Threat score"
    ymax = 1.0
    #ymax = 0.79
    if BS:
       dat = bs_l
       ymax = 5.0
       ylab = "Bias score"
       tit = "Bias score"
       ax1.hlines( xmin=0, xmax=1800, y=1.0, ls='dashed', lw=0.5)

    if AVE:
       dat = np.mean( dat, axis=1 )
       c_l = [ 'k', 'r', 'b', 'g', 'y', 'cyan', 'gray', 'p' ]
       lw_l = [ 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]
       alp = 1.0
       alp_l = [ 1.0, alp, alp, alp, alp, alp, alp, alp  ]

    for n in range( INFO["NEXP"]) :

       ax1.plot( t_l, dat[n,:], lw=lw_l[n], color=c_l[n], 
                 linestyle='solid', label=INFO["LAB"+str(n+1)],
                 alpha=alp_l[n] )

    ax1.set_xlim(0, 1800)

    plt.legend( fontsize=12 )



    ax1.set_ylim(0, ymax )

    dxsec = 300
    xsecs = np.arange(dxsec, 1800, dxsec)
    ax1.vlines( x=xsecs, ymin=0, ymax=ymax, ls='dashed', lw=0.5)

    ax1.set_xticks( xsecs )

    ax1.set_xlabel("Forecast time (s)", fontsize=12 )
    ax1.set_ylabel( ylab , fontsize=12 )

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
    sys.exit()





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
       ymax = 5.0
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


EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG02VG02_M2"
LAB2 = "H2V2"

EXP3 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG04_M2"
LAB3 = "H4V4"
EXP4 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG06VG06_M2"
LAB4 = "H6V6"
EXP5 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG08VG08_M2"
LAB5 = "H8V8"

nexp = 5
tmax = 61


EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG01_M2"
LAB2 = "H4V1"

EXP3 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG04_M2"
LAB3 = "H4V4"

EXP4 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG06_M2"
LAB4 = "H4V6"

EXP5 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG08_M2"
LAB5 = "H4V8"




nexp = 4

EXP2 = "D4_500m_TEST_DEFAULT_0515_MEAN_THIN_HG04VG04_M2"
LAB2 = "H4V4"

EXP3 = "D4_500m_TEST_DEFAULT_0604_MEAN"
LAB3 = "NEW (H4V4+small HLOC)"

EXP4 = "D4_500m_TEST_DEFAULT_0604_MEAN_PQV"
LAB4 = "NEW_PQV"

EXP4 = "D4_500m_TEST_DEFAULT_0604_MEAN_PQV_ALL"
EXP4 = "D4_500m_TEST_DEFAULT_0604_MEAN_PQV_ALL_mem_ref"
EXP4 = "D4_500m_TEST_DEFAULT_0604_MEAN_PQV_ALL_QC5"

EXP4 = "D4_500m_TEST_DEFAULT_0604_MEAN_HLOC4km"
LAB4 = "NEW (H4V4+4km HLOC)"

nexp = 2
EXP1 = "D4_500m_TEST_DEFAULT_0604_MEAN"
LAB1 = "DEFAULT"

EXP2 = "D4_500m_TEST_DEFAULT_0604_MEAN_HLOC4km_0630"
LAB2 = "NEW (H4V4+4km HLOC)"

EXP2 = "D4_500m_TEST_DEFAULT_0604_MEAN_VRMORE"
LAB2 = "+VR"

EXP2 = "D4_500m_TEST_DEFAULT_0604_MEAN_THIN4_M2"
LAB2 = "New thinning"

EXP1 = "D4_500m_TEST_DEFAULT_0708_NOBS100_NEAR2"
LAB1 = "NOBS100_NEAR2"

EXP2 = "D4_500m_TEST_DEFAULT_0708_NOBS100"
LAB2 = "NOBS100"

EXP1 = "D4_500m_TEST_DEFAULT_0604_MEAN"
LAB1 = "DEFAULT"

EXP6 = None
LAB6 = None
EXP7 = None
LAB7 = None
EXP8 = None
LAB8 = None


nexp = 3
#

EXP2 = "D4_500m_0730_M4_NEAR4_HT4_VT4_HLOC2km"
LAB2 = "NEAR4_HT4_VT4_HLOC2km"

EXP3 = "D4_500m_0730_M4_NEAR4_HT2_VT2_HLOC2km"
LAB3 = "NEAR4_HT2_VT2_HLOC2km"

EXP4 = "D4_500m_0730_M4_NEAR4_HT4_VT4_HLOC4km"
LAB4 = "NEAR4_HT4_VT4_HLOC4km"

EXP5 = "D4_500m_0730_M4_NEAR4_HT2_VT2_HLOC4km"
LAB5 = "NEAR4_HT2_VT2_HLOC4km"


nexp = 5

EXP1 = "D4_500m_TEST_DEFAULT_0604_MEAN"
LAB1 = "DEFAULT"


#EXP2 = "D4_500m_0730_M4_NEAR4_HT4_VT4_HLOC4km"
#LAB2 = "NEAR4_HT4_VT4_HLOC4km"

EXP2 = "D4_500m_0730_M4_NEAR4_HT4_VT4_HLOC4km_MELT"
LAB2 = "NEAR4_HT4_VT4_HLOC4km_MELT"

EXP3 = "D4_500m_0730_M4_NEAR4_HT2_VT2_HLOC4km_MELT"
LAB3 = "NEAR4_HT2_VT2_HLOC4km_MELT"

EXP4 = "D4_500m_0730_M4_NEAR4_HT4_VT4_HLOC2km_MELT"
LAB4 = "NEAR4_HT4_VT4_HLOC2km_MELT"


#EXP5 = "D4_500m_0730_M4_NEAR4_HT2_VT2_HLOC2km_MELT"
#LAB5 = "NEAR4_HT2_VT2_HLOC2km_MELT"

EXP5 = "D4_500m_0730_M4_NEAR4_HT2_VT2_HLOC4km_MELT"
LAB5 = "NEAR4_HT2_VT2_HLOC2km_MELT"


EXP7 = "D4_500m_0730_ORGTHIN_HT4_VT4_HLOC2km"
LAB7 = "ORGTHIN_HT4_VT4_HLOC2km"

EXP4 = "D4_500m_0730_M4_NEAR4_HT4_VT4_HLOC4km_MELT_NOBS100"
LAB4 = "D4_500m_0730_M4_NEAR4_HT4_VT4_HLOC4km_MELT_NOBS100"


nexp = 7
EXP2 = "D4_500m_0730_M4_NEAR4_HT2_VT1_HLOC4km_MELT"
LAB2 = "NEAR4_HT2_VT1_HLOC2km_MELT"

EXP3 = "D4_500m_0730_M4_NEAR4_HT4_VT1_HLOC4km_MELT"
LAB3 = "NEAR4_HT4_VT1_HLOC2km_MELT"

EXP4 = "D4_500m_0730_M4_NEAR4_HT6_VT1_HLOC4km_MELT"
LAB4 = "NEAR4_HT6_VT1_HLOC2km_MELT"

EXP5 = "D4_500m_0730_M4_NEAR4_HT1_VT2_HLOC4km_MELT"
LAB5 = "NEAR4_HT1_VT2_HLOC4km_MELT"

EXP6 = "D4_500m_0730_M4_NEAR4_HT1_VT4_HLOC4km_MELT"
LAB6 = "NEAR4_HT1_VT4_HLOC4km_MELT"

EXP7 = "D4_500m_0730_M4_NEAR4_HT1_VT6_HLOC4km_MELT"
LAB7 = "NEAR4_HT1_VT6_HLOC4km_MELT"

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


#time = datetime(2019,9,3,2,50,0)
#time = datetime(2019, 6, 10, 8, 20, 0)
stime = datetime( 2019, 8, 24, 15, 0, 30)
#stime = datetime( 2019, 8, 24, 15, 30, 30)
etime = datetime( 2019, 8, 24, 16, 0, 0)
#etime = datetime( 2019, 8, 24, 15, 50, 0)

#etime = datetime( 2019, 8, 24, 15, 10, 0)

#etime = datetime( 2019, 8, 24, 15, 30, 0)
#stime = datetime( 2019, 8, 24, 15, 0, 30)

#stime = datetime( 2019, 8, 24, 15, 0, 30)
#etime = datetime( 2019, 8, 24, 15, 30, 0)

main( INFO, stime=stime, etime=etime, theight=theight, thrs_dbz=thrs_dbz )

