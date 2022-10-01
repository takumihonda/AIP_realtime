import os
import sys
import numpy as np 
from netCDF4 import Dataset, default_fillvals
from datetime import datetime, timedelta
import pathlib
import glob
import tarfile

quick = False
quick = True

path = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825"
path_a = "/data_ballantine02/miyoshi-t/amemiya/SCALE-LETKF-rt-archive/result/ope/d4_500m"

ft_sec = 30*60
ft = timedelta( seconds=ft_sec )

# JST
stime = datetime( 2020, 8,  24, 16, 0 )
etime = datetime( 2020, 9,  7,  9, 0 )

dtsec = 30
dt = timedelta( seconds=dtsec )

# Figure x range
stime_ = datetime( 2020, 8, 24, 12, 0 )
stime_ = datetime( 2020, 8, 25, 0, 0 )
etime_ = datetime( 2020, 9, 7,  0, 0 )

# JMA
rmin = 30.0 # mm/h
#rmin = 20.0 # mm/h
#rmin = 10.0 # mm/h

def get_JMA_rmax( rmin=30.0 ):

    of = "JMA_rmax_rmin{0:.0f}.npz".format( rmin )
    data = np.load( of, allow_pickle=True )
    
    rmax_l = data["rmax_l"]
    rarea_l = data["rarea_l"]
    time_l = data["time_l"]
    
    # UTC2JST
    time_l += timedelta( hours=9 )

    return( rmax_l, rarea_l, time_l )

def untar():
    tfile_l = glob.glob( os.path.join( path_a, "202009*/dafcst/nc*.tar.gz" ) )
    tfile_l += glob.glob( os.path.join( path_a, "202008[2,3]*/dafcst/nc*.tar.gz" ) )

    #print( tfile_l )
    for tf in tfile_l:
        print( tf )
        with tarfile.open( tf, 'r' ) as tar:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner) 
                
            
            safe_extract(tar, os.path.join(path,"tmp"))


def get_leadtime( file_l ):

    ofn = os.path.join( path, "result_leadtime/leadtime.npz" )

    try:   
       data = np.load( ofn, allow_pickle=True ) 
       return( data["ftime_l"], data["lt_l"] )

    except:
       print( "Get time stamps now" )
       ftime_l, lt_l = prep_lt()
   
       for j, fn in enumerate( file_l ):
           if j % 1000 == 0:
              print( j, len(file_l), fn )
           p = pathlib.Path( fn )
           nc = Dataset( fn, "r", format="NETCDF4" )
           nc.set_auto_mask( False )
           var_ = nc.variables["Reflectivity"][:]
           vmax = np.max( var_ )
           if vmax == default_fillvals['f4']:
              print( "hit" )
              continue
       
           fn_time_ =  os.path.splitext(os.path.basename(fn) )[0] 
           fn_time = datetime( year   = int( fn_time_[0:4] ), 
                               month  = int( fn_time_[4:6] ),
                               day    = int( fn_time_[6:8] ), 
                               hour   = int( fn_time_[9:11] ), 
                               minute = int( fn_time_[11:13] ), 
                               second = int( fn_time_[13:15] ) )
           f_time = datetime.fromtimestamp( p.stat().st_mtime ) 
   
           dt_ = ( fn_time - stime ).total_seconds() 
           idx_ = int( dt_ / dtsec )
           lt_l[idx_] = ( fn_time + ft - f_time ).total_seconds()
   
       ftime_l = np.array( ftime_l )
       lt_l = np.array( lt_l )

       np.savez( ofn, ftime_l=ftime_l, lt_l=lt_l ) 
       return( ftime_l, lt_l )

def prep_lt():

    ftime_l = []
    lt_l = []
    
    
    time = stime
    while time <= etime:
    
       ftime_l.append( time ) # initial time
       lt_l.append( np.nan )
    
       time += dt

    return( ftime_l, lt_l )

def plot( ftime_l, lt_l ):

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':1.5, 'edgecolor':'w' }

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.dates as mdates
    
    
    fig, ((ax1, ax2)) = plt.subplots( 1, 2, figsize=( 16, 6.5 ),
                           gridspec_kw={
                           'width_ratios': [3, 1],
                           } )


    ##fig.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.95, )
    #fig.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, )
    fig.subplots_adjust( left=0.05, bottom=0.1, right=0.99, top=0.95, 
                         wspace=0.2 )
    
    # lead time
    lt_l = lt_l / 60.0 # minute
#    print( len(lt_l) )
#    print( lt_l[10000:10100] )
    ax1.plot( ftime_l, lt_l, color='k', lw=1.0 )

    ymin_ = 0
    ymax_ = 30
    dy_ = 3
    ax1.set_ylim( ymin_, ymax_ )
    
    ylab = "Forecast lead time (min)"
    ax1.set_ylabel( ylab, fontsize=12 )
    
    xlab = "Forecast initial time (JST)"
    ax1.set_xlabel( xlab, fontsize=12 )
    
    ylevs = np.arange( ymin_, ymax_+dy_, dy_ )
    ax1.set_yticks( ylevs )

    xlabs = []
    time_ = stime_
    while time_ <= etime_:
       xlabs.append( time_ )
       time_ += timedelta( hours=12 )
    
    ax1.hlines( y=ylevs, xmin=stime_, xmax=etime_, color='gray', ls='dashed', lw=0.5)
    
    tit = "Lead time of 30-min forecasts"
    ax1.text( 0.5, 1.01, tit,
              fontsize=12, transform=ax1.transAxes,
              horizontalalignment='center',
              verticalalignment='bottom' )
    
    ax1.set_xlim( stime_, etime_ )

    ax1.xaxis.set_major_locator( mdates.HourLocator(interval=24) )
    ax1.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M\n%m/%d') )
 
    xticks = []
    time_ = stime_
    while time_ <= etime_:
       xticks.append( time_ )
       time_ += timedelta( hours=24 )
   
    ax1.set_xticks( xticks )

    ax1.vlines( x=ftime_l[ np.isnan( lt_l ) ], ymin=ymin_, ymax=ymax_, color='gray', ls='dashed', lw=0.01, alpha=0.5 )


    ax2_ = ax1.twinx()

    ymin2_ = 0.0
    ymax2_ = 18.0
    ylevs2 = [ 0, 1, 2, 3, 4, 5, 6, ]

    ymin2_ = 0.0
    ymax2_ = 80.0
    ylevs2 = [ 0, 8, 16, 24, 32, 40 ]

    rmin = 30.0
    rmin_l = [ 1.0, 20.0, ]
    cc_l = [ "cyan", "b", ]
    for k, rmin in enumerate( rmin_l ):

        rmax_l, rarea_l, jtime_l = get_JMA_rmax( rmin=rmin )
        ax2_.plot( jtime_l, rarea_l/100.0, color=cc_l[k], lw=0.5, 
            label='>={0:.0f}mm h$^{{-1}}$'.format( rmin ) )

        ax2_.set_ylim( ymin2_, ymax2_ )
        ax2_.set_xlim( stime_, etime_ )
        ax2_.set_yticks( ylevs2 )
        ax2_.tick_params( axis='y', labelsize=8 ) 
        ax2_.yaxis.label.set_color( 'b' )
        ax2_.set_xlim( stime_, etime_ )

    ax2_.legend( bbox_to_anchor=( 0.05, 0.1), 
                 loc='lower left', fontsize=9, ).get_frame().set_alpha( 1.0 )

    #ylab2 = 'Precipitation area from JMA radar (x10$^2$km$^2$)\n(where >{0:.0f}mm h$^{{-1}}$)'.format( rmin )
    ylab2 = 'Precipitation area from JMA radar (x10$^2$km$^2$)'
    ax2_.set_ylabel( ylab2, fontsize=9 )
    ax2_.xaxis.set_major_locator( mdates.HourLocator(interval=24) )
    ax2_.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M\n%m/%d') )
    ax2_.tick_params(axis='y', colors='b')
    ax2_.yaxis.set_label_coords( 1.03, 0.15 )


    # plot histogram
    
    ymin = 0.0
    ymax = 30.0
    dy = 3.0
 
#    ft_l = 30.0 - lt_l
#    ft_l[ ( ft_l > 30.0 ) | ( ftime_l < stime_ ) | ( ftime_l > etime_) ] = np.nan 
#    ft_l = 30.0 - ft_l # lead time
    ft_l = np.copy( lt_l )
    ft_l[ ( ft_l <= 0.0 ) | ( ftime_l < stime_ ) | ( ftime_l > etime_) ] = np.nan
    

    fac = 0.001
    wgt_l = np.ones( ft_l.shape ) * fac   
   
    ylabs = np.arange( ymin, ymax+dy, dy )
    
    bins = int( ( ymax - ymin ) / 0.5 )
   
    y, x, _ = ax2.hist( ft_l[ ~np.isnan(ft_l) ], range=(ymin,ymax),
                 bins=bins, rwidth=0.8, orientation="horizontal", 
                 weights=wgt_l[ ~np.isnan(ft_l) ] )
 
#    xmax = x.max()
    xmax = 19000 * fac
    
    # power of 10
    pow10 = int( np.log10( 1.0/fac ) ) 
    xlab = r'Count (x10$^' + str( pow10 ) + '$)' 
    ax2.set_xlabel( xlab, fontsize=12 )
    
    #ylab = "Process time for 30-min forecast (min)"
    ax2.set_ylabel( ylab, fontsize=12 )
    
    ax2.set_ylim( ymin, ymax )
    ax2.set_yticks( ylabs )
    ax2.set_xticks( np.arange( 0, xmax, 2000*fac ) )
    ax2.set_xlim( x.min(), xmax )
    
#    yobs = 10.0 / 60.0
#    ax2.hlines( y=yobs, xmin=0, xmax=xmax, color='r', ls='solid', lw=1.0)
#    
#    arrow_dict = dict(arrowstyle = "->", color = "r", )
#    
#    ax2.annotate("Obs\nreceived(10s)", xy=(8, yobs), size=10, xytext=( 10, yobs+1.5 ),
#                color="r", arrowprops=arrow_dict )
    
    ylabs_ = np.arange( ymin, ymax, 1 )
    ax2.hlines( y=ylabs_, xmin=0, xmax=xmax, color='gray', ls='dotted', lw=0.3)
 

    total = len( ft_l[ ~np.isnan(ft_l) ] )
    total_nan = len( ft_l[ np.isnan(ft_l) ] )
    rdtime = timedelta( seconds=30*total )

    period = "Total: {0:}\n(active for {1:})\nFrom {2:}\nto {3:}".format( total, rdtime,
                                    stime.strftime('%H:%M:%SJST %m/%d'),
                                    etime.strftime('%H:%M:%SJST %m/%d'))
    ax2.text( 0.95, 0.35, period,
              fontsize=10, transform=ax2.transAxes,
              ha='right', va="top" )
    
    #tit2 = "Process time for 30-min forecast (min)"
    tit2 = "Forecast lead time (min)"
    ax2.text( 0.5, 1.01, tit2,
              fontsize=12, transform=ax2.transAxes,
              horizontalalignment='center',
              verticalalignment='bottom' )
    
    pnum_l = [ "(a)", "(b)" ]
    ax_l = [ ax1, ax2 ] 
    for i, ax in enumerate( ax_l ):
        ax.text( 0.01, 0.99, pnum_l[i],
                  fontsize=11, transform=ax.transAxes,
                  ha='left',
                  va='top', 
                  bbox=bbox )
    


    ofig = "2p_realtime_leadtime.png"
    
    print( ofig )
    if quick:
       plt.show()
    else:
#       odir = "."
#       os.makedirs( odir, exist_ok=True)
       plt.savefig( ofig,
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')
    

def main():

    #untar()
    file_l = glob.glob( os.path.join( path, "tmp/*.nc" ) )
    ftime_l, lt_l = get_leadtime( file_l, )

    plot( ftime_l, lt_l )
############

main()
