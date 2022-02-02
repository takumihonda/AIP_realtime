import os
import sys
from datetime import datetime, timedelta

import numpy as np

data_path = "../../dat4figs_JAMES/Fig10"
os.makedirs( data_path, exist_ok=True )

USE_ARCH_DAT = True
#USE_ARCH_DAT = False

quick_hist = False
quick_bar = True
#quick_bar = False


def d4_computation_time_nparray( top='' ):
    dirs = [ f.name for f in os.scandir( top ) ] #if f.is_file() ]

    path_l = []
    ftimes = []
    ctimes = []

    # Prepare file path list
    for dir_ in dirs:
        path_l.append( os.path.join( top, dir_, ) )
 
    scale_l = []

    # Get computation time for SCALE
    for path in path_l:
          
        if not os.path.isfile( path ):
           break

        with open( path ) as f:
             lines = f.readlines()

             for l in lines:
                 if '[Info:fcst] End forecast' in l:
                    data = l.split()
                    try:
                       ftimes.append( float( data[7] ) )
                    except:
                       print( "Failed", data )
                 elif '[Info:DA]' in l:
                    data = l.split()
                    try:
                       ctimes.append( float( data[6] ) )
                    except:
                       print( "Failed", data )

                 elif '##### TIMER' in l:
                    data = l.split()
                    try:
                       tit_ =  data[3] 

                       dat_ =  float( data[5] ) 

                       if tit_ == 'SCALE':
                          scale_l.append( dat_ ) 

                    except:
                       print( "Failed", data )
    
    scale_l = np.array( scale_l )

    key_l = [ "SCALE", "READ_OBS", 
               "OBS_OPERATOR",
               "INITIALIZE",
               "INITIALIZE_OTHERS",
               "INIT_LETKF",
               "PROCESS_OBS",
               "SET_GRID",
               "READ_GUES",
               "GUES_MEAN",
               "WRITE RESTART/GRADS(GUES)",
               "DAS_LETKF",
               "ANAL_MEAN",
               "WRITE_ANAL",
               "DEALLOCATE",
               "WRITE RESTART/GRADS(ANAL)",
               "OTHERS",
               "FINALIZE",
               "JIT_GET",
              ]
   
    # prepare nan array
    iarray = np.zeros( scale_l.shape )
    iarray[:] = np.nan
    DETAIL = {}
    for key in key_l:
        if key == 'SCALE':
           DETAIL[key] = scale_l
        else:
           DETAIL[key] = np.copy( iarray )


    # Get computation time for all
    i = -1
    for path in path_l:
          
        if not os.path.isfile( path ):
           break

        with open( path ) as f:
             lines = f.readlines()

             for l in lines:
                 if '##### TIMER' in l:
                    data = l.split()
                    try:
                       tit_ =  data[3] 
                       tit4_ =  data[4] 

                       dat_ =  float( data[5] ) 

                       if tit_ == 'SCALE':
                          i += 1

                       if tit_ == "WRITE":
                          dat_ =  float( data[6] ) 
                          if tit4_ == "RESTART/GRADS(ANAL)":
                             tit_ = "WRITE RESTART/GRADS(ANAL)"
                          elif tit4_ == "RESTART/GRADS(GUES)":
                             tit_ = "WRITE RESTART/GRADS(GUES)"
                         

                       i_ = i
                       if i_ < 0:
                          i_ = 0

                       if tit_ in DETAIL:
                          DETAIL[tit_][i_] =  dat_ 
                       else:
                          DETAIL["OTHERS"][i_] = dat_ 

                    except:
                       print( "Failed", data )

                 elif '......jitdt_read_toshiba:jitget:' in l:
                    data = l.split()
                    try:
                       tit_ =  "JIT_GET"
                       dat_ =  float( data[1] ) 
                       DETAIL[tit_][i] = dat_  
                    except:
                       print( "Failed", data )


    return( ftimes, ctimes, DETAIL )

def d4_computation_time( top='', ctmax=600 ):
    dirs = [ f.name for f in os.scandir( top ) ] #if f.is_file() ]

    ftimes = []
    ctimes = []
    path_l = []

    init = []
    init_others = []

    init_letkf = []

    scale = []
    others = []
    read_obs = []
    obsope = []

    process_obs = []
    set_grid = []
    read_gues = []
    gues_mean = []
    write_restartg = []

    das_letkf = []
    anal_mean = []
    write_anal = []
    deallocate = []
    write_restarta = []
    others = []
    finalize = []

    jitget = []

    DETAIL = { "SCALE": scale,
               "READ_OBS":read_obs,
               "OBS_OPERATOR": obsope,
               "INITIALIZE": init,
               "INITIALIZE_OTHERS": init_others,
               "INIT_LETKF": init_letkf,
               "PROCESS_OBS": process_obs,
               "SET_GRID": set_grid,
               "READ_GUES": read_gues, 
               "GUES_MEAN": gues_mean,
               "WRITE RESTART/GRADS(GUES)": write_restartg,
               "DAS_LETKF": das_letkf,
               "ANAL_MEAN": anal_mean,
               "WRITE_ANAL": write_anal,
               "DEALLOCATE": deallocate,
               "WRITE RESTART/GRADS(ANAL)": write_restarta,
               "OTHERS": others,
               "FINALIZE": finalize,
               "JIT_GET": jitget,
              }

    # Prepare file path list
    for dir_ in dirs:
        fname = 'job.o' #[ f.name for f in os.scandir( os.path.join( top, dir_ ) ) ] #if f.is_file() ]
        path_l.append( os.path.join( top, dir_, fname ) )
 
    # Get computation time
    for path in path_l:
          
        if not os.path.isfile( path ):
           break

        with open( path ) as f:
             lines = f.readlines()

             for l in lines:
                 if '[Info:fcst] End forecast' in l:
                    data = l.split()
                    try:
                       ftimes.append( float( data[7] ) )
                    except:
                       print( "Failed", data )

                 elif '[Info:DA]' in l:
                    data = l.split()
                    try:
                       ctimes.append( float( data[6] ) )
                    except:
                       print( "Failed", data )

                 elif '##### TIMER' in l:
                    data = l.split()
                    try:
                       tit_ =  data[3] 
                       tit4_ =  data[4] 

                       dat_ =  float( data[5] ) 

                       if tit_ == "WRITE":
                          dat_ =  float( data[6] ) 
                          if tit4_ == "RESTART/GRADS(ANAL)":
                             tit_ = "WRITE RESTART/GRADS(ANAL)"
                          elif tit4_ == "RESTART/GRADS(GUES)":
                             tit_ = "WRITE RESTART/GRADS(GUES)"
                         


                       if tit_ in DETAIL:
                          DETAIL[tit_].append( dat_ ) 
                       else:
                          DETAIL["OTHERS"].append( dat_ ) 

                    except:
                       print( "Failed", data )

                 elif '......jitdt_read_toshiba:jitget:' in l:
                    data = l.split()
                    try:
                       tit_ =  "JIT_GET"
                       dat_ =  float( data[1] ) 
                       DETAIL[tit_].append( dat_ ) 
                    except:
                       print( "Failed", data )

    for key in DETAIL.keys():
        DETAIL[key] = np.array( DETAIL[key] )

    return( ftimes, ctimes, DETAIL )

def plot_hist( key="", dat=np.array([]) ):

    import matplotlib.pyplot as plt
    from scipy import stats
    
    xmin = 0
    xmax = 60
    

    # Scott's choise
    #h = 3.5 * np.std( dat, ddof=1 ) / np.power( dat.size, 1.0/3.0)
    #bins = int( ( xmax - xmin ) / h )
    # Square-root choice
    bins = int( np.sqrt( dat.size ) )
    
    fig, ax = plt.subplots( 1, 1, figsize=(6,4) )
    fig.subplots_adjust( left=0.15, bottom=0.15, right=0.95, top=0.92, )
    
    rn, rbins, rpatches = ax.hist( dat, range=(xmin, xmax), bins=bins, alpha=0.6 )
    
    imode = np.argmax( rn )
    mode = np.mean( rbins[imode:imode+2] )
    mean = np.mean( dat )
    #print( len(rn), len(rbins), mode )
    
    lw = 1.0
    ymin = 0.0
    ymax = 4000 #dat_.size
    ls = 'dashed'
    color = 'b'
    ax.vlines( x=mode, ymin=ymin, ymax=ymax,
               linewidths=lw, linestyles=ls, color=color )
    
    color = 'k'
    ax.vlines( x=mean, ymin=ymin, ymax=ymax,
               linewidths=lw, linestyles=ls, color=color )

    text_ = 'Mean:{0:.3f} s\nMode:{1:.3f} s\nN={2:}'.format( mean, mode, dat.size )
    ax.text( 0.99, 0.99, text_,
              fontsize=12, transform=ax.transAxes,
              ha='right',
              va='top' )
    
    tit_ = key
    ax.text( 0.5, 1.01, tit_,
              fontsize=12, transform=ax.transAxes,
              ha='center',
              va='bottom' )
    
    ax.set_xlim( xmin, xmax )
    ax.set_ylim( ymin, ymax )
    
    xlab = 'Computation time (s)'
    ylab = 'Frequency'
    ax.set_xlabel( xlab, fontsize=11)
    ax.set_ylabel( ylab, fontsize=11)
 
    key_ =  key.replace( ' ', '_' ).replace( '/', '_' ) #.replace( '(', '_' ).replace( ')')
    ofig = 'png/1p_d4_{0:}.png'.format( key_ )

    print( ofig )
    if quick_hist:
       plt.show()
    else:
       plt.savefig( ofig,
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')
       
    return( mode, mean )

def plot_bar_2p( dic={}, ftimes=np.array([]) ):

    import matplotlib.pyplot as plt
    fig, ( ax1,ax2 ) = plt.subplots( 1, 2, figsize=(6,4) )
#    fig.subplots_adjust( left=0.15, bottom=0.05, right=0.5, top=0.92, )
    fig.subplots_adjust( left=0.15, bottom=0.06, right=0.95, top=0.92,
                          wspace=0.3, hspace=0.05)

    ax1.set_xlim( 0, 2.0 )
    width1 = 0.8

    #c_l = [ 'firebrick', 'dodgerblue', 'limegreen', 'gold' ]
    #c_l = [ 'dodgerblue', 'firebrick', 'forestgreen', 'goldenrod' ]
    c_l = [ 'dodgerblue', 'firebrick', 'gray', 'goldenrod', 'k' ]
    #c_l = [ 'cyan', 'magenta', 'y', 'k' ]
    acm = 0.0
    for i, key in enumerate( dic.keys() ):
        lab = key
        if lab == 'OBS':
           lab = 'Obs pre-\nprocessing'
        elif lab == 'DATA TRANSFER':
           lab = 'Memory copy'
        elif lab == 'JIT-DT':
           continue
        ax1.bar( 1.0, dic[key], bottom=acm,
               label=lab, color=c_l[i], width=width1 )
        acm += dic[key]

#    ax.legend( fontsize=12, bbox_to_anchor=(1.01, 1.00) )
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend( reversed(handles), reversed(labels), bbox_to_anchor=(1.01, 1.00),
               fontsize=12 )
    ax1.set_ylabel( 'Computation time (s)', fontsize=12 )
    #ax.set_xlim( 0, 1.0 )
    yticks = np.arange( 0, 22, 2 )
    ax1.set_ylim( 0, 20.0 )
    ax1.set_yticks( yticks )

    ax1.tick_params( axis='x', which='both', 
                     bottom=False, top=False, 
                     labelbottom=False )

    ax1.hlines( xmin=0, xmax=2, y=np.arange( 4, 20, 4 ), lw=1.0, linestyle='dashed', 
                color='gray', alpha=0.5 )


    ax2.set_ylim( 0, 151.0 )
    ax2.set_xlim( 0, 2.0 )
    ax2.hlines( xmin=0, xmax=2, y=[60, 120], lw=1.0, linestyle='dashed', 
                color='gray', alpha=0.5 )
    width2 = 0.8
    ax2.bar( 1, np.mean(ftimes), label="30-min forecast", width=width2,
             color='dodgerblue' )
    print( "std:", np.std( ftimes, ddof=1 ), len( ftimes ) )

    ax2.tick_params( axis='x', which='both', 
                     bottom=False, top=False, 
                     labelbottom=False )

    ax_l = [ ax1, ax2 ]

    tit_l = [ "Data assimilation",
              "30-min forecast" ]
    pnum_l = [ "(a)", "(b)" ]
    for i, ax in enumerate( ax_l ):
       ax.text( 0.5, 1.01, tit_l[i],
                 fontsize=12, transform=ax.transAxes,
                 ha='center',
                 va='bottom' )

       ax.text( 0.0, 1.01, pnum_l[i],
                 fontsize=10, transform=ax.transAxes,
                 ha='left',
                 va='bottom' )

    ofig = 'pdf/Fig10.pdf'
    print( ofig )
    if quick_bar:
       plt.show()
    else:
       plt.savefig( ofig,
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')


def plot_bar_2p_scale( dic={}, ftimes=np.array([]), dic2={} ):

    import matplotlib.pyplot as plt
    fig, ( ax1,ax2 ) = plt.subplots( 1, 2, figsize=(6,4) )
#    fig.subplots_adjust( left=0.15, bottom=0.05, right=0.5, top=0.92, )
    fig.subplots_adjust( left=0.15, bottom=0.06, right=0.95, top=0.92,
                          wspace=0.3, hspace=0.05)

    ax1.set_xlim( 0, 3.0 )
    width1 = 0.8

    #c_l = [ 'firebrick', 'dodgerblue', 'limegreen', 'gold' ]
    #c_l = [ 'dodgerblue', 'firebrick', 'forestgreen', 'goldenrod' ]
    c_l = [ 'dodgerblue', 'firebrick', 'gray', 'goldenrod', 'k' ]
    #c_l = [ 'cyan', 'magenta', 'y', 'k' ]
    acm = 0.0
    for i, key in enumerate( dic.keys() ):
        lab = key
        if lab == 'OBS':
           lab = 'Obs pre-\nprocessing'
        elif lab == 'DATA TRANSFER':
           lab = 'Memory copy'
        elif lab == 'JIT-DT':
           continue
        ax1.bar( 1.0, dic[key], bottom=acm,
               label=lab, color=c_l[i], width=width1 )

        acm += dic[key]

    acm2 = 0.0
    for i, key in enumerate( dic2.keys() ):
        lab = key
        if lab == 'OBS':
           lab = 'Obs pre-\nprocessing'
        elif lab == 'DATA TRANSFER':
           lab = 'Memory copy'
        elif lab == 'JIT-DT':
           continue
        print( "check", dic2[key] )
        ax1.bar( 2.0, dic2[key], bottom=acm2,
               label=None, color=c_l[i], width=width1 )

        acm2 += dic[key]



#    ax.legend( fontsize=12, bbox_to_anchor=(1.01, 1.00) )
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend( reversed(handles), reversed(labels), bbox_to_anchor=(1.01, 1.00),
               fontsize=12 )
    ax1.set_ylabel( 'Computation time (s)', fontsize=12 )
    #ax.set_xlim( 0, 1.0 )
    yticks = np.arange( 0, 22, 2 )
    ax1.set_ylim( 0, 20.0 )
    ax1.set_yticks( yticks )

    ax1.tick_params( axis='x', which='both', 
                     bottom=False, top=False, 
                     labelbottom=False )

    ax1.hlines( xmin=0, xmax=2, y=np.arange( 4, 20, 4 ), lw=1.0, linestyle='dashed', 
                color='gray', alpha=0.5 )


    ax2.set_ylim( 0, 151.0 )
    ax2.set_xlim( 0, 2.0 )
    ax2.hlines( xmin=0, xmax=2, y=[60, 120], lw=1.0, linestyle='dashed', 
                color='gray', alpha=0.5 )
    width2 = 0.8
    ax2.bar( 1, np.mean(ftimes), label="30-min forecast", width=width2,
             color='dodgerblue' )
    print( "std:", np.std( ftimes, ddof=1 ), len( ftimes ) )

    ax2.tick_params( axis='x', which='both', 
                     bottom=False, top=False, 
                     labelbottom=False )

    ax_l = [ ax1, ax2 ]

    tit_l = [ "Data assimilation",
              "30-min forecast" ]
    pnum_l = [ "(a)", "(b)" ]
    for i, ax in enumerate( ax_l ):
       ax.text( 0.5, 1.01, tit_l[i],
                 fontsize=12, transform=ax.transAxes,
                 ha='center',
                 va='bottom' )

       ax.text( 0.0, 1.01, pnum_l[i],
                 fontsize=10, transform=ax.transAxes,
                 ha='left',
                 va='bottom' )

#    ofig = 'png/2p_d4_bar_scale.png'
    print( ofig )
    if quick_bar:
       plt.show()
    else:
       plt.savefig( ofig,
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')




def plot_bar( dic={} ):

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots( 1, 1, figsize=(5,5) )
    fig.subplots_adjust( left=0.15, bottom=0.05, right=0.5, top=0.92, )

    #c_l = [ 'firebrick', 'dodgerblue', 'limegreen', 'gold' ]
    #c_l = [ 'dodgerblue', 'firebrick', 'forestgreen', 'goldenrod' ]
    c_l = [ 'dodgerblue', 'firebrick', 'gray', 'goldenrod', 'k' ]
    #c_l = [ 'cyan', 'magenta', 'y', 'k' ]
    acm = 0.0
    for i, key in enumerate( dic.keys() ):
        lab = key
        if lab == 'OBS':
           lab = 'Obs pre-\nprocessing'
        elif lab == 'DATA TRANSFER':
           lab = 'Memory copy'
        elif lab == 'JIT-DT':
           continue
        ax.bar( '', dic[key], bottom=acm,
               label=lab, color=c_l[i] )
        acm += dic[key]

#    ax.legend( fontsize=12, bbox_to_anchor=(1.01, 1.00) )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend( reversed(handles), reversed(labels), bbox_to_anchor=(1.01, 1.00),
               fontsize=13 )
    ax.set_ylabel( 'Computation time (s)', fontsize=12 )
    #ax.set_xlim( 0, 1.0 )
    yticks = np.arange( 0, 32, 2 )
    ax.set_ylim( 0, 31.0 )
    ax.set_yticks( yticks )

    ofig = 'png/1p_d4_bar.png'
    print( ofig )
    if quick_bar:
       plt.show()
    else:
       plt.savefig( ofig,
                    bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
       plt.close('all')

####

SUM = { "SCALE": 0.0,
        "LETKF": 0.0,
        "OBS": 0.0,
#        "DATA TRANSFER": 0.0,
        "JIT-DT": 0.0,
      }

fn_sum = '{0:}/SUM.npz'.format( data_path, )
fn_ftimes = '{0:}/ftimes.npz'.format( data_path, )
if not USE_ARCH_DAT:
   
   
   top = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/log_from_amemiya/d4_500m/exp'
   top = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/amemiya/d4_500m'
   
   top_test = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/data/D4_500m_TEST_DEFAULT_0708_NOBS100_NEAR_HV4/exp/3008084_cycle_20190824150000'
   
   #dtime_max = 1000
   
   ftimes, ctimes, DETAIL = d4_computation_time_nparray( top=top,  )
   ftimes_test, ctimes_test, DETAIL_test = d4_computation_time_nparray( top=top_test,  )
   
   #print( DETAIL["DAS_LETKF"][0:5], DETAIL["WRITE_ANAL"][0:5])
   
   #ftimes, ctimes, DETAIL = d4_computation_time( top=top,  )
   ctimes = np.array( ctimes )
   print( '{0:} average: {1:} (N: {2:})'.format( "cycle", np.nanmean( ctimes ), len(ctimes) ) )
   print( '{0:} average: {1:} (N: {2:})'.format( "fcst ", np.mean( ftimes ), len(ftimes) ) )
   
   print("")
   DETAIL_MODE = { }
   DETAIL_MODE_test = { }
   
   
   min_read_obs = 1.0
   max_read_obs = 30.0
   read_obs_ = DETAIL["READ_OBS"]
   
   
   dat_jit = DETAIL['JIT_GET']
   dat_jit[ ( read_obs_ < min_read_obs ) | ( read_obs_ > max_read_obs )] = np.nan
   dat_jit_ = dat_jit[ ~np.isnan(dat_jit) ]
   for key in DETAIL.keys():
       DETAIL[key][ ( read_obs_ < min_read_obs ) | ( read_obs_ > max_read_obs )] = np.nan
       time_ = np.nanmean( DETAIL[key] )
   
       dat = DETAIL[key]
       dat_ = dat[ ~np.isnan(dat) & ~np.isnan( dat_jit ) ]
       num = len( dat_ ) 
   
       if key == "READ_OBS":
          dat_ -= dat_jit_
   
       print( "#### ", key, time_, num, np.nanmax( DETAIL[key] ), np.nanmin( DETAIL[key] ) )
   
       if num > 100:
          mode_, mean_ = plot_hist( key=key, dat=dat_ )
   
          #DETAIL_MODE[key] = mode_
          DETAIL_MODE[key] = mean_
   
       else:
          print( 'Not plot ', key)
   
   
   read_obs_test = DETAIL_test["READ_OBS"]
   #dat_jit_test = DETAIL_test['JIT_GET']
   #dat_jit_test[ ( read_obs_test < min_read_obs ) | ( read_obs_test > max_read_obs )] = np.nan
   #dat_jit_test = dat_jit_test[ ~np.isnan(dat_jit_test) ]
   for key in DETAIL_test.keys():
       DETAIL_test[key][ ( read_obs_test < min_read_obs ) | ( read_obs_test > max_read_obs )] = np.nan
       time_ = np.nanmean( DETAIL_test[key] )
   
       dat = DETAIL_test[key]
       print( key, dat )
       #dat_ = dat[ ~np.isnan(dat) & ~np.isnan( dat_jit_test ) ]
       dat_ = dat[ ~np.isnan(dat) ]
       num = len( dat_ ) 
   
   #    if key == "READ_OBS":
   #       dat_ -= dat_jit_
   
       print( "#### ", key, time_, num, np.nanmax( DETAIL_test[key] ), np.nanmin( DETAIL_test[key] ) )
   
       if num > 100:
          mode_, mean_ = plot_hist( key=key, dat=dat_ )
   
          DETAIL_MODE_test[key] = mean_
   
       else:
          print( 'Not plot ', key)
   
   
   
   
   
   for key in DETAIL_MODE.keys():
       print( key )
       if key == "SCALE":
          SUM["SCALE"] += DETAIL_MODE[key]   
       elif key == "READ_OBS":
          SUM["OBS"] += DETAIL_MODE[key]   
   #    elif key == "READ_GUES" or key == "WRITE_ANAL":
   #       SUM["DATA TRANSFER"] += DETAIL_MODE[key]   
       elif key == "JIT_GET":
          SUM["JIT-DT"] += DETAIL_MODE[key]   
       else:
          SUM["LETKF"] += DETAIL_MODE[key]   
   
   
   SUM_test = { "SCALE": 0.0,
                "LETKF": 0.0,
                "OBS": 0.0,
               "JIT-DT": 0.0,
              } 
   
   for key in DETAIL_MODE_test.keys():
       if key == "SCALE":
          SUM_test["SCALE"] += DETAIL_MODE_test[key]   
       elif key == "READ_OBS":
          SUM_test["OBS"] += DETAIL_MODE_test[key]   
   #    elif key == "READ_GUES" or key == "WRITE_ANAL":
   #       SUM["DATA TRANSFER"] += DETAIL_MODE[key]   
       elif key == "JIT_GET":
          SUM_test["JIT-DT"] += DETAIL_MODE_test[key]   
       else:
          SUM_test["LETKF"] += DETAIL_MODE_test[key]   
   
   
   np.savez( fn_sum, **SUM, ftimes=ftimes )
   np.savez( fn_ftimes, ftimes=ftimes )

else:

   with np.load( fn_sum, allow_pickle=True ) as npz:
      for key in SUM.keys():
          SUM[key] = npz[key]
   ftimes = np.load( fn_ftimes, allow_pickle=True )['ftimes']


print( SUM )
#print( DETAIL_MODE )
#print( SUM_test )
#print( DETAIL_MODE_test )
#sys.exit()
#plot_bar( dic=SUM )
plot_bar_2p( dic=SUM, ftimes=ftimes )
#plot_bar_2p_scale( dic=SUM, dic2=SUM_test, ftimes=ftimes )

