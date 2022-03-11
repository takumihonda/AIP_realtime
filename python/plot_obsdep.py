import os
import numpy as np

from datetime import datetime
from tools_AIP import read_dat, set_cmap_pawr_scat

import matplotlib.pyplot as plt

quick = True

def main( time=datetime( 2019, 8, 24, 15, 10, 0 ), top='', exp_l=[], texp_l=[],
          time0=datetime( 2019, 8, 24, 15, 10, 0 ), height=3000.0, dz=300.0,
          otyp='z', nvar='ob' ):

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    if otyp == 'z':
       otyp_ = 4001
    elif otyp == 'cz':
       otyp_ = 4004
    elif otyp == 'vr':
       otyp_ = 4002

    if nvar == 'ob':
       DIF = True
    else:
       DIF = False

    unit, cmap, levels, norm, vmin, vmax, tvar = set_cmap_pawr_scat( otyp=otyp, DIF=DIF )

    #dat = ombs

    fig = plt.figure( figsize=(10, 5) )
    fig.subplots_adjust( left=0.08, bottom=0.1, right=0.93, top=0.9,
                         wspace=0.1, hspace=0.2)
    ax_l = []
    xfig = 2
    yfig = 1
    for i in range( xfig ):
        ax_l.append( fig.add_subplot( yfig,xfig,i+1, projection=None ) )
#    ax_l.append( fig.add_subplot( 1,1,1, projection=None ) )

    size = 5.0   

    xlab = 'X (km)'
    ylab = 'Y (km)'

    xmin = -60
    xmax = 60
    ymin = -60
    ymax = 60

#    ymin = 0
#    ymax = 60
#    xmin = -30
#    xmax = 30

    dx = 0.5
    dy = 0.5

    x1d = np.arange( xmin, xmax+dx, dx )
    y1d = np.arange( ymin, ymax+dy, dy )

    x2d, y2d = np.meshgrid( x1d, y1d )
    dist2d = np.sqrt( np.square( x2d ) + np.square( y2d ) )


    for i, ax in enumerate( ax_l ):
        fn = os.path.join( top, exp_l[i], time0.strftime('%Y%m%d%H%M%S'),
                           'obsdep', 'obsdep_{0:}'.format( time.strftime('%Y%m%d-%H%M%S.dat')) )
    
        print( fn )
        elms, lons, lats, levs, dats, ombs, qcs, omas = read_dat( fn )
    
        lats = lats[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
        lons = lons[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
        dats = dats[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
        ombs = ombs[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
    
        if nvar == 'ob':
           dat = ombs
           tit_ = 'O–B'
        elif nvar == 'obs':
           dat = dats
           tit_ = 'Obs'
        elif nvar == 'b':
           dat = -ombs + dats
           tit_ = 'B'
        print( "debug ", np.min(dat ), np.max( dat ) )
    
        SHADE = ax.scatter( lons*0.001, lats*0.001, c=dat, s=size,
                     norm=norm, 
                     cmap=cmap,
                     vmin=vmin,
                     vmax=vmax,
                    )
 
        ax.set_xlim( xmin, xmax )
        ax.set_ylim( ymin, ymax )
   
        CONT = ax.contour( x2d, y2d, dist2d, 
                          levels=np.arange( 10, 60, 10 ), 
                          colors='gray',
                          linestyles='dashed',
                          linewidths=1.0,
                         )
        ax.clabel( CONT, fontsize=10, fmt='%.0f km', inline=True )

        ax.text( 0.5, 1.01, texp_l[i],
                va='bottom', 
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=12, )

        if i == 0:
           ax.set_ylabel( ylab, fontsize=12 )
        ax.set_xlabel( xlab, fontsize=12 )

        if i == 1:
           pos = ax.get_position()
           cb_width = 0.006
           cb_height = pos.height*0.9
           ax_cb = fig.add_axes( [ pos.x1+0.003, pos.y0, 
                                   cb_width, cb_height] )
           cb = plt.colorbar( SHADE, cax=ax_cb, orientation='vertical',  
                              ticks=levels[::1], 
                              extend='both' )

           cvtime = time.strftime( '%H:%M:%S' )
           ax.text( 1.0, 1.01, 'Valid at {0:}'.format( cvtime ),
                   va='bottom', 
                   ha='right',
                   transform=ax.transAxes,
                   color='k', fontsize=10, )
   
           ax.text( 0.99, 0.01, r'Z={0:.1f}$\pm${1:.1f} km'.format( height*0.001, dz*0.001 ), 
                   va='bottom', 
                   ha='right',
                   bbox=bbox,
                   transform=ax.transAxes,
                   color='k', fontsize=10, )


    tit = '{0:} {1:}'.format( tvar,  tit_ )
    fig.suptitle( tit, fontsize=14 )

    ofig = "2p_obsdep_ac_{0:}_{1:}_h{2:}km_dz{3:.3f}km_v{4:}.png".format( 
                    exp_l[0],
                    exp_l[1],
                    height*0.001, 
                    dz*0.001, 
                    time.strftime('%Y%m%d%H%M%S'),
                     )

    print( ofig )
    if not quick:
       opath = "png/attenuation"
       os.makedirs( opath, exist_ok=True )
       ofig = os.path.join(opath, ofig)
       plt.savefig(ofig,bbox_inches="tight", pad_inches = 0.1)
       plt.clf()
    else:
       plt.show()


######################

top = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test'
exp = 'd4_500m_ac'

exp_l = [ 'd4_500m_ac', 'd4_500m_ac_vr' ]
texp_l = [ "TEST", "TEST (VR only)" ]

exp_l = [ 'd4_500m_ac', 'd4_500m_ac' ]
texp_l = [ "TEST", "TEST" ]

#exp_l = [ 'd4_500m_ac_z', 'd4_500m_ac_vr' ]
#texp_l = [ "TEST (Z only)", "TEST (VR only)" ]

time0 = datetime( 2019, 8, 24, 15, 0 )
time = datetime( 2019, 8, 24, 15, 30, 0 )

time0 = datetime( 2019, 8, 19, 13, 0 )
time = datetime( 2019, 8, 19, 13, 30, 0 )

otyp = 'vr'
#otyp = 'z'
height = 3000
#height = 10000
#height = 7000
#height = 1000
height = 500
height = 1000
#height = 2000
#height = 1500
#height = 2500
#height = 3500
#height = 6500
#height = 8500

nvar = 'obs'
#nvar = 'ob'
#nvar = 'b'

main( time=time, top=top, exp_l=exp_l, time0=time0, otyp=otyp, texp_l=texp_l,
      height=height, nvar=nvar )
