import os
import numpy as np

from datetime import datetime
from tools_AIP import read_dat, set_cmap_pawr_scat

import matplotlib.pyplot as plt
import scipy.stats as st

def main( time=datetime( 2019, 8, 24, 15, 10, 0 ), top='', exp_l=[], texp_l=[],
          time0=datetime( 2019, 8, 24, 15, 10, 0 ), height_l=[3000.0], dz=300.0,
          otyp='z', nvar='ob' ):

    bbox = { 'facecolor':'w', 'alpha':0.95, 'pad':0.5,
             'edgecolor':'w' }

    if otyp == 'z':
       otyp_ = 4001
       oerr = 5.0
       gqc = 10.0
    elif otyp == 'cz':
       otyp_ = 4004
    elif otyp == 'vr':
       otyp_ = 4002
       oerr = 3.0
       gqc = 5.0

    if nvar == 'ob':
       DIF = True
    else:
       DIF = False

    unit, cmap, levels, norm, vmin, vmax, tvar = set_cmap_pawr_scat( otyp=otyp, DIF=DIF )

    #dat = ombs

    fig = plt.figure( figsize=(10, 5) )
    fig.subplots_adjust( left=0.08, bottom=0.1, right=0.98, top=0.9,
                         wspace=0.15, hspace=0.2)
    ax_l = []
    xfig = 2
    yfig = 1
    for i in range( xfig ):
        ax_l.append( fig.add_subplot( yfig,xfig,i+1, projection=None ) )
#    ax_l.append( fig.add_subplot( 1,1,1, projection=None ) )

    size = 5.0   

    xlab = 'X (km)'
    ylab = 'Y (km)'

    xmin = -25
    xmax = 25

    for i, ax in enumerate( ax_l ):

        height = height_l[i]

        fn = os.path.join( top, exp_l[i], time0.strftime('%Y%m%d%H%M%S'),
                           'obsdep', 'obsdep_{0:}'.format( time.strftime('%Y%m%d-%H%M%S.dat')) )
    
        print( fn )
        elms, lons, lats, levs, dats, ombs, qcs, omas = read_dat( fn )
    
#        lats = lats[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
#        lons = lons[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
        dats = dats[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
        ombs = ombs[ ( np.abs( levs - height ) < dz ) & ( elms == otyp_ ) ]
    
        if nvar == 'ob':
           dat = ombs
           tit_ = 'O–B'
        elif nvar == 'obs':
           dat = dats
           tit_ = 'Obs'
        print( "debug ", np.min(dat ), np.max( dat ) )
      

#        dat = np.random.normal( 0.0, oerr, 100000 )

        sigma = np.std( dat )
        nsample = len( dat )
        h = 3.5 * sigma / np.power( nsample, 1.0/3.0 )
        bins = int( ( xmax - xmin ) / h )

        _, _, _ = ax.hist( dat, range=(xmin, xmax), bins=bins, alpha=0.6 )

        mu = np.mean( dat )
        xs = np.linspace( xmin, xmax, 100 )
        ys = st.norm.pdf( xs, mu, sigma ) * nsample * h
        ax.plot( xs, ys, color='b', lw=1.0 )

        ymin, ymax = ax.get_ylim()
        ax.vlines( x=[ -oerr*gqc, oerr*gqc], ymin=ymin, ymax=ymax,
                   colors='gray', ls='dashed', lw=1.0 )
        ax.set_ylim( ymin, ymax )
    
        ax.text( 0.5, 1.01, texp_l[i],
                va='bottom', 
                ha='center',
                transform=ax.transAxes,
                color='k', fontsize=12, )

        if i == 1:
           cvtime = time.strftime( '%H:%M:%S' )
           ax.text( 1.0, 1.01, 'Valid at {0:}'.format( cvtime ),
                   va='bottom', 
                   ha='right',
                   transform=ax.transAxes,
                   color='k', fontsize=10, )
   
        ax.text( 0.99, 0.9, r'Z={0:.1f}$\pm${1:.1f} km'.format( height*0.001, dz*0.001 ), 
                va='top', 
                ha='right',
                bbox=bbox,
                transform=ax.transAxes,
                color='k', fontsize=10, )


    tit = '{0:} {1:}'.format( tvar,  tit_ )
    fig.suptitle( tit, fontsize=14 )


    plt.show()

top = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/scale-5.4.3/OUTPUT/realtime2021/test'
exp = 'd4_500m_ac'

exp_l = [ 'd4_500m_ac', 'd4_500m_ac' ]
texp_l = [ "TEST", "TEST" ]

time0 = datetime( 2019, 8, 24, 15, 0 )
time = datetime( 2019, 8, 24, 15, 30, 0 )
time = datetime( 2019, 8, 24, 15, 20, 0 )
time = datetime( 2019, 8, 24, 15, 2, 0 )

otyp = 'vr'
#otyp = 'z'

height_l = [ 500, 1500 ]
dz = 500.0

height_l = [ 500, 5000 ]
dz = 5000.0

height_l = [ 1000, 3000 ]
height_l = [ 2000, 6000 ]
#height_l = [ 9000, 11000 ]
dz = 2000.0

nvar = 'ob'

main( time=time, top=top, exp_l=exp_l, time0=time0, otyp=otyp, texp_l=texp_l,
      height_l=height_l, nvar=nvar, dz=dz )
