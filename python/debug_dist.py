import numpy as np

def dist( lon1=130.0, lat1=30.0, lon2=130.0, lat2=30.0 ):
    lon1 = lon1 * np.pi / 180.0
    lat1 = lat1 * np.pi / 180.0

    lon2 = lon2 * np.pi / 180.0
    lat2 = lat2 * np.pi / 180.0

    cosd = np.array( [ np.sin( lat1 ) * np.sin( lat2 ) + np.cos( lat1 ) * np.cos( lat2 ) * np.cos( lon2 - lon1 ) ] )

    cosd[ cosd > 1.0 ] = 1.0
    cosd[ cosd < -1.0 ] = -1.0

    dist = np.arccos( cosd[0] ) * 6371.3e3

    
#    dlon = lon1 - lon2 
#    dlat = lat1 - lat2 
#    dist = 2.0 * np.arcsin( np.sqrt( np.square( np.sin( dlat*0.5 ) ) + np.cos( lat1 ) * np.cos( lat2 ) * np.square( np.sin( dlon*0.5 ) ) ) ) * 6371.3e3
#    print( dist )

    return( dist )

if __name__ == "__main__":


    lon1 = 140.19128
    lon2 = 139.0267 
    
    lat1 = 36.33150465274322
    #lat1 = 35.38768582820228
    lat2 = lat1

#    lon1 = 140.19128
#    lon2 = lon1
#    lat1 = 36.33150465274322
#    lat2 = 35.38768582820228

    dist_ = dist( lon1=lon1, lat1=lat1, 
                  lon2=lon2, lat2=lat2 )

    gx = 260
    print( "dist:{0:.2f}km, num grid:{1:}, grid spacing:{2:.2f}km".format( dist_*0.001, gx,  dist_/(gx-1)*0.001 ) )
#    read_mask_full()
