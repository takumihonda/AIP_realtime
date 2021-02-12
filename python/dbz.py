import numpy as np

def dbz2rain( dbz ):
    B_ = 111.1
    beta_ = 1.664
    const_ = np.log10( B_ )
    #return( np.power( 10, ( dbz*0.1 - const_ ) / beta_ ) )
    return( np.power( np.power( 10, dbz*0.1) / 200.0, 1.0/1.6 ) ) 

def main():

    
    dbz = 10.0
    dbz = 20.0
    
    dbz_l = np.arange( 10, 70, 10 )
    for dbz in dbz_l:
    
       z = np.power( 10.0, dbz/10.0 )
    
       print( "DBZ:{0:.1f}, Rainfall:{1:.1f}, Z:{2:.1f}".format(dbz, dbz2rain( dbz ), z ))


main()
