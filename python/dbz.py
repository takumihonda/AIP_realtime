import numpy as np

dbz = 10.0
dbz = 20.0

dbz_l = np.arange( 0, 30, 5 )
for dbz in dbz_l:

   z = np.power( 10.0, dbz/10.0 )

   print( "DBZ:{0:}, Z:{1:}".format(dbz, z ))
