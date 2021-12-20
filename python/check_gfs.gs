

fn = '/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/GFS/grads/20190824120000/mean/atm0.25.ctl'

'reinit'
'open 'fn

'set grads off'

'set gxout shade2'

lons = 120
lone = 160
lats = 20
late = 60
lev = 1000
lev = 500

'set lon 'lons' 'lone
'set lat 'lats' 'late

'set lev 'lev

'd tmpprs'
*'d tmpprs(lev=1000)-tmpprs(lev=500)'

'cbarn'

