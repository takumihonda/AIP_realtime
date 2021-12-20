'reinit'

'set display color white'
'c'

'set grads off'

'set mpdset hires'

'open /data11/honda/SCALE-LETKF/JMA/ctl/2020jmaradar.ctl'

lons = 139.027 
lone = 140.191
lats = 35.3877 
late = 36.3315

zlev = 2000

clevs = "0.5 1 5 10 20 30 50 80"

time = "15:30Z31JUL2020"
'set time 'time

'set lon 'lons' 'lone
'set lat 'lats' 'late

'set lev 'zlev

'set clevs 'clevs

'set gxout shade2'
'd rain'

'set grid off'
'set xlab off'
'set ylab off'
'set frame off'

'cbar '

