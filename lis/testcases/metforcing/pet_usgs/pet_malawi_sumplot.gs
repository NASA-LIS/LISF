
  prompt 'which filenum:'
   pull filenum

  'rgbset.gs'
  'set parea 1.0 10.0 0.8 7.9'
  'set gxout grfill'
  'set mpdset hires'
  'set grads off'

  'set lon 32 36.2'
  'set lat -18 -9'
  'set ylopts 1 0 0.2'
  'set xlopts 1 0 0.2'
  'set xlint 1 '
  'set ylint 1 '

* Kristi's interval/color range:
*  'set clevs 0 10 25 50 100 150 200 250 300 350 400 450 '
*  'set ccols  59 59 57 55 53 52 51 21 22 23 25 27 29 '

* Brad's color range:
  'set clevs   0 10 50 100 150 200 250 300 400 500 1000 '
  'set ccols  0 16 62 32 35  38   42  45 48   22  25   27'

  'set t 1'
  'd sum(pet.'filenum'*86400,t=1,t=60)'
  'cbarn'

*  'draw title USGS PET (1.0 deg) -- April 1, 2011 (00Z)'

  'set strsiz 0.18'
  'draw string 8.6 7.6 (mm)'
