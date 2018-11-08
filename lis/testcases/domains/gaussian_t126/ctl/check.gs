'open input.ctl'
'open lmask_gfs_T126.ctl'

'set gxout grfill'

'd landmask.1-landmask.2'
'draw title landmask diff'
'cbar'
'printim plot_landmask.png png x1000 y800 white'
'c'

'reinit'

'open input.ctl'
'open vtype_gfs_T126.ctl'

'set gxout grfill'

'd landcover1.1-landcover1.2'
'draw title landcover1 diff'
'cbar'
'printim plot_landcover1.png png x1000 y1000 white'
'c'

'd landcover2.1-landcover2.2'
'draw title landcover2 diff'
'cbar'
'printim plot_landcover2.png png x1000 y1000 white'
'c'

'd landcover3.1-landcover3.2'
'draw title landcover3 diff'
'cbar'
'printim plot_landcover3.png png x1000 y1000 white'
'c'

'd landcover4.1-landcover4.2'
'draw title landcover4 diff'
'cbar'
'printim plot_landcover4.png png x1000 y1000 white'
'c'

'd landcover5.1-landcover5.2'
'draw title landcover5 diff'
'cbar'
'printim plot_landcover5.png png x1000 y1000 white'
'c'

'd landcover6.1-landcover6.2'
'draw title landcover6 diff'
'cbar'
'printim plot_landcover6.png png x1000 y1000 white'
'c'

'd landcover7.1-landcover7.2'
'draw title landcover7 diff'
'cbar'
'printim plot_landcover7.png png x1000 y1000 white'
'c'

'd landcover8.1-landcover8.2'
'draw title landcover8 diff'
'cbar'
'printim plot_landcover8.png png x1000 y1000 white'
'c'

'd landcover9.1-landcover9.2'
'draw title landcover9 diff'
'cbar'
'printim plot_landcover9.png png x1000 y1000 white'
'c'

'd landcover10.1-landcover10.2'
'draw title landcover10 diff'
'cbar'
'printim plot_landcover10.png png x1000 y1000 white'
'c'

'd landcover11.1-landcover11.2'
'draw title landcover11 diff'
'cbar'
'printim plot_landcover11.png png x1000 y1000 white'
'c'

'd landcover12.1-landcover12.2'
'draw title landcover12 diff'
'cbar'
'printim plot_landcover12.png png x1000 y1000 white'
'c'

'reinit'

'open input.ctl'
'open elev_gfs_T126.ctl'

'set gxout grfill'

'd elevation.1-elevation.2'
'draw title elevation diff'
'cbar'
'printim plot_elevation.png png x1000 y800 white'
'c'

'reinit'

'run tex.gs'
'draw title texture diff'
'cbar'
'printim plot_texture.png png x1000 y800 white'
'c'

'reinit'

'open input.ctl'
'open maxsnowalb_gfs_T126.ctl'

'set gxout grfill'

'd mxsnalbedo.1-mxsnalbedo.2'
'draw title mxsnalbedo diff'
'cbar'
'printim plot_mxsnalbedo.png png x1000 y800 white'
'c'

'reinit'

'open input.ctl'
'open tbot_gfs_T126.ctl'

'set gxout grfill'

'd tbot.1-tbot.2'
'draw title tbot diff'
'cbar'
'printim plot_tbot.png png x1000 y800 white'
'c'

