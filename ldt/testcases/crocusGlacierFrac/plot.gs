'reinit'
'xdfopen ldt_orig.xdf'
'xdfopen ldt_cro.xdf'

*'set stat on'
'set gxout grfill'
*'set mpdset mres'
'set grads off'
'set x 1 685'

'set vpage 0 5.5 4.25 8.5'
'd tbot.1'
'cbar'
'draw title ldt-orig'

'set vpage 5.5 11 4.25 8.5'
'set grads off'
'set gxout grfill'
'd tbot.2'
'cbar'
'draw title ldt-cro'

'set vpage 0 5.5 0 4.25'
'set grads off'
'set gxout grfill'
'd tbot.1-tbot.2'
'cbar'
'draw title ldt-orig - ldt-cro'

'printim ./fig/tbot.png x1200 y900 white'
*'set vpage 5.5 11 0 4.25'
*'set grads off'
*'set gxout scatter'
*'d tbot.1;tbot.2'
*'draw title ldt-orig vs. ldt-cro'
*'draw xlab ldt-orig'
*'draw ylab ldt-cro'

*say 'Hit enter to continue'
*pull args





