" ex -s -S config_mo.vim <file>
1s/#latex: //
/#latex: }
s/#latex: //
mark a
'a,$s/#latex://
'a,$s/BEGIN_DESCRIPTION//
'a,$s/END_DESCRIPTION//
'a,$s/BEGIN_CONFIG/\\begin{verbatim}/
'a,$s/END_CONFIG/\\end{verbatim}/
wq
