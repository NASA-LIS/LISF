" ex -s -S config_dev.vim <file>
1s/#latex: //
/#latex: }
s/#latex: //
mark a
'a,$s/#latex://
'a,$s/BEGIN_DESCRIPTION//
'a,$s/END_DESCRIPTION//
'a,$s/BEGIN_CONFIG/\\begin{Verbatim}[frame=single]/
'a,$s/END_CONFIG/\\end{Verbatim}/
"'a,$s/BEGIN_DEVELOPMENT_ONLY//
"'a,$s/END_DEVELOPMENT_ONLY//
wq
