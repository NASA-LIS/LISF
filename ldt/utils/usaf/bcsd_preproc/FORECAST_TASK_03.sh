#!/bin/sh

source /usr/share/modules/init/sh
ulimit -s unlimited
module load python/GEOSpyD/Ana2019.10_py3.7

currentmon1=${1}
currentyear=${2}

python code_library/nmme_reorg_f.py ${currentmon1} ${currentyear}
