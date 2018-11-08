#!/bin/sh

#srcdir=.
#docdir=.
#
#echo $srcdir

list=`find ./ -name "*.F90"`

for name in $list
do 
 echo " filename $name "
 mv $name $name.old
 sed -e "s/LIS/LPS/1" $name.old > $name
# sed -e "s/ESMF_ConfigMod/ESMF_Mod/1" $name.old > $name
# mv $name.old $name
 rm -f $name.old
# filename=${name%.*}
# echo " filename1 $filename "
# filename1=${name##*/}
# filename1=${filename1%.*}
# filename2=$filename1".tex"
# perl /data1/sujay/protex $name > $filename2
done
