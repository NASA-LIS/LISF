#!/usr/local/bin/bash
echo $1
exec<"global_name.txt"
while read line
do 
    #echo $line
    line1="\b"$line"\b"
    line2=$line"_412_l"
    sed -i "s/$line1/$line2/g" $1
done
