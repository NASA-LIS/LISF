#!/bin/sh 
if [ $# -ne 2 ]
then 
  echo "Usage: rm_comments.sh <input lvt.config> <output lvt.config>"
  exit
fi
grep -v "#" $1 > $2
echo 'done!'