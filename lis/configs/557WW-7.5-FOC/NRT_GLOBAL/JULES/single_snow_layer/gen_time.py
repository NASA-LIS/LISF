#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

from datetime import datetime, timedelta
f = open("loobos_time.txt", "w")
delta = timedelta(seconds=1800)
t0 = datetime(1997, 01, 01, 0, 0, 0)

N = 365*24*2
for k in range(0,N):
  t = t0 + k * delta
  str = "%04d %02d %02d %02d %02d \n" %(t.year, t.month, t.day, t.hour, t.minute)
  f.write(str)
f.close() 

