dset ^forcing_%h2.bin
options little_endian
options template
undef -9999
ydef 1152 linear -89.921875 0.15625
xdef 1536 linear 0.117187 0.234378
tdef 4 linear 00Z11apr2016 6hr
zdef 13 levels 100000 97500 95000 92500 90000 85000 80000 75000 70000 65000 60000 55000 50000
vars 6
sp       1  99 surface pressure
t       13  99 temperature
gh      13  99 geopotenital height
r       13  99 relative humidity
u10      1  99 u-wind
v10      1  99 v-wind
endvars
