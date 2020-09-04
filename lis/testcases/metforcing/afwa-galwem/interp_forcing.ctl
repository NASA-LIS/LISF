dset ^20160411_CY.%h2_FH.000_DF.GR2_lis.bin
options little_endian
options template
undef -9999
ydef 600 linear   -59.875 0.25
xdef 1440 linear -179.875 0.25
tdef 4 linear 00Z11apr2016 6hr
zdef 13 levels 100000 97500 95000 92500 90000 85000 80000 75000 70000 65000 60000 55000 50000
vars 5
sp       1  99 surface pressure
t       13  99 temperature
gh      13  99 geopotenital height
r       13  99 relative humidity
wspd     1  99 wind speed
endvars
