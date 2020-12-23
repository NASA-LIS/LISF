DSET ^OUTPUT/SOS_current_2010.d01.bil
*options sequential
*options big_endian
options yrev
TITLE WRSI output
UNDEF -9999.0
*LIS-WRSI:
XDEF  294  LINEAR    22.05  0.1
YDEF  348  LINEAR   -11.65  0.1
ZDEF 1 LINEAR 1 1
TDEF 1 LINEAR 06:00Z10sep2009 10dy
VARS 1
SOS     1  -1,40,4 ** Start-of-season [in dekads]
ENDVARS
