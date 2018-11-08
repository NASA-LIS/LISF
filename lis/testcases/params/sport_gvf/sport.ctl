*DSET ^gvf_SPORT_3KM.%y4%m2%d2.1gd4r
DSET ^gvf_SPORT_3KM.%y4%m2%ch.1gd4r
OPTIONS template
options big_endian
UNDEF -9999
XDEF 12000 LINEAR   -179.985  0.03
YDEF  5000 LINEAR    -59.985  0.03
ZDEF 1 LINEAR 1 1
TDEF 31 LINEAR 00Z01may2013 1dy
* use chsub when comparing raw SPoRT GVF data to LIS using realtime mode
* because in realtime mode, LIS reads GVF data two hours behind; e.g.,
* at hour 3, LIS reads GVF stamped at hour 1.
chsub 3 3 01
chsub 4 4 02
chsub 5 5 03
chsub 6 6 04
chsub 7 7 05
chsub 8 8 06
chsub 9 9 07
chsub 10 10 08
chsub 11 11 09
chsub 12 12 10
chsub 13 13 11
chsub 14 14 12
chsub 15 15 13
chsub 16 16 14
chsub 17 17 15
chsub 18 18 16
chsub 19 19 17
chsub 20 20 18
chsub 21 21 19
chsub 22 22 20
chsub 23 23 21
chsub 24 24 22
chsub 25 25 23
chsub 26 26 24
chsub 27 27 25
chsub 28 28 26
chsub 29 29 27
chsub 30 30 28
chsub 31 31 29
VARS 1
gvf      1  99  ** description 
ENDVARS
