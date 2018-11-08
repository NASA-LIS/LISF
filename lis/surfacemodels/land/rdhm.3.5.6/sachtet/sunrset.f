      subroutine sunrset(month,day,rlat,rise,set,sin_term,
     +                                          cos_term,reldst)
c
c -----------------------------------------------------------------
c  subroutine estimates a local time of sunrise and sunset
c  All equations from Dingman book, 2002
c  rlat - latitude
c  tnoon - local noon time
c  rise - local sun rise time
c  set - local sun set time
c  sin_term - sinus term in sun angle calculation
c  cos_term - cosine term in sun angle calculation
c  decline - declination of the sun
c  reldst - sun-earth eccentricity correction factor
c  Assumed that noon time = 12.0 for CONUS  
c -----------------------------------------------------------------       
c
      parameter (TNOON = 12.0)
      parameter (ROTV = 0.2618)
      integer ndays(12)/0,31,59,90,120,151,181,212,243,273,304,334/
      
      rlatr=rlat*3.14/180.
      nday=ndays(month)+day
      dangl=2*3.14*(nday-1)/365.
      decline=(0.006918-0.399912*cos(dangl)+0.070257*sin(dangl)-
     *         0.006758*cos(2*dangl)+0.000907*sin(2*dangl)-
     *         0.002697*cos(3*dangl)+0.00148*sin(3*dangl))
      reldst=1.00011+0.034221*cos(dangl)+0.00128*sin(dangl)+
     +           0.000719*cos(2*dangl)+0.000077*sin(2*dangl)
      z=acos(-tan(decline)*tan(rlatr))
      dts=z/ROTV
      rise=TNOON-dts
      set=TNOON+dts
      sin_term=sin(rlatr)*sin(decline)
      cos_term=cos(rlatr)*cos(decline)

      return
      end
  
