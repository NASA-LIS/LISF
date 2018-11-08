      subroutine vapoprs(sfctmp,sfcprs,q2sat,q2)
      parameter ( T0 = 273.16 )
c
c  estimate the specific humidity, q2 in kg/kg, (Noah refer as the vapor 
c  mixing ratio) from air temperature, sfctmp in K deg., and pressure, 
c  sfcprs in Pa, q2sat is the saturation mixing ratio 
c  First, Popov's empirical relationship (adjusted to SnowMIP2 data) is used
c  to estimate vapor pressure:
c       e2 = 446.02*EXP(0.0579*sfctmp), in Pa
c  Then, it's converted into mixin ratio:
c       q2 = 0.622*e2/(sfcprs-(1-0.622)*e2), in kg/kg
c
c  calculate saturation mixing ratio
cfews      es=esat(sfctmp)
cfews      q2sat=0.622*es/(SFCPRS-(1.-0.622)*es)
      
c  then estimate actual mixing ratio
      e2 = 446.02*EXP(0.0579*(sfctmp-T0))
      q2 = 0.622*e2/(SFCPRS-(1.-0.622)*e2)
      if(q2 .lt. 0.1E-5) q2=0.1E-5
      if(q2 .ge. q2sat) q2=q2sat*0.99

      return
      end
                   
