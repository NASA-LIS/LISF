      subroutine soltmxmn(month,day,dtsec,dtlocal,tmax,tmin,rlat,
     +                          rsmin,rsmax,rgl,xlai,sfsl,solardt)
c
c--------------------------------------------------------------------------------
c  Estimation incomming short wave solar radiation from Bristow & Campbell (1984)
c  relationship for daily total radiation: Rg=So * Ktrs, in cal/(cm2*day)
c  where Ktrs is the atmosphere transmissivity, dimensionless:
c     Ktrs = A(1-EXP(-B*(tmax-tmin)**C)); A, B, & C parameters: A=0.75, B=0.004
c     (summer, May-October) or =0.01 (winter, Noevmber-April), and C=2.4
c  So is the clear sky solar irradiance on the horisontal plane (Dingman):
c     So = 2*RISC *Eo *(cos(lat)cos(dlt)sin(rot*Trs)/rot+sin(lat)sin(dlt)Trs)
c     RISC is the solar constant = 117.54 cal/(cm2*hr), Eo is the eccentricity 
c     correction, lat is the latitude, dlt is the sun declinantion, rot is the 
c     angular velocity of the earth's rotation =0.2618 radian/hr, and Trs is 
c     a half of the daylight time in hrs
c  Estmated daily radiation Rg then converted into time interval instantenious 
c  radiation values using sun angle at a specific time interval
c   solardt() - time interval average radiation, W/m2
c   rgx() - time interval factors to distribute hourly average radiation   
c   tnoon - local noon time
c   dtlocal - difference between local and Z time
c   sfsl - is daily sum of canopy resistance factors to solar radiation which
c          will be used to distribute daily ET demand 
c--------------------------------------------------------------------------------
c
      parameter (A = 0.75, BSM = 0.004, BWT = 0.01, C = 2.4)
      parameter (ROTV = 0.2618, RISC = 117.54)
      parameter (TNOON = 12.0)
      
      real solardt(144),rgx(144)
      integer dtlocal
      
c  calculate atmospheric transmissivity, radj
      if(month .ge. 5 .and. month .le. 10) then
       B=BSM
      else
       B=BWT
      endif
      radj=A*(1.0-EXP(-B*(tmax-tmin)**C))
      
c get sun rise/set      
      call sunrset(month,day,rlat,rise,set,sin_term,
     +                                     cos_term,reldst)

c  getting hourly average daily solar radiation, cal/(cm2*hr)
      trs=TNOON-rise
      trmis=2*(cos_term*sin(ROTV*trs)/ROTV+sin_term*trs)
      solar=RISC*reldst*trmis*radj

c loop throuugh day by each time interval       
      ntnoon=TNOON*3600/dtsec+0.001
      ndhr=24*3600/dtsec+0.001
      srg=0.0
      ind=0
      rgx(ndhr)=0.0
      do i=1,ntnoon
       time=i*dtsec/3600.
       stime=time-TNOON
       if(time .le. rise) then
        rgx(i)=0.0
        rgx(ndhr-i)=0.0
       else
        rgx(i)=sin_term+cos_term*cos(ROTV*stime)
        rgx(ndhr-i)=rgx(i)
       endif
       srg=srg+2*rgx(i)
      enddo
      srg=srg-rgx(ntnoon)

      ist=24-dtlocal
      sfsl=0.0 
      do i=1,ndhr
cc  time interval radiation in cal/(cm2*dt)
       solarx=solar*rgx(i)/srg
cc  time interval radiation in W/m2
       if(i .le. ist) then
        k=i+dtlocal
       else
        k=i-ist
       endif
       solardt(k)=41840*solarx/dtsec
cc  calculate daily total radiation term of canopy resistance
       fx=0.55*2.0*solardt(k)/(rgl*xlai)
       sfsl=sfsl+(rsmin/rsmax+fx)/(1.0+fx)
      enddo

      return
      end
        
       
      
      
