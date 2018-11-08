      subroutine smc2sac(dsm,dsmd,et,xztwc,xzfwc,xzfwc1,xztwh,xzfwh,
     +            xzfwh1,edmnd,e1,e2,red,xztwm,xzfwm,xzfwm1)

c -------------------------------------------------------------
c  Recalculate soil layer moisture change from evaporation 
c  into SAC upper or lower storages. Redistributes changes
c  to tension storages first, and left to free storage
c dsm is SMC water change to distribute to SAC zone
c dsmd is the zone water excess if any 
c -------------------------------------------------------------       
      e1=0.0
      e2=0.0
      red=0.0
      dx=0.0
      dx1=0.0
      dsmd=0.0
      dsmx=dsm

      if(dsmx .gt. 0.0) then
c add water to sac-sma zone
       xztwc=xztwc+dsmx
       if(xztwc .gt. xztwm) then
        dx=xztwc-xztwm
        xztwc=xztwm
        if(xzfwc1 .eq. -1.) then
c upper zone free water
         xzfwc=xzfwc+dx
         if(xzfwc .gt. xzfwm) then
          dsmd=xzfwc-xzfwm
          xzfwc=xzfwm
         endif
         xzfwh=xzfwh+dx-dsmd
         xztwh=xztwh+dsmx-dx
        else
c lower zone free water
         xzfwc1=xzfwc1+dx
         if(xzfwc1 .gt. xzfwm1) then
          dx1=xzfwc1-xzfwm1
          xzfwc1=xzfwm1
          xzfwc=xzfwc+dx1
          if(xzfwc .gt. xzfwm) then
           dsmd=xzfwc-xzfwm
           xzfwc=xzfwm
          endif
         endif
         xzfwh1=xzfwh1+dx-dx1
         xzfwh=xzfwh+dx1-dsmd  
         xztwh=xztwh+dsmx-dx
        endif 
       else
        xztwh=xztwh+dsmx
       endif
      else       
c subtract water from sac-sma zone
       xztwh=xztwh+dsmx
       if(xztwh .ge. 0.0) then
        xztwc=xztwc+dsmx
        e1=et
        red=edmnd-e1
       else 
        alloc=dsmx-xztwh
        dsmx=xztwh
        xztwc=xztwc+alloc
        xztwh=0.0
        xzfwh=xzfwh+dsmx
        e1=-alloc
        red=edmnd-e1
        if(xzfwh .ge. 0.0) then
         xzfwc=xzfwc+dsmx
         e2=edmnd-e1
         red=0.0
        else
         alloc=dsmx-xzfwh
         dsmx=xzfwh
         xzfwc=xzfwc+alloc
         xzfwh=0.0
         e2=-alloc
         red=red-e2
         if(xzfwh1 .ne. -1) then
          xzfwh1=xzfwh1+dsmx
          if(xzfwh1 .ge. 0.0) then
           xzfwc1=xzfwc1+dsmx
          else
           alloc=dsmx-xzfwh1
           dsmx= xzfwh1
           xzfwc1=xzfwc1+alloc
           xzfwh1=0.0
           dsmd=dsmx
          endif
         endif  
        endif 
       endif
      endif

      return
      end 

