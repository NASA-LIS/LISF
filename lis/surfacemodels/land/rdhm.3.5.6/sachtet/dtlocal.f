      function offsetlocal(xhrap,yhrap)

c--------------------------------------------------------------------
c  Estimation of offset of time comparing to Z-time      
c  IDAY - day of year. Not much difference for different days. Assume 
c         constant value =361 when offset close to annual average
c  IHR  - hour, offset practically constant for all hours. 
c         Assume constant = 12
c  ITZONE - zone time in hours; uses Greenvich meridian zone time 
c  xhrap, yhrap - HRAP coordinates
c--------------------------------------------------------------------
 
      parameter (ITZONE = 0, IHR = 12, IDAY = 361)
      
c  Convert HRAP to LatLon
      call sbllgd(rlon,rlat,1,xhrap,yhrap,0,istat)
      if(istat .eq. -1) then
       write(*,*) 'ERROR:',wron hrap to LatLon transform,xhrap,yhrap
       stop
      endif

c  Estimate offset in hours
      gama=2*3.14*(IDAY-1+(IHR-12)/24.0)/365.0
      egtime=229.18*(0.000075+0.001868*cos(gama)-0.032077*sin(gama)-
     +               0.014615*cos(2*gama)-0.040849*sin(2*gama))
      dtime=(egtime-4*rlon+60*ITZONE)
      offsetlocal=-dtime/60.0

      return
      end
