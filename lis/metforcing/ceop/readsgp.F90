!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readsgp
! \label{readsgp}
! 
! !REVISION HISTORY:
!
!  08 Dec 2002: Sujay Kumar; Initial Specification 
!
! !INTERFACE:      
subroutine readsgp(n,tair,qair,swdown,lwdown,u,v,psurf,pcp, &
     errorcode)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick, LIS_date2time
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  use ceop_forcingMod, only : ceop_struc
  
  implicit none
! !ARGUMENTS:
  integer,  intent (IN)     :: n 
  real    , intent (INOUT)  :: psurf(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: tair(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: qair(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: u(ceop_struc(n)%nstns),v(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: pcp(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: swdown(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: lwdown(ceop_struc(n)%nstns)
  integer , intent (INOUT)  :: errorcode
! !DESCRIPTION: 
!
!  Reads the CEOP station data for SGP. The format of 
!  the data is ASCII. 
!
!  The arguments are: 
! 
!  \begin{description}
!  \item[n]    index of the nest \newline
!  \item[order] flag indicating which data to 
!               be read (order=1, read the previous 
!               instance, order=2, read the next instance)
!  \item[tair] array for storing the 2m air temperature for all stations \newline
!  \item[qair] array for storing the 2m relative humidity for all stations \newline
!  \item[swdown] array for storing the downward shortwave radiation
!              for all stations \newline
!  \item[lwdown] array for storing the downward longwave radiation
!              for all stations \newline
!  \item[u] array for storing the u component of wind for all stations \newline
!  \item[v] array for storing the v component of wind for all stations \newline
!  \item[psurf] array for storing the surface pressure for all stations \newline
!  \item[pcp] array for storing the precipitation for all stations \newline
!  \end{description}
!EOP
  integer :: f,i
  integer :: yr,mo,da,hr,mn
  real :: lat,lon
  real :: elev(ceop_struc(n)%nstns)
  real :: dewp,rh,w,wdir,snowd,qle,qh,swout,lwout,netrad,skint
  real :: co2,par1,par2
  real*8  :: listime, ctime
  real    :: lisgmt, cgmt
  integer :: lisdoy, cdoy
  real*8 :: time1

  errorcode = -1
  open(222,file=ceop_struc(n)%ceopdir,status='old')

!  if(order.eq.1) then 
     do 
        do i=1,ceop_struc(n)%nstns
           read(222,333) yr,mo,da,hr,mn,&
                yr,mo,da,hr,mn,lat,lon,elev(i),psurf(i),&
                tair(i),dewp,rh,qair(i),w,wdir,u(i),v(i),pcp(i),&
                snowd,qle,qh,swdown(i),swout,lwdown(i),lwout,&
                netrad,skint,co2,par1,par2
!           if(tair(i).ne.-999.99) tair(i) = tair(i)+8.1*elev(i)/1000.0
           if(pcp(i).ne.-999.99) pcp(i) = pcp(i)/3600
!           write(LIS_logunit,*) yr, mo, da, hr, mn, tair(i),qair(i),swdown(i),lwdown(i)!,u(i),v(i),&
!                psurf(i),pcp(i)
           
           call LIS_date2time(ctime,cdoy,cgmt,yr,mo,da,hr,mn,0)
           call LIS_date2time(listime,lisdoy, lisgmt, LIS_rc%yr, LIS_rc%mo,&
                LIS_rc%da, LIS_rc%hr,LIS_rc%mn, LIS_rc%ss)

!           write(LIS_logunit,*) yr,LIS_rc%yr
!           write(LIS_logunit,*) mo,LIS_rc%mo
!           write(LIS_logunit,*) da,LIS_rc%da
!           write(LIS_logunit,*) hr,LIS_rc%hr
!           write(LIS_logunit,*) mn,LIS_rc%mn
!           stop
!           write(LIS_logunit,*) dewp,rh,qair(i),w,wdir,u(i),v(i),pcp(i)
!           write(LIS_logunit,*) snowd,qle,qh,swdown(i),swout,lwdown(i),lwout
!           write(LIS_logunit,*) netrad,skint,co2,par1,par2
!           stop
        enddo
        if(ctime.gt.listime) then
           exit
        endif

        if((LIS_rc%yr.eq.yr).and. &
                (LIS_rc%mo.eq.mo).and.&
                (LIS_rc%da.eq.da).and.&
                (LIS_rc%hr.eq.hr).and.&
                (LIS_rc%mn.eq.mn)) then 
           errorcode = 0 
           exit
        endif
     enddo
     close(222)
     if(errorcode.eq.0) &
          write(LIS_logunit,*) 'Successfully read the CEOP station data..'
#if 0 
  else
     do i=1,ceop_struc(n)%nstns
        read(222,333) yr,mo,da,hr,mn,&
             yr,mo,da,hr,mn,lat,lon,elev(i),psurf(i),&
             tair(i),dewp,rh,qair(i),w,wdir,u(i),v(i),pcp(i),&
             snowd,qle,qh,swdown(i),swout,lwdown(i),lwout,&
             netrad,skint,co2,par1,par2
        pcp = pcp/3600
        yr1 = LIS_rc%yr  !current time
        mo1 = LIS_rc%mo
        da1 = LIS_rc%da
        hr1 = LIS_rc%hr
        mn1 = LIS_rc%mn
        ss1 = LIS_rc%ss
        ts1 = 60*60
        call LIS_tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   
        if((yr1.ne.yr).or.&
             (mo1.ne.mo).or.&
             (da1.ne.da).or.&
             (hr1.ne.hr).or.&
             (mn1.ne.mn)) then 
           write(LIS_logunit,*) 'Time mismatch.., Stopping..'
           write(LIS_logunit,*) yr,yr1
           write(LIS_logunit,*) mo,mo1
           write(LIS_logunit,*) da,da1
           write(LIS_logunit,*) hr,hr1
           write(LIS_logunit,*) mn,mn1
           call LIS_endrun
        endif
!        if(tair(i).ne.-999.99) tair(i) = tair(i)+8.1*elev(i)/1000.0
     enddo
  endif
#endif

333  format(I4,4(1X,I2),1X,I4,4(1X,I2),43X,F11.5,F12.5,12F8.2,&
         11F9.2)    

end subroutine readsgp
