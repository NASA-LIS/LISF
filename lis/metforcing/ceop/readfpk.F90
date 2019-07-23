!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readfpk
! \label{readfpk}
! 
! !REVISION HISTORY:
!
!  08 Dec 2002: Sujay Kumar; Initial Specification 
!
! !INTERFACE:      
subroutine readfpk(n,order,ftn,tair,qair,swdown,lwdown,u,v,psurf,pcp)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_timeMgrMod, only : LIS_tick
  use LIS_logMod, only : LIS_logunit, LIS_endrun
  use ceop_forcingMod, only :ceop_struc
  
  implicit none
! !ARGUMENTS:
  integer,  intent (IN)     :: n 
  integer , intent (IN)     :: ftn
  integer , intent (IN)     :: order
  real    , intent (INOUT)  :: psurf(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: tair(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: qair(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: u(ceop_struc(n)%nstns),v(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: pcp(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: swdown(ceop_struc(n)%nstns)
  real    , intent (INOUT)  :: lwdown(ceop_struc(n)%nstns)
!
! !DESCRIPTION: 
!
!  Reads the CEOP station data for Fort Peck, MT. The format of 
!  the data is ASCII. 
!
!  The arguments are: 
! 
!  \begin{description}
!  \item[n]    index of the nest \newline
!  \item[order] flag indicating which data to 
!               be read (order=1, read the previous 
!               instance, order=2, read the next instance)
!  \item[ftn] unit number used to open the file \newline
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

  integer :: yr,mo,da,hr,mn
  real :: lat,lon,elev
  real :: dewp,rh,w,wdir,snowd,qle,qh,swout,lwout,netrad,skint
  real :: co2,par1,par2
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  real :: gmt1, ts1               
  real*8 :: time1
! there is only one station..
  if(order.eq.1) then 
     do 
        read(ftn,333) yr,mo,da,hr,mn,&
             yr,mo,da,hr,mn,lat,lon,elev,psurf,&
             tair(1),dewp,rh,qair(1),w,wdir,u(1),v(1),pcp(1),&
             snowd,qle,qh,swdown(1),swout,lwdown(1),lwout,&
             netrad,skint,co2,par1,par2
        pcp(1) = pcp(1)/3600
        if((LIS_rc%yr.eq.yr).and. &
             (LIS_rc%mo.eq.mo).and.&
             (LIS_rc%da.eq.da).and.&
             (LIS_rc%hr.eq.hr).and.&
             (LIS_rc%mn.eq.mn)) then
           
           exit
        endif
     enddo
  else
     read(ftn,333) yr,mo,da,hr,mn,&
          yr,mo,da,hr,mn,lat,lon,elev,psurf,&
          tair(1),dewp,rh,qair(1),w,wdir,u(1),v(1),pcp(1),&
          snowd,qle,qh,swdown(1),swout,lwdown(1),lwout,&
          netrad,skint,co2,par1,par2
     pcp(1) = pcp(1)/3600
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
  endif

333  format(I4,4(1X,I2),1X,I4,4(1X,I2),43X,F11.5,F12.5,12F8.2,&
         11F9.2)    

end subroutine readfpk
