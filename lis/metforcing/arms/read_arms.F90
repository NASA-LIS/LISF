!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_arms
! \label{read_arms}
! 
! !REVISION HISTORY:
!
!  28 Jun 2009: Sujay Kumar; Initial Specification 
!
! !INTERFACE:      
subroutine read_arms(n, currTime, findex,order) 
! !USES:
  use ESMF
  use LIS_logMod, only      : LIS_logunit, LIS_verify
  use LIS_timeMgrMod, only  : LIS_calendar, LIS_time2date
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_coreMod, only     : LIS_rc,LIS_domain
  use LIS_logMod,  only     : LIS_getNextUnitNumber, LIS_releaseUnitNumber, &
       LIS_endrun
  use LIS_constantsMod, only : LIS_CONST_TKFRZ, LIS_CONST_STEBOL
  use arms_forcingMod, only : arms_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: n
  real*8,  intent(in)           :: currTime
  integer, intent(in)           :: findex
  integer, intent(in)           :: order
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  the Walnut Gulch station data (ASCII), transforms into LIS forcing 
!  parameters and interpolates to the LIS domain. Please note that the
!  routine applies the met forcing (except precip) uniformly over the 
!  entire domain. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[currTime]
!    time at which forcing data should be read
!  \item[findex]
!   index of the supplemental forcing
!  \item[metdata]
!   output forcing array that stores the spatially interpolated
!   forcing data. 
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[normalize\_stnwts](\ref{normalize_stnwts}) \newline
!    renormalizes the station weights accounting for
!    missing data
!  \item[interp\_stndata](\ref{interp_stndata}) \newline
!    spatially interpolates the station data onto the 
!    LIS grid.
!  \end{description}
!EOP  
  integer :: yr, mo, da, doy, hr, mn
  integer :: i,j,c,r
  real    :: gmt
  integer :: yr1, mo1, da1, hr1, mn1
  integer :: ios, status
  integer :: offset
  integer :: ftn
  real    :: tmp, rh, radswg, lwgdwn, u, dummy,ps
  logical :: file_exists, readflag
  type(ESMF_Time) :: dataTime, stTime
  real    :: varfield(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer :: stnid(arms_struc(n)%nstns)
  real    :: pcp(arms_struc(n)%nstns)
  real    :: E_a, Esat, emissi_a

  if(arms_struc(n)%startRead) then 
     inquire(file=trim(arms_struc(n)%armsfile), exist=file_exists)

     if(.not.file_exists) then 
        write(LIS_logunit,*) 'ARMS file ',trim(arms_struc(n)%armsfile), ' not found'
        write(LIS_logunit,*) 'Program stopping ....'
        call LIS_endrun()
     endif

     ftn = LIS_getNextUnitNumber()
     open(ftn,file=trim(arms_struc(n)%armsfile),form='formatted',status='old')     
     readflag = .true. 
     
     i = 1
     do while (readflag) 
        read(ftn,fmt=300, iostat=ios) yr, doy, hr, tmp, rh, u, ps, radswg, lwgdwn, dummy
        if(ios.ne.0) then 
           readflag = .false. 
        endif
        
        if(readflag) then 
           call doytoCalday(yr,doy, mo, da)
           call ESMF_TimeSet(dataTime, yy=yr,mm=mo,dd=da,h =hr,calendar=LIS_calendar,&
                rc=status)
           call LIS_verify(status,'error in ESMF_TimeSet in read_arms')
           
           offset =nint((dataTime - arms_struc(n)%startDate)/arms_struc(n)%timeStep)+1

           arms_struc(n)%t_a(offset) = tmp+LIS_CONST_TKFRZ
           call RHtoq(tmp,ps, rh,Esat,arms_struc(n)%q_a(offset))
           arms_struc(n)%swd(offset) = radswg
!------------------------------------------------------------------------------------
!        equation from Brutsaert (1975; WRR)
!------------------------------------------------------------------------------------
           E_a = (rh/100.0)*(Esat/100.0)
           emissi_a = 1.24 * (E_a/arms_struc(n)%t_a(offset))**(1.0/7.0) 
           
           arms_struc(n)%lwd(offset) = emissi_a*LIS_CONST_STEBOL*arms_struc(n)%t_a(offset)**4
           arms_struc(n)%uwind(offset) = u
           arms_struc(n)%psurf(offset) = ps*100.0
           
        endif

     enddo
     call LIS_releaseUnitNumber(ftn)
!------------------------------------------------------------------------------------
!        Read the precip data
!------------------------------------------------------------------------------------
     ftn = LIS_getNextUnitNumber()
     write(LIS_logunit,*) 'Reading ARMS pcp file ',trim(arms_struc(n)%armspcpfile)
     open(ftn,file=trim(arms_struc(n)%armspcpfile),form='formatted',status='old')     
     readflag = .true. 
     read(ftn,332) (stnid(j),j=1,arms_struc(n)%nstns)
     do while(readflag) 
        read(ftn,fmt=333,iostat=ios) mo, da, yr, hr, mn, (pcp(j),j=1,arms_struc(n)%nstns)

        call ESMF_TimeSet(dataTime, yy=yr,mm=mo,dd=da,h =hr,m=mn, calendar=LIS_calendar,&
             rc=status)
        call LIS_verify(status,'error in ESMF_TimeSet in read_arms')
        
        offset =nint((dataTime - arms_struc(n)%startDate)/arms_struc(n)%timeStep)+1
        
        arms_struc(n)%precip(offset, :) = pcp(:)
        if(ios.ne.0) then 
           readflag = .false. 
        endif
        
     enddo

332 format(16X,85I7)
333 format(I2,1X,I2,1X,I4,1X,I2,1X,I2,85F7.2)     

     call LIS_releaseUnitNumber(ftn)

     arms_struc(n)%startRead = .false. 
  endif

  call LIS_time2date(currTime, doy,gmt, yr1, mo1,da1, hr1,mn1)

  call ESMF_TimeSet(stTime, yy=yr1, mm=mo1,dd=da1, h=hr1, m=mn1, &
       calendar=LIS_calendar, rc=status)
  call LIS_verify(status, 'error in stTime set ')

  offset = nint((stTime - arms_struc(n)%startDate)/arms_struc(n)%timeStep)+1

  call normalize_stnwts( arms_struc(n)%precip(offset,:),arms_struc(n)%nstns, &
       LIS_rc%lnc(n)*LIS_rc%lnr(n),-9999.0, arms_struc(n)%stnwt)  
  call interp_stndata(arms_struc(n)%stnwt, -9999.0,  arms_struc(n)%precip(offset,:), &
       varfield, LIS_rc%lnc(n)*LIS_rc%lnr(n), arms_struc(n)%nstns)

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(varfield(c+(r-1)*LIS_rc%lnc(n)).lt.0) then 
              print*, 'problem',varfield(c+(r-1)*LIS_rc%lnc(n))
              print*, 'pcp ', arms_struc(n)%precip(offset,:)
              stop
           endif
           if(varfield(c+(r-1)*LIS_rc%lnc(n)).lt.1E-12) then 
              varfield(c+(r-1)*LIS_rc%lnc(n)) = 0.0 !eliminate floating pt noise
           endif
        endif
     enddo
  enddo

!  open(100,file='precip.bin',form='unformatted')
!  write(100) varfield
!  close(100)
!  stop

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              arms_struc(n)%metdata1(1,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%T_a(offset)
              arms_struc(n)%metdata1(2,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Q_a(offset)
              arms_struc(n)%metdata1(3,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Swd(offset)
              arms_struc(n)%metdata1(4,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Lwd(offset)
              arms_struc(n)%metdata1(5,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Uwind(offset)
              arms_struc(n)%metdata1(6,LIS_domain(n)%gindex(c,r)) = 0.0
              arms_struc(n)%metdata1(7,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%psurf(offset)
              arms_struc(n)%metdata1(8,LIS_domain(n)%gindex(c,r)) = &
                   varfield(c+(r-1)*LIS_rc%lnc(n))/3600.0
           elseif(order.eq.2) then 
              arms_struc(n)%metdata2(1,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%T_a(offset)
              arms_struc(n)%metdata2(2,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Q_a(offset)
              arms_struc(n)%metdata2(3,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Swd(offset)
              arms_struc(n)%metdata2(4,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Lwd(offset)
              arms_struc(n)%metdata2(5,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%Uwind(offset)
              arms_struc(n)%metdata2(6,LIS_domain(n)%gindex(c,r)) = 0.0
              arms_struc(n)%metdata2(7,LIS_domain(n)%gindex(c,r)) =  &
                   arms_struc(n)%psurf(offset)
              arms_struc(n)%metdata2(8,LIS_domain(n)%gindex(c,r)) = &
                   varfield(c+(r-1)*LIS_rc%lnc(n))/3600.0
           endif

        endif
     end do
  enddo


300 format(I4,I5,I6,7F8.2)

end subroutine read_arms

!BOP
! 
! !ROUTINE: doytoCalday
! 
! !INTERFACE: 
subroutine doytoCalday(yr,doy, mo, da)

  implicit none
!
! !ARGUMENTS: 
  
  integer, intent(in) :: yr
  integer, intent(in) :: doy
  integer, intent(out) :: mo
  integer, intent(out) :: da
!
! !DESCRIPTION: 
!  This subroutine converts the day of the year to the calendar
!  month and day
! 
!EOP
  logical :: leapyear
  integer :: i 
  integer :: days1(12), days2(12), days

  data days1 /0, 31,59,90,120,151,181,212,243,273,304,334/
  data days2 /0, 31,60,91,121,152,182,213,244,274,305,335/

  if(mod(yr,4).eq.0.and.mod(yr,100).ne.0) then !leap year
     leapyear = .true.
  else
     leapyear = .false. 
  endif

  do i = 1, 12
     if(leapyear) then
        days = days2(i)
     else
        days = days1(i)
     endif
     if(days.lt.doy) then 
        mo = i
     endif
  enddo
  
  if(leapyear) then 
     da = doy - days2(mo)
  else
     da = doy - days1(mo)
  endif
  
end subroutine doytoCalday

!BOP
! !ROUTINE: RHtoq
! 
! !INTERFACE:
subroutine RHtoq(TC,ph,RHp,Es,q) 
!
! !Revision History:
! 03/06/1998  Pablo Grunmann  Initial specification for Noah LSM as 'qdatap.f'
! 06/24/1998  Pablo Grunmann  Eliminated subroutine SVP
!  04 Dec 04  Matthew Garcia  Edit/update for incorporation in LISv3.2
!  24 May 05  Matthew Garcia  Edit/update for LISv4.1 and calling arguments
!
  implicit none
! !INPUT ARGUMENTS: 
! In:    TC     Temperature (C)
!        ph     Pressure (hPa)
!        RHp    Relative humidity (%)
  real, intent(in) :: TC,ph,RHp
! !OUTPUT ARGUMENTS:
! Out:   q      Specific humidity (Kg/Kg)
!        Es     Saturation vapor pressure (hPa)
  real, intent(out) :: q,Es
! 
! !DESCRIPTION:  OBTAIN SPECIFIC HUMIDITY (q) FROM RELATIVE HUMIDITY 
!           AND GIVEN PRESSURE AND TEMPERATURE.
!
!           FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989: 'A 
!           SHORT COURSE IN CLOUD PHYSICS', PERGAMON PRESS, 3rd ED.
! 
!EOP
!
! Local Variables
  real :: T,p,RH
  real :: EP,LW,qsat
!
! Parameters
  real, parameter :: eps = 0.622 ! (Water)/(dry air) molecular mass ratio, epsilon
  real, parameter :: LVH2O = 2.501000E+6 ! latent heat of vaporization for water
  real, parameter :: CPV = 1870.0 ! specific heat of water vapor
  real, parameter :: RV = 461.5 ! universal gas constant for water vapor
  real, parameter :: CW = 4187.0 ! specific heat of water liquid
  real, parameter :: Es_0 = 611.2 ! saturation vapor pressure (Pa) at water triple point
  real, parameter :: T_0 = 273.16 ! temperature at water triple point
!
! convert temperature to the value in Kelvin
  T = TC + 273.15
!
! convert pressure to the value in Pascals
  p = ph * 100.0
!
! convert RH to the fractional value
  RH = RHp / 100.0
!
! calculate the saturation vapor pressure (Es, in Pascal) 
!  at temperature T using Clausius-Clapeyron equation
  LW = LVH2O - (CW - CPV) * (T - T_0)
  Es = Es_0 * EXP(LW * (1/T_0 - 1/T) / RV)  
!
! convert RH to q (Rogers, pg. 17)
  EP = (p * Es * RH) / (p - Es * (1. - RH))
  q = eps * EP / (p - (1. - eps) * EP)
!
! calculate saturation mixing ratio (qsat)
  qsat = 0.622 * Es / (p - (1. - 0.622) * Es)
!
! check values against minimum and qsat     
  if (q.lt.(0.1E-5)) q = 0.1E-5
  if (q.ge.qsat) q = qsat * 0.99
!
  return
end subroutine RHtoq
