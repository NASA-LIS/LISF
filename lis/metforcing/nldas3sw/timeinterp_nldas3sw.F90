!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: timeinterp_nldas3sw
! \label{timeinterp_nldas3sw}
!
! !REVISION HISTORY:
! 27 Dec 2024: David Mocko, Initial Specification
!                           (derived from timeinterp_nldas20.F90)
!
! !INTERFACE:
subroutine timeinterp_nldas3sw(n,findex)
! !USES:
  use ESMF
  use LIS_FORC_AttributesMod
  use LIS_coreMod,        only  : LIS_rc,LIS_domain
  use LIS_metforcingMod,  only  : LIS_forc,LIS_FORC_Base_State
  use LIS_timeMgrMod,     only  : LIS_tick,LIS_time2date
  use LIS_logMod,         only  : LIS_logunit,LIS_verify,LIS_endrun
  use nldas3sw_forcingMod, only : nldas3sw_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Linear interpolation available for SWdown, but should not be
!  used in LIS.  Simply read in hourly 4-km CERES SWdown data,
!  and write hourly 1-km NLDAS-3 SWdown data in the lis.config.
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[LIS\_tick](\ref{LIS_tick}) \newline
!    advances or retracts time by the specified amount
!  \end{description}
!
!EOP
  real    :: wt1,wt2,gmt1,gmt2,tempbts
  integer :: t
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss

  integer          :: status
  type(ESMF_Field) :: swdField
  real, pointer    :: swd(:)

!________________________________________

  btime = nldas3sw_struc(n)%nldas3swtime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)

  tempbdoy = bdoy
  tempgmt1 = gmt1
  tempbyr = byr
  tempbmo = bmo
  tempbda = bda
  tempbhr = bhr
  if (tempbhr.eq.24) tempbhr = 0
  tempbmn = bmn
  tempbss = 0
  tempbts = 0
  call LIS_tick(newtime1,tempbdoy,tempgmt1,tempbyr,tempbmo,tempbda,    &
                tempbhr,tempbmn,tempbss,tempbts)

  btime = nldas3sw_struc(n)%nldas3swtime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  tempbdoy = bdoy
  tempgmt2 = gmt2
  tempbyr = byr
  tempbmo = bmo
  tempbda = bda
  tempbhr = bhr
  if (tempbhr.eq.24) tempbhr = 0
  tempbmn = bmn
  tempbss = 0
  tempbts = 0
  call LIS_tick(newtime2,tempbdoy,tempgmt2,tempbyr,tempbmo,tempbda,    &
                tempbhr,tempbmn,tempbss,tempbts)

!=== Interpolate Data in time
  wt1 = (nldas3sw_struc(n)%nldas3swtime2 - LIS_rc%time) /              &
        (nldas3sw_struc(n)%nldas3swtime2 - nldas3sw_struc(n)%nldas3swtime1)
  wt2 = 1.0 - wt1

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),                    &
       LIS_FORC_SWdown%varname(1),swdField,rc=status)
  call LIS_verify(status,                                              &
       "[ERR] Enable SWdown in the forcing variables list")

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  do t = 1,LIS_rc%ntiles(n)
     if ((nldas3sw_struc(n)%metdata1(1,t).ne.LIS_rc%udef).and.         &
             (nldas3sw_struc(n)%metdata2(1,t).ne.LIS_rc%udef)) then
        swd(t) = (nldas3sw_struc(n)%metdata1(1,t)*wt1) +               &
                 (nldas3sw_struc(n)%metdata2(1,t)*wt2)
     else
        swd(t) = 0.0
     endif
  enddo

end subroutine timeinterp_nldas3sw

