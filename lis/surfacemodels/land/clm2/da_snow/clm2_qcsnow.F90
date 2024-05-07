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
! !ROUTINE: clm2_qcsnow
! \label{clm2_qcsnow}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 06Oct2008: Gabrielle De Lannoy: limit the updated state to realistic values
! 04Nov2008: Gabrielle De Lannoy: now only for 2 state vars (total snow water, depth)
!
! !INTERFACE:
subroutine clm2_qcsnow(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clm2_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation - limit outcome to min/max in attribute-file
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: h2osnoField
  type(ESMF_Field)       :: snowdpField

  integer                :: t
  integer                :: status

  real, pointer          :: h2osno(:)  ! clm2.0 total snow water (kg m^-2)
  real, pointer          :: snowdp(:)  ! clm2.0 total snow depth

  real                   :: h2osnomin,h2osnomax
  real                   :: snowdpmin,snowdpmax

! -----------------------------------------------------------------------

  call ESMF_StateGet(LSM_State,"Total Snow Water",h2osnoField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Total Snow Depth",snowdpField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(h2osnoField,localDE=0,farrayPtr=h2osno,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snowdpField,localDE=0,farrayPtr=snowdp,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(h2osnoField,"Max Value",h2osnomax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(h2osnoField,"Min Value",h2osnomin,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(snowdpField,"Max Value",snowdpmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(snowdpField,"Min Value",snowdpmin,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     if(h2osno(t).lt.h2osnomin) then 
        h2osno(t) = h2osnomin
     endif
     if(h2osno(t).gt.h2osnomax) then
        h2osno(t) = h2osnomax
     endif
     if(snowdp(t).lt.snowdpmin) then
        snowdp(t) = snowdpmin
     endif
     if(snowdp(t).gt.snowdpmax) then
        snowdp(t) = snowdpmax
     endif
    !Additional qc's to ensure consistency..

     if(h2osno(t) .le. 0.0) snowdp(t) = 0.0
     if(snowdp(t) .le. 0.0) h2osno(t) = 0.0
         
  enddo

end subroutine clm2_qcsnow

