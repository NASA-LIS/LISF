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
! !ROUTINE: clsmf25_qcsnow
! \label{clsmf25_qcsnow}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine clsmf25_qcsnow(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use clsmf25_lsmMod
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: swe1Field, swe2Field, swe3Field
  integer                :: t
  integer                :: status
  real, pointer          :: swe1(:), swe2(:), swe3(:)
  real                   :: swe1max,swe1min, swe2max,swe2min, swe3max,swe3min
 

  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 1",swe1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 2",swe2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 3",swe3Field,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(swe1Field,"Max Value",swe1max,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(swe1Field,"Min Value",swe1min,rc=status)
  call LIS_verify(status)


  call ESMF_AttributeGet(swe2Field,"Max Value",swe2max,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(swe2Field,"Min Value",swe2min,rc=status)
  call LIS_verify(status)


  call ESMF_AttributeGet(swe3Field,"Max Value",swe3max,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(swe3Field,"Min Value",swe3min,rc=status)
  call LIS_verify(status)


  call ESMF_FieldGet(swe1Field,localDE=0,farrayPtr=swe1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe2Field,localDE=0,farrayPtr=swe2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe3Field,localDE=0,farrayPtr=swe3,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(swe1(t).lt.swe1min) then 
        swe1(t) = swe1min
     endif
     if(swe1(t).gt.swe1max) then 
        swe1(t) = swe1max
     endif
     if(swe2(t).lt.swe2min) then 
        swe2(t) = swe2min
     endif
     if(swe2(t).gt.swe2max) then 
        swe2(t) = swe2max
     endif
     if(swe3(t).lt.swe3min) then 
        swe3(t) = swe3min
     endif
     if(swe3(t).gt.swe3max) then 
        swe3(t) = swe3max
     endif
  enddo

end subroutine clsmf25_qcsnow

