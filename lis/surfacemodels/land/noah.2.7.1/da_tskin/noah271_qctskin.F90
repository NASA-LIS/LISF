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
! !ROUTINE: noah271_qctskin
! \label{noah271_qctskin}
!
! !REVISION HISTORY:
! 13Nov2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_qctskin(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify, LIS_logunit
  use noah271_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  QC's the tskin related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
!  type(ESMF_Field)       :: avgsurftField
  type(ESMF_Field)       :: soiltField
  integer                :: t
  integer                :: status
  real                   :: avgsurft(LIS_rc%ntiles(n))
  real, allocatable          :: soilt(:)
!  real                   :: tskinmax1
!  real                   :: tskinmin1
 
!  call ESMF_StateGet(LSM_State,"Skin Temperature",avgsurftField,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGet(avgsurftField,localDE=0,farrayPtr=avgsurft,rc=status)
!  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"Soil Temperature",soiltField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(soiltField,localDE=0,farrayPtr=soilt,rc=status)
  call LIS_verify(status)

!  call ESMF_AttributeGet(avgsurftField,"Max Value",tskinmax1,rc=status)
!  call ESMF_AttributeGet(avgsurftField,"Min Value",tskinmin1,rc=status)
!  call LIS_verify(status)
  
#if 0 
  do t=1,LIS_rc%ntiles(n)
     avgsurft(t) = noah271_struc(n)%noah(t)%t1
  enddo

  call noah271_calc_tskin(n,soilt,avgsurft)
#endif

  do t=1,LIS_rc%ntiles(n)
     
!     if(t.eq.1) then 
!        write(LIS_logunit,*) 'qc ',t,soilt(t), avgsurft(t)
!     endif
!     if(avgsurft(t).gt.tskinmax1) avgsurft(t) = tskinmax1
!     if(avgsurft(t).lt.tskinmin1) avgsurft(t) = tskinmin1
     if(soilt(t).gt.400) soilt(t) = noah271_struc(n)%noah(t)%stc(1)
     if(soilt(t).lt.200) soilt(t) = noah271_struc(n)%noah(t)%stc(1)
!     noah271_struc(n)%noah(t)%t1 = avgsurft(t)

  enddo

end subroutine noah271_qctskin

