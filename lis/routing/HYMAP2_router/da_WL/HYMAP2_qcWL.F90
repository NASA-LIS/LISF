!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: HYMAP2_qcWL
! \label{HYMAP2_qcWL}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP2_qcWL(n, Routing_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use HYMAP2_routingMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
!
! !DESCRIPTION:
!
!  QCs the water level related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: rivstoField
  integer                :: t
  integer                :: status
  real, pointer          :: rivsto(:)
  real                   :: rivstomax
  real                   :: rivstomin

#if 0 
  call ESMF_StateGet(Routing_State,"River storage",rivstoField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for River storage failed in HYMAP2_qcWL")
 
  call ESMF_FieldGet(rivstoField,localDE=0,farrayPtr=rivsto,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for River storage failed in HYMAP2_qcWL")

  call ESMF_AttributeGet(rivstoField,"Max Value",rivstomax,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in HYMAP2_qcWL")

  call ESMF_AttributeGet(rivstoField,"Min Value",rivstomin,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in HYMAP2_qcWL")

  do t=1,HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n)
     if(rivsto(t).gt.rivstomax) rivsto(t) = rivstomax
     if(rivsto(t).lt.rivstomin) rivsto(t) = rivstomin
     print*, 'here in h2_qcWL'
  enddo

#endif
end subroutine HYMAP2_qcWL

