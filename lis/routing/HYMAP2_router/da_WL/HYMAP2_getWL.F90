!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: HYMAP2_getWL
! \label{HYMAP2_getWL}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial specification
!
! !INTERFACE:
subroutine HYMAP2_getWL(n, Routing_State)

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
!  Returns the water level DA related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[Routing\_State] ESMF State container for Routing state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: rivstoField
  integer                :: t,i,m
  integer                :: status
  real, pointer          :: rivsto(:)
  character*100          :: lsm_state_objs(4)

  call ESMF_StateGet(Routing_State,"River storage",rivstoField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in HYMAP2_getWL')

  call ESMF_FieldGet(rivstoField,localDE=0,farrayPtr=rivsto,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for rivsto in HYMAP2_getWL')

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m
        rivsto(t) = HYMAP2_routing_struc(n)%rivsto(i,m)
     enddo
  enddo

end subroutine HYMAP2_getWL

