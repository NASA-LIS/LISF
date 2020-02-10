!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: HYMAP2_updateWL
!  \label{HYMAP2_updateWL}
!
! !REVISION HISTORY:
!  07 Nov 2019: Sujay Kumar, Initial Specification
!
! !INTERFACE:
subroutine HYMAP2_updateWL(n, Routing_State, Routing_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use HYMAP2_routingMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: Routing_State
  type(ESMF_State)       :: Routing_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: rivstoField
  type(ESMF_Field)       :: rivstoIncrField
  real, pointer          :: rivsto(:)
  real, pointer          :: rivstoIncr(:)
  integer                :: t,i,m
  integer                :: status

  call ESMF_StateGet(Routing_State,"River storage",rivstoField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: River storage failed in HYMAP2_updateWL")

  call ESMF_FieldGet(rivstoField,localDE=0,farrayPtr=rivsto,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: River storage failed in HYMAP2_updateWL")

  call ESMF_StateGet(Routing_Incr_State,"River storage",rivstoIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: River storage failed in HYMAP2_updateWL")

  call ESMF_FieldGet(rivstoIncrField,localDE=0,farrayPtr=rivstoIncr,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: River storage failed in HYMAP2_updateWL")

  do t=1,HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n)
     rivsto(t) = rivsto(t) + rivstoIncr(t)
  enddo
end subroutine HYMAP2_updateWL

