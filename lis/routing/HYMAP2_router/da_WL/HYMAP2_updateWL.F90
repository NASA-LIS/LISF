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

  type(ESMF_Field)       :: sfcelevField
  type(ESMF_Field)       :: sfcelevIncrField
  real, pointer          :: sfcelev(:)
  real, pointer          :: sfcelevIncr(:)
  integer                :: t,i,m
  integer                :: status

  call ESMF_StateGet(Routing_State,"Surface elevation",sfcelevField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Surface elevation failed in HYMAP2_updateWL")

  call ESMF_FieldGet(sfcelevField,localDE=0,farrayPtr=sfcelev,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Surface elevation failed in HYMAP2_updateWL")

  call ESMF_StateGet(Routing_Incr_State,"Surface elevation",sfcelevIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Surface elevation failed in HYMAP2_updateWL")

  call ESMF_FieldGet(sfcelevIncrField,localDE=0,farrayPtr=sfcelevIncr,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Surface elevation failed in HYMAP2_updateWL")

  do t=1,HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n)
     sfcelev(t) = sfcelev(t) + sfcelevIncr(t)
!     if(t.ge.25001.and.t.le.25020) then 
!        print*,'upd ',t,sfcelev(t),sfcelevIncr(t)
!     endif
  enddo
end subroutine HYMAP2_updateWL

