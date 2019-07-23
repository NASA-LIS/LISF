!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_coop
! \label{timeinterp_coop}
! 
!
! !REVISION HISTORY:
! 13 Jul 2011: Yuqiong Liu; Initial Specification
!
! !INTERFACE:
subroutine timeinterp_coop(n, findex)
! !USES:
  use ESMF
  use LIS_logMod,  only        : LIS_logunit
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
  use LIS_FORC_AttributesMod
  use LIS_logMod,         only : LIS_verify
  use coop_forcingMod,    only : coop_struc
    
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

! !DESCRIPTION:
! 
! Temporally interpolates the SCAN forcing data to the 
! model timestep. All variables except precipitation is 
! linearly interpolated. 
! 
!EOP
  
  integer :: t
  integer :: index1
  integer            :: status
  type(ESMF_Field)   :: pcpField
  real,pointer       :: pcp(:)

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enables Rainf in forcing variables list')

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((coop_struc(n)%metdata2(1,index1).ne.coop_struc(n)%undef)) then  
        pcp(t)= coop_struc(n)%metdata2(1,index1)
     endif
  enddo
end subroutine timeinterp_coop
