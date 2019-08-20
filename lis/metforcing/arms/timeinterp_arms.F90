!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_arms
! \label{timeinterp_arms}
! 
!
! !REVISION HISTORY:
! 28 Jun 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine timeinterp_arms(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
  use LIS_FORC_AttributesMod
  use LIS_logMod,         only : LIS_verify
  use arms_forcingMod, only : arms_struc
    
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
! !DESCRIPTION:
! 
! Temporally interpolates the ARMS forcing data to the 
! model timestep. All variables except precipitation is 
! linearly interpolated. 
! 
!EOP
  
  real :: wt1,wt2
  integer :: t
  integer :: index1
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  wt1 = (arms_struc(n)%armstime2-LIS_rc%time) / & 
       (arms_struc(n)%armstime2-arms_struc(n)%armstime1)
  wt2 = 1.0 - wt1
  
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1),vField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),psurfField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1),cpcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable CRainf in the forcing variables list')


  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if((arms_struc(n)%metdata1(1,index1).ne.LIS_rc%udef) .and. &
          (arms_struc(n)%metdata2(1,index1).ne.LIS_rc%udef)) then  
        tmp(t) = &
             wt1 * arms_struc(n)%metdata1(1,index1) +  & 
             wt2 *arms_struc(n)%metdata2(1,index1)
     endif
     if((arms_struc(n)%metdata1(2,index1).ne.LIS_rc%udef) .and. &
          (arms_struc(n)%metdata2(2,index1).ne.LIS_rc%udef)) then  
        q2(t) = &
             wt1 * arms_struc(n)%metdata1(2,index1) +  & 
             wt2 *arms_struc(n)%metdata2(2,index1)
     endif
     if((arms_struc(n)%metdata1(3,index1).ne.LIS_rc%udef) .and. &
          (arms_struc(n)%metdata2(3,index1).ne.LIS_rc%udef)) then  
        swd(t) = &
             wt1 * arms_struc(n)%metdata1(3,index1) +  & 
             wt2 *arms_struc(n)%metdata2(3,index1)
     endif
     if((arms_struc(n)%metdata1(4,index1).ne.LIS_rc%udef) .and. &
          (arms_struc(n)%metdata2(4,index1).ne.LIS_rc%udef)) then  
        lwd(t) = &
             wt1 * arms_struc(n)%metdata1(4,index1) +  & 
             wt2 *arms_struc(n)%metdata2(4,index1)
     endif
     if((arms_struc(n)%metdata1(5,index1).ne.LIS_rc%udef) .and. &
          (arms_struc(n)%metdata2(5,index1).ne.LIS_rc%udef)) then  
        uwind(t) = &
             wt1 * arms_struc(n)%metdata1(5,index1) +  & 
             wt2 *arms_struc(n)%metdata2(5,index1)
     endif
     if((arms_struc(n)%metdata1(6,index1).ne.LIS_rc%udef) .and. &
          (arms_struc(n)%metdata2(6,index1).ne.LIS_rc%udef)) then  
        vwind(t) = &
             wt1 * arms_struc(n)%metdata1(6,index1) +  & 
             wt2 *arms_struc(n)%metdata2(6,index1)
     endif
     if((arms_struc(n)%metdata1(7,index1).ne.LIS_rc%udef) .and. &
          (arms_struc(n)%metdata2(7,index1).ne.LIS_rc%udef)) then  
        psurf(t) = &
             wt1 * arms_struc(n)%metdata1(7,index1) +  & 
             wt2 *arms_struc(n)%metdata2(7,index1)
     endif
     
     if((arms_struc(n)%metdata2(8,index1).ne.LIS_rc%udef)) then  
        pcp(t)= arms_struc(n)%metdata2(8,index1)
     endif
  enddo

end subroutine timeinterp_arms
