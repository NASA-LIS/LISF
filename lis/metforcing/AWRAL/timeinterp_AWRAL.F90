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
!
! !ROUTINE: timeinterp_AWRAL
! \label{timeinterp_AWRAL}
! 
! !REVISION HISTORY:
!  30 Jan 2017:   Sujay Kumar: Initial version
! !INTERFACE:
  subroutine timeinterp_AWRAL(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod,        only : LIS_rc, LIS_domain
    use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
    use LIS_FORC_AttributesMod
    use LIS_logMod,         only : LIS_verify
    use AWRAL_forcingMod,    only : AWRAL_struc

    implicit none

! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: findex ! AF index

!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP
    integer          :: c, v, index1
    integer          :: status
    type(ESMF_Field) :: tmpField,q2Field,uField,swdownField,swdirField,pcpField
    real, pointer    ::  tmp(:),q2(:),uwind(:),swdown(:),swdir(:),pcp(:)
! ___________________________________________________
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdownField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdirect%varname(1),swdirField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdirect in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdownField,localDE=0,farrayPtr=swdown,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdirField,localDE=0,farrayPtr=swdir,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

   !    reminder what the indexing is for all fields
   !    1'tat',    &
   !    2'avpt',    &
   !    3'rgt',    &
   !    4'radcskyt',    &
   !    5'u2t',    &
   !    6'pt'     /)


   do c = 1,LIS_rc%ntiles(n)
        index1 = LIS_domain(n)%tile(c)%index
        ! test that precip is above or equal to 0 just in case forcing undef not set to LIS undef
        if( AWRAL_struc(n)%metdata2(1,1,index1) .ne.LIS_rc%udef.and. &
              AWRAL_struc(n)%metdata2(1,6,index1) >= 0.0 ) then
            tmp(c) = AWRAL_struc(n)%metdata2(1,1,index1)
            q2(c) = AWRAL_struc(n)%metdata2(1,2,index1)
            swdown(c) = AWRAL_struc(n)%metdata2(1,3,index1)
            swdir(c) = AWRAL_struc(n)%metdata2(1,4,index1) 
            uwind(c) = AWRAL_struc(n)%metdata2(1,5,index1)
            pcp(c) = AWRAL_struc(n)%metdata2(1,6,index1)
         else
            tmp(c) = LIS_rc%udef
            q2(c) = LIS_rc%udef
            swdown(c) = LIS_rc%udef
            swdir(c) = LIS_rc%udef
            uwind(c) = LIS_rc%udef
            pcp(c) = LIS_rc%udef
         endif
   enddo

end subroutine timeinterp_AWRAL
  
