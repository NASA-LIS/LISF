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
! !ROUTINE: timeinterp_stg4
! \label{timeinterp_stg4}
! 
! !REVISION HISTORY:
!  25May2006: Kristi Arsenault; Initial implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
! !INTERFACE:
  subroutine timeinterp_stg4(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod,        only : LIS_rc, LIS_domain
    use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
    use LIS_FORC_AttributesMod
    use LIS_logMod,         only : LIS_verify
    use stg4_forcingMod,    only : stg4_struc

    implicit none

! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: findex

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
    integer          :: c, index1
    integer          :: status
    type(ESMF_Field) :: pcpField, cpcpField
    real, pointer    :: pcp(:), cpcp(:)
    real,  allocatable :: ratio(:)
! ___________________________________________________

    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
         trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status, 'Error: enable Rainf in forcing variables list')
    
    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)

!-- Convert precipitation sum to rate:
!    if( LIS_FORC_CRainf%selectOpt .ne. 1 ) then

   !- Assign STAGE IV rainfall product to LIS variable for use in LSM (rate of mm/s)
      do c = 1,LIS_rc%ntiles(n)
         index1 = LIS_domain(n)%tile(c)%index
         if( stg4_struc(n)%metdata2(1,index1) .ne.LIS_rc%udef.and. &
             stg4_struc(n)%metdata2(1,index1) >= 0.0 ) then
            pcp(c) = stg4_struc(n)%metdata2(1,index1) / 3600.0  ! mm/s
         else
            pcp(c) = LIS_rc%udef
         endif
      enddo

#if 0
!-- Applying a convective rainfall amount to observed precip field: 
    elseif( LIS_FORC_CRainf%selectOpt == 1 ) then

      allocate(ratio(LIS_rc%ntiles(n)))

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex), &
           trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
           rc=status)
      call LIS_verify(status, 'Error: enable CRainf in forcing variables list')
    
      call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
      call LIS_verify(status)

!------------------------------------------------------------------------
! Compute ratio between convective model precip (field(9)) and total 
!  model precip (field(8)), so that it can be applied to the observed 
!  STAGE IV precipitation
!------------------------------------------------------------------------
      do c = 1,LIS_rc%ntiles(n)
         if( pcp(c) .ne. 0.0 .and. pcp(c) .ne. LIS_rc%udef .and.  &
             cpcp(c) .ne. LIS_rc%udef) then
           ratio(c) = cpcp(c) / pcp(c) 
           if (ratio(c) .gt. 1.0) ratio(c) = 1.0
           if (ratio(c) .lt. 0.0) ratio(c) = 0.0
           cpcp(c) = ratio(c) * pcp(c)           
         else
           cpcp(c) = 0.0
         endif
      enddo
      deallocate(ratio)
   endif
#endif

end subroutine timeinterp_stg4
  
