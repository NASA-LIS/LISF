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
! !ROUTINE: timeinterp_stg2
! \label{timeinterp_stg2}
! 
! !REVISION HISTORY:
!  25May2006: Kristi Arsenault; Initial implementation
!  10Oct2006: Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
! !INTERFACE:
  subroutine timeinterp_stg2(n,findex)
! !USES:
    use ESMF
    use LDT_coreMod,        only : LDT_rc, LDT_domain
    use LDT_metforcingMod,  only : LDT_FORC_Base_State, LDT_forc
    use LDT_FORC_AttributesMod
    use LDT_logMod,         only : LDT_verify
    use stg2_forcingMod,    only : stg2_struc

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
    real, allocatable :: ratio(:)
! ________________________________________________________
 
    call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
         trim(LDT_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LDT_verify(status, 'Error: enable Rainf in forcing variables list')

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LDT_verify(status)
    
!-- Convert precipitation sum to rate:
    do c = 1,LDT_rc%ntiles(n)
       index1 = LDT_domain(n)%tile(c)%index
       if( LDT_forc(n,findex)%metdata2(1,index1) .ne. LDT_rc%udef .and. &
          LDT_forc(n,findex)%metdata2(1,index1) >= 0.0 ) then
         pcp(c) = LDT_forc(n,findex)%metdata2(1,index1) / 3600.0  ! mm/s
       else
         pcp(c) = LDT_rc%udef
       endif
    enddo

#if 0
!-- Applying a convective rainfall amount to observed precip field: 
    if( LDT_FORC_CRainf%selectOpt == 1 ) then

      allocate(ratio(LDT_rc%ntiles(n)))

      call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
           trim(LDT_FORC_CRainf%varname(1)),cpcpField,&
           rc=status)
      call LDT_verify(status, 'Error: enable CRainf in forcing variables list')
    
      call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
      call LDT_verify(status)

!------------------------------------------------------------------------
! Compute ratio between convective model precip (field(9)) and total 
!  model precip (field(8)), so that it can be applied to the observed 
!  STAGE II precipitation
!------------------------------------------------------------------------
      do c = 1,LDT_rc%ntiles(n)
         if( pcp(c) .ne. 0.0 .and. pcp(c) .ne. LDT_rc%udef .and.  &
              cpcp(c) .ne. LDT_rc%udef) then
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

end subroutine timeinterp_stg2
