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
! !ROUTINE: timeinterp_HiMATGMU
! \label{timeinterp_HiMATGMU}
! 
! !REVISION HISTORY:
!  28 July 2017: Sujay Kumar; Initial implementation
!
! !INTERFACE:
  subroutine timeinterp_HiMATGMU(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod,        only : LIS_rc, LIS_domain
    use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
    use LIS_FORC_AttributesMod
    use LIS_logMod,         only : LIS_verify
    use HiMATGMU_forcingMod,    only : HiMATGMU_struc

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

    do c = 1,LIS_rc%ntiles(n)
       index1 = LIS_domain(n)%tile(c)%index
       if( HiMATGMU_struc(n)%metdata2(1,index1) .ne.LIS_rc%udef.and. &
            HiMATGMU_struc(n)%metdata2(1,index1) >= 0.0 ) then
          pcp(c) = HiMATGMU_struc(n)%metdata2(1,index1) / 3600.0  ! mm/s
       else
          pcp(c) = LIS_rc%udef
       endif
    enddo

end subroutine timeinterp_HiMATGMU
  
