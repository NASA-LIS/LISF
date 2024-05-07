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
! !ROUTINE: timeinterp_mrms_grib
! \label{timeinterp_mrms_grib}
! 
! !REVISION HISTORY:
!  25May2006: Kristi Arsenault; Initial implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  13 Feb 2015: Jonathan Case: Modified for MRMS QPE
!  07 Sep 2017; Jessica Erlingis: Modified for MRMS operational QPE
!
! !INTERFACE:
   subroutine timeinterp_mrms_grib(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod,        only : LIS_rc, LIS_domain
    use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
    use LIS_FORC_AttributesMod
    use LIS_logMod,         only : LIS_verify
    use mrms_grib_forcingMod,    only : mrms_grib_struc

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
    integer :: index1
    integer :: c
    integer          :: status
    type(ESMF_Field) :: pcpField, cpcpField
    real, pointer    :: pcp(:), cpcp(:)
    real,  allocatable :: ratio(:)

    allocate(ratio(LIS_rc%ntiles(n)))

    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status, 'Error: enable Rainf in forcing variables list')
    
    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
         rc=status)
    call LIS_verify(status, 'Error: enable CRainf in forcing variables list')
       
    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)
    
    call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
    call LIS_verify(status)

!------------------------------------------------------------------------
! Compute ratio between convective model precip (field(9)) and total 
!  model precip (field(8)), so that it can be applied to the observed 
!  MRMS precipitation
!------------------------------------------------------------------------
    do c = 1,LIS_rc%ntiles(n)
       if (pcp(c) .ne. 0.0 .and.  &
            pcp(c) .ne. LIS_rc%udef .and.  &
            cpcp(c) .ne. LIS_rc%udef) then
          
          ratio(c) = cpcp(c) / pcp(c) 
          
          if (ratio(c) .gt. 1.0) ratio(c) = 1.0
          if (ratio(c) .lt. 0.0) ratio(c) = 0.0
          
       else
          ratio(c) = 0.0
       endif
    enddo

!-- Assign MRMS precipitation product to LIS variable for use in LSM (rate of mm/s)
    do c = 1,LIS_rc%ntiles(n)
       index1 = LIS_domain(n)%tile(c)%index
       if ( mrms_grib_struc(n)%metdata2(1,index1) .ne. LIS_rc%udef.and. &
            mrms_grib_struc(n)%metdata2(1,index1) >= 0.0 ) then
          pcp(c) = mrms_grib_struc(n)%metdata2(1,index1) / 3600.0  ! mm/s
          cpcp(c) = ratio(c) * pcp(c)           
! J.Case (2/13/2015) -- ensure that small negative values are reset to 0.
          if (pcp(c)  < 0)  pcp(c) = 0.0          
          if (cpcp(c) < 0) cpcp(c) = 0.0
          else
             pcp(c) = LIS_rc%udef
             cpcp(c) = LIS_rc%udef          
       endif
    enddo
    deallocate(ratio)

end subroutine timeinterp_mrms_grib  
