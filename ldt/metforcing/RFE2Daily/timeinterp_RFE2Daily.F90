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
! !ROUTINE: timeinterp_RFE2Daily
! \label{timeinterp_RFE2Daily}
!
! !REVISION HISTORY:
!
!  2 June 2010: Soni Yatheendradas; Initial LDT version for FEWSNET
! 16 June 2016: KR Arsenault; Updated code to more correctly handle in LIS
!
! !INTERFACE:
subroutine timeinterp_RFE2Daily(n,findex)
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc,LDT_domain
  use LDT_metforcingMod, only :  LDT_FORC_Base_State, LDT_forc
  use LDT_FORC_AttributesMod
  use LDT_logMod, only : LDT_verify
  use RFE2Daily_forcingMod, only : RFE2Daily_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model
!  timestep (evenly distributed over sub-daily timesteps for now). 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the supplemental forcing source
!  \item[metdata1]
!    previous forcing values
!  \item[metdata2]
!    next forcing values
!  \end{description}
!
!EOP

    integer          :: t, index1
    integer          :: status
    type(ESMF_Field) :: pcpField, cpcpField
    real, pointer    :: pcp(:), cpcp(:)
    real,  allocatable :: ratio(:)

    call ESMF_StateGet(LDT_FORC_Base_State(n,findex),trim(LDT_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LDT_verify(status, 'Error: enable Rainf in forcing variables list')

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LDT_verify(status)

!-- Convert precipitation sum to rate:
    do t=1,LDT_rc%ntiles(n)
       index1 = LDT_domain(n)%tile(t)%index
     ! Leave out undefined and NODATA value grid cells if after upscaling:
       if( LDT_forc(n,findex)%metdata2(1,index1) .ge. 0.0 ) then
          pcp(t) = LDT_forc(n,findex)%metdata2(1,index1) / 86400.0 ! KRA: mm/s for entire day
       else
          pcp(t) = LDT_rc%udef 
       endif
    enddo

#if 0
!-- Applying a convective rainfall amount to observed precip field: 
    if( LDT_FORC_CRainf%selectOpt .eq. 1 ) then

      allocate(ratio(LDT_rc%ntiles(n)))

      call ESMF_StateGet(LDT_FORC_Base_State(n,findex),trim(LDT_FORC_CRainf%varname(1)),cpcpField,&
           rc=status)
      call LDT_verify(status, 'Error: enable CRainf in forcing variables list')

      call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
      call LDT_verify(status)

!------------------------------------------------------------------------
! Compute ratio between convective model precip and total model precip
! so that it can be applied to the observed global precip
!------------------------------------------------------------------------
      do t = 1,LDT_rc%ntiles(n)
         if( pcp(t) .ne. 0.0 .and.  &
             pcp(t) .ne. LDT_rc%udef .and.  &
             cpcp(t) .ne. LDT_rc%udef) then
           ratio(t) = cpcp(t) / pcp(t)
           if (ratio(t) .gt. 1.0) ratio(t) = 1.0
           if (ratio(t) .lt. 0.0) ratio(t) = 0.0
           cpcp(t) = ratio(t) * pcp(t)
         else
           cpcp(t) = 0.0
         endif
      enddo
      deallocate(ratio)
    endif
#endif

end subroutine timeinterp_RFE2Daily
