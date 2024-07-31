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
! 2 June 2010: Soni Yatheendradas; Initial LIS version for FEWSNET
!
! !INTERFACE:
subroutine timeinterp_RFE2Daily(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use LIS_metforcingMod, only :  LIS_FORC_Base_State, LIS_forc
  use LIS_FORC_AttributesMod
  use LIS_logMod, only : LIS_verify
  use RFE2Daily_forcingMod, only : RFE2Daily_struc
  use LIS_forecastMod, only : LIS_get_iteration_index

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

    integer          :: t,index1
    integer          :: mfactor, m, k, kk
    integer          :: status
    type(ESMF_Field) :: pcpField, cpcpField
    real, pointer    :: pcp(:), cpcp(:)

    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status, 'Error: enable Rainf in forcing variables list')

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)

!-- Convert precipitation sum to rate:
   mfactor = LIS_rc%nensem(n)/RFE2Daily_struc(n)%nIter

   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
       t = m + (k-1)*mfactor
       index1 = LIS_domain(n)%tile(t)%index
       kk = LIS_get_iteration_index(n, k, index1, mfactor)
       ! Leave out undefined and NODATA value grid cells if after upscaling:
       if( RFE2Daily_struc(n)%metdata2(kk,1,index1) .ge. 0. ) then
          pcp(t) = RFE2Daily_struc(n)%metdata2(kk,1,index1) / 86400.0 ! KRA: mm/s for entire day
       endif
     enddo
   enddo

end subroutine timeinterp_RFE2Daily
