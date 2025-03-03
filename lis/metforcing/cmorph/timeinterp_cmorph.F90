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
! !ROUTINE: timeinterp_cmorph
! \label{timeinterp_cmorph}
! 
! !REVISION HISTORY: 
!  20 Jan 2006:   Yudong Tian: Initial Implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
! !INTERFACE:
subroutine timeinterp_cmorph(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc,LIS_domain
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,  only : LIS_FORC_Base_State, LIS_forc
  use LIS_logMod,         only : LIS_verify
  use cmorph_forcingMod,  only : cmorph_struc
  use LIS_forecastMod, only : LIS_get_iteration_index

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
!   index of the nest
!  \end{description}
!
!EOP

  integer          :: t,index1
  integer          :: mfactor, m, k, kk
  integer          :: status
  type(ESMF_Field) :: pcpField, cpcpField
  real, pointer    :: pcp(:), cpcp(:)
  
  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),&
       pcpField,rc=status)
  call LIS_verify(status, 'Error: enable Rainf in forcing variables list')

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)
  
  ! Convert precipitation sum to rate:
   mfactor = LIS_rc%nensem(n)/cmorph_struc(n)%nIter

   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)

        if (cmorph_struc(n)%metdata2(kk,1,index1) .ge. 0.0) then
           pcp(t) = cmorph_struc(n)%metdata2(kk,1,index1) / 3600.0  ! mm/s
        else
           pcp(t) = LIS_rc%udef
        endif
     enddo
   end do

end subroutine timeinterp_cmorph
  
