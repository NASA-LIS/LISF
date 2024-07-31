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
  use LDT_coreMod,        only : LDT_rc,LDT_domain
  use LDT_FORC_AttributesMod
  use LDT_metforcingMod,  only : LDT_FORC_Base_State, LDT_forc
  use LDT_logMod,         only : LDT_verify
  use cmorph_forcingMod,  only : cmorph_struc

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

  real    :: ratio(LDT_rc%ntiles(n))
  integer :: t,index1
  integer          :: status
  type(ESMF_Field) :: pcpField, cpcpField
  real, pointer    :: pcp(:), cpcp(:)
  
  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),trim(LDT_FORC_Rainf%varname(1)),&
       pcpField,rc=status)
  call LDT_verify(status, 'Error: enable Rainf in forcing variables list')

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LDT_verify(status)
  
!- Convert precipitation sum to rate:
   do t=1,LDT_rc%ntiles(n)
      index1 = LDT_domain(n)%tile(t)%index
      if (LDT_forc(n,findex)%metdata2(1,index1) .ge. 0.0) then
         pcp(t) = LDT_forc(n,findex)%metdata2(1,index1) / 3600.0  ! mm/s
      else
         pcp(t) = LDT_rc%udef
      endif
   enddo

#if 0
!-- Applying a convective rainfall amount to observed precip field: 
   if( LDT_FORC_CRainf%selectOpt == 1 ) then

     call ESMF_StateGet(LDT_FORC_Base_State(n,findex),trim(LDT_FORC_CRainf%varname(1)),&
          cpcpField,rc=status)
     call LDT_verify(status, 'Error: enable CRainf in forcing variables list')
   
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LDT_verify(status)

!------------------------------------------------------------------------
! Compute ratio between convective model precip and total model precip
! so that it can be applied to the observed global precip
!------------------------------------------------------------------------
     do t=1,LDT_rc%ntiles(n)
       if (pcp(t) .ne. 0.0 .and.  &
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
  endif
#endif

end subroutine timeinterp_cmorph
  
