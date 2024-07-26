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
! !ROUTINE: timeinterp_cmap
! \label{timeinterp_cmap}
! 
! !REVISION HISTORY: 
!  20 Jan 2006:   Yudong Tian: Initial Implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  02 Dec 2014:  KR Arsenault: Corrected undefined mask values
!
! !INTERFACE:
  subroutine timeinterp_cmap(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod,       only : LIS_rc, LIS_domain
    use LIS_metforcingMod, only : LIS_FORC_Base_State, LIS_forc
    use LIS_FORC_AttributesMod
    use LIS_logMod,        only : LIS_verify
    use cmap_forcingMod,   only : cmap_struc

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
    integer           :: t,index1
    integer           :: status
    type(ESMF_Field)  :: pcpField, cpcpField
    real, pointer     :: pcp(:), cpcp(:)
    real, allocatable :: ratio(:)


    call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
         trim(LIS_FORC_Rainf%varname(1)),&
         pcpField,rc=status)
    call LIS_verify(status, 'Error: enable Rainf in forcing variables list')
    
    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)

!-- Convert precipitation sum to rate:
!    if( LIS_FORC_CRainf%selectOpt .ne. 1 ) then

      do t=1,LIS_rc%ntiles(n)
         index1 = LIS_domain(n)%tile(t)%index
!         if (cmap_struc(n)%metdata2(1,index1) .ne. -1.0) then
         if( cmap_struc(n)%metdata2(1,index1) .ge. 0. ) then
            pcp(t) = cmap_struc(n)%metdata2(1,index1) !/ 3600.0
         else
            pcp(t) = LIS_rc%udef
         endif
      enddo

#if 0
!-- Applying a convective rainfall amount to observed precip field: 
    elseif( LIS_FORC_CRainf%selectOpt == 1 ) then

      allocate(ratio(LIS_rc%ntiles(n)))

      call ESMF_StateGet(LIS_FORC_Base_State(n,findex),&
           trim(LIS_FORC_CRainf%varname(1)),&
           cpcpField,rc=status)
      call LIS_verify(status, 'Error: enable CRainf in forcing variables list')

      call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
      call LIS_verify(status)

!------------------------------------------------------------------------
! Compute ratio between convective model precip and total model precip
! so that it can be applied to the observed global precip
!------------------------------------------------------------------------
      do t = 1,LIS_rc%ntiles(n)
         if (pcp(t) .ne. 0.0 .and.  & 
             pcp(t) .ne. LIS_rc%udef .and.  & 
             cpcp(t) .ne. LIS_rc%udef) then
            ratio(t) = cpcp(t) / pcp(t)
            if (ratio(t) .gt. 1.0) ratio(t) = 1.0
            if (ratio(t) .lt. 0.0) ratio(t) = 0.0
         else
            ratio(t) = 0.0
         endif
      enddo

      do t=1,LIS_rc%ntiles(n)
         index1 = LIS_domain(n)%tile(t)%index
!         if( cmap_struc(n)%metdata2(1,index1) .ne. -1.0 ) then
         if( cmap_struc(n)%metdata2(1,index1) .ge. 0. ) then
            pcp(t) = cmap_struc(n)%metdata2(1,index1) !/ 3600.0
            cpcp(t) = ratio(t) * pcp(t)
         else
            pcp(t) = LIS_rc%udef
            cpcp(t) = LIS_rc%udef
         endif
      enddo
      deallocate(ratio)
    endif
#endif

  end subroutine timeinterp_cmap
