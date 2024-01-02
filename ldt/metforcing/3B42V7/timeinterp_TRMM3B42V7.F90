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
! !ROUTINE: timeinterp_TRMM3B42V7
! \label{timeinterp_TRMM3B42V7}
!
! !REVISION HISTORY: 
!  20 Jan 2006:   Yudong Tian: Initial Implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  19 Apr 2013: Jonathan Case; Included comparison against undefined value,
!                              set small -ve values to 0
!  21 Jun 2013: Soni Yatheendradas; Based on new TRMM3B42V6 code, which is
!               based on changes from earlier code to avoid (a) alternate file skip,
!               (b) jump to previous day TRMM, and (c) absence of rain rate weighting
! !INTERFACE: 
 subroutine timeinterp_TRMM3B42V7(n,findex)
! !USES:
    use ESMF
    use LDT_coreMod, only           : LDT_rc, LDT_domain
    use LDT_FORC_AttributesMod
    use LDT_metforcingMod,  only    : LDT_FORC_Base_State, LDT_forc
    use LDT_logMod,         only    : LDT_verify, LDT_logunit 
    use TRMM3B42V7_forcingMod, only : TRMM3B42V7_struc

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

    integer :: t, index1

    integer          :: status
    type(ESMF_Field) :: pcpField, cpcpField
    real, pointer    :: pcp(:), cpcp(:)
    real             :: wt1, wt2       ! SY
    real*8           :: BreakpointTime ! SY
    real,  allocatable :: ratio(:)


    ! SY: Begin calculating weights, note that weights are NOT linear interpolation-based
    if ((TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepEnd .EQ. &
         TRMM3B42V7_struc(n)%TRMM3B42V7yr_TStepStart) .AND. &
        (TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepEnd .EQ. &
         TRMM3B42V7_struc(n)%TRMM3B42V7mo_TStepStart) .AND. &
        (TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepEnd .EQ. &
         TRMM3B42V7_struc(n)%TRMM3B42V7da_TStepStart) .AND. &
        (TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepEnd .EQ. &
         TRMM3B42V7_struc(n)%TRMM3B42V7hr_TStepStart)) then
      wt1 = 0.5
      wt2 = 0.5
    else
      BreakpointTime = (TRMM3B42V7_struc(n)%TRMM3B42V7time_TStepStart+&
                        TRMM3B42V7_struc(n)%TRMM3B42V7time_TStepEnd)/2.0
      wt1 = (abs(BreakpointTime-TRMM3B42V7_struc(n)%LDT_timeAtTStepStart)) / &
            (TRMM3B42V7_struc(n)%LDT_timeAtTStepEnd- &
             TRMM3B42V7_struc(n)%LDT_timeAtTStepStart) ! SY
      wt2 = 1.0 - wt1
    endif
    ! SY: End calculating weights, note that weights are NOT linear interpolation-based
     
    call ESMF_StateGet(LDT_FORC_Base_State(n,findex),trim(LDT_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LDT_verify(status, 'Error: enable Rainf in forcing variables list')

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LDT_verify(status)

!- Convert precipitation sum to rate:
!   if( LDT_FORC_CRainf%selectOpt .ne. 1 ) then

     do t = 1,LDT_rc%ntiles(n)
        index1 = LDT_domain(n)%tile(t)%index
        ! SY: Start modification
        if ( ( LDT_forc(n,findex)%metdata2(1,index1) .ne. LDT_rc%udef .and. &
               LDT_forc(n,findex)%metdata2(1,index1) >= 0.0 ) .OR. &
             ( LDT_forc(n,findex)%metdata1(1,index1) .ne. LDT_rc%udef .and. &
               LDT_forc(n,findex)%metdata1(1,index1) >= 0.0 ) ) then
           if ( LDT_forc(n,findex)%metdata2(1,index1) .eq. LDT_rc%udef .or. &
                LDT_forc(n,findex)%metdata2(1,index1) < 0.0 ) then
             pcp(t) = LDT_forc(n,findex)%metdata1(1,index1) / 3600.0      ! mm/s
           elseif ( LDT_forc(n,findex)%metdata1(1,index1) .eq. LDT_rc%udef .or. &
                LDT_forc(n,findex)%metdata1(1,index1) < 0.0 ) then
             pcp(t) = LDT_forc(n,findex)%metdata2(1,index1) / 3600.0      ! mm/s
           else
             pcp(t) = (wt1 * LDT_forc(n,findex)%metdata1(1,index1) + wt2 * &
                       LDT_forc(n,findex)%metdata2(1,index1)) / 3600.0    ! mm/s
           endif
         ! SY: End modification
        else
           pcp(t) = LDT_rc%udef
        endif
      ! J.Case -- ensure that small negative values are reset to 0.
        if ( pcp(t) /= LDT_rc%udef ) then
           if (pcp(t)  < 0)  pcp(t) = 0.0
        endif
     enddo

#if 0
!-- Applying a convective rainfall amount to observed precip field: 
   elseif( LDT_FORC_CRainf%selectOpt == 1 ) then

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
        if (pcp(t) .ne. 0.0 .and.  &
             pcp(t) .ne. LDT_rc%udef .and.  &
             cpcp(t) .ne. LDT_rc%udef) then
           ratio(t) = cpcp(t) / pcp(t)
           if (ratio(t) .gt. 1.0) ratio(t) = 1.0
           if (ratio(t) .lt. 0.0) ratio(t) = 0.0
        else
           ratio(t) = 0.0
        endif
     enddo
 
     do t = 1,LDT_rc%ntiles(n) 
        index1 = LDT_domain(n)%tile(t)%index
 ! J.Case (4/22/2013) -- make consistent with Stg4/NMQ forcing examples.
 ! SY: Further modified to allow for averaging over 2 TRMM datasets
 !       if (LDT_forc(n,findex)%metdata2(1,index1) .ge. 0.0) then ! J.Case
 !       if ( LDT_forc(n,findex)%metdata2(1,index1) .ne. LDT_rc%udef .and. & ! J.Case ! SY
 !            LDT_forc(n,findex)%metdata2(1,index1) >= 0.0 ) then ! J.Case ! SY
 ! SY: Start modification
        if ( ( LDT_forc(n,findex)%metdata2(1,index1) .ne. LDT_rc%udef .and. &
               LDT_forc(n,findex)%metdata2(1,index1) >= 0.0 ) .OR. &
             ( LDT_forc(n,findex)%metdata1(1,index1) .ne. LDT_rc%udef .and. &
               LDT_forc(n,findex)%metdata1(1,index1) >= 0.0 ) ) then
           if ( LDT_forc(n,findex)%metdata2(1,index1) .eq. LDT_rc%udef .or. &
                LDT_forc(n,findex)%metdata2(1,index1) < 0.0 ) then
             pcp(t) = LDT_forc(n,findex)%metdata1(1,index1) / 3600.0       ! mm/s
           elseif ( LDT_forc(n,findex)%metdata1(1,index1) .eq. LDT_rc%udef .or. &
                LDT_forc(n,findex)%metdata1(1,index1) < 0.0 ) then
             pcp(t) = LDT_forc(n,findex)%metdata2(1,index1) / 3600.0   ! mm/s
           else
             pcp(t) = (wt1 * LDT_forc(n,findex)%metdata1(1,index1) &
                    + wt2 * LDT_forc(n,findex)%metdata2(1,index1)) / 3600.0       ! mm/s
           endif
 ! SY: End modification
           !pcp(t) = LDT_forc(n,findex)%metdata2(1,index1) / 3600.0       ! mm/s ! SY
           cpcp(t) = ratio(t) * pcp(t)
        else
           pcp(t) = LDT_rc%udef
           cpcp(t) = LDT_rc%udef
        endif
 ! J.Case -- ensure that small negative values are reset to 0.
        if ( pcp(t) /= LDT_rc%udef ) then
           if (pcp(t)  < 0)  pcp(t) = 0.0
        endif
        if ( cpcp(t) /= LDT_rc%udef ) then
           if (cpcp(t) < 0) cpcp(t) = 0.0
        endif
     enddo
     deallocate(ratio)

   endif
#endif

  end subroutine timeinterp_TRMM3B42V7
