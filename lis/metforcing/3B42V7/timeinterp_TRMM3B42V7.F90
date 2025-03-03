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
    use LIS_coreMod, only           : LIS_rc, LIS_domain
    use LIS_FORC_AttributesMod
    use LIS_metforcingMod,  only    : LIS_FORC_Base_State, LIS_forc
    use LIS_logMod,         only    : LIS_verify, LIS_logunit 
    use TRMM3B42V7_forcingMod, only : TRMM3B42V7_struc
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

   integer          :: t, index1
   integer          :: status
   type(ESMF_Field) :: pcpField, cpcpField
   real, pointer    :: pcp(:), cpcp(:)
   real             :: wt1, wt2       
   real*8           :: BreakpointTime ! SY

   integer          :: m, mfactor, kk, k
! _______________________________________

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
     wt1 = (abs(BreakpointTime-TRMM3B42V7_struc(n)%LIS_timeAtTStepStart)) / &
           (TRMM3B42V7_struc(n)%LIS_timeAtTStepEnd- &
            TRMM3B42V7_struc(n)%LIS_timeAtTStepStart) ! SY
     wt2 = 1.0 - wt1
   endif
    
   ! Read in precip:
   call ESMF_StateGet(LIS_FORC_Base_State(n,findex),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
        rc=status)
   call LIS_verify(status, 'Error: enable Rainf in forcing variables list')

   call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
   call LIS_verify(status)

   mfactor = LIS_rc%nensem(n)/TRMM3B42V7_struc(n)%nIter
   do k=1,LIS_rc%ntiles(n)/mfactor
     do m=1,mfactor
        t = m + (k-1)*mfactor
        index1 = LIS_domain(n)%tile(t)%index
        kk = LIS_get_iteration_index(n, k, index1, mfactor)

        ! J.Case (4/22/2013) -- make consistent with Stg4/NMQ forcing examples.
        ! SY: Further modified to allow for averaging over 2 TRMM datasets
 !       if (TRMM3B42V7_struc(n)%metdata2(kk,1,index1) .ge. 0.0) then ! J.Case
 !       if ( TRMM3B42V7_struc(n)%metdata2(kk,1,index1) .ne. LIS_rc%udef .and. & ! J.Case ! SY
 !           TRMM3B42V7_struc(n)%metdata2(kk,1,index1) >= 0.0 ) then ! J.Case ! SY
        ! SY: Start modification
        if ( ( TRMM3B42V7_struc(n)%metdata2(kk,1,index1) .ne. LIS_rc%udef .and. &
             TRMM3B42V7_struc(n)%metdata2(kk,1,index1) >= 0.0 ) .OR. &
           ( TRMM3B42V7_struc(n)%metdata1(kk,1,index1) .ne. LIS_rc%udef .and. &
             TRMM3B42V7_struc(n)%metdata1(kk,1,index1) >= 0.0 ) ) then
           if ( TRMM3B42V7_struc(n)%metdata2(kk,1,index1) .eq. LIS_rc%udef .or. &
               TRMM3B42V7_struc(n)%metdata2(kk,1,index1) < 0.0 ) then
             pcp(t) = TRMM3B42V7_struc(n)%metdata1(kk,1,index1) / 3600.0       ! mm/s
           elseif ( TRMM3B42V7_struc(n)%metdata1(kk,1,index1) .eq. LIS_rc%udef .or. &
                   TRMM3B42V7_struc(n)%metdata1(kk,1,index1) < 0.0 ) then
             pcp(t) = TRMM3B42V7_struc(n)%metdata2(kk,1,index1) / 3600.0   ! mm/s
           else
             pcp(t) = (wt1 * TRMM3B42V7_struc(n)%metdata1(kk,1,index1) &
                    + wt2 * TRMM3B42V7_struc(n)%metdata2(kk,1,index1)) / 3600.0   ! mm/s
           endif
           ! SY: End modification
        else
           pcp(t) = LIS_rc%udef
        endif
        ! J.Case -- ensure that small negative values are reset to 0.
        if ( pcp(t) /= LIS_rc%udef ) then
           if (pcp(t)  < 0)  pcp(t) = 0.0
        endif
     enddo
   enddo

 end subroutine timeinterp_TRMM3B42V7
