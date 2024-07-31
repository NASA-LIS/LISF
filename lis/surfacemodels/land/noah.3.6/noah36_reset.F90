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
! !ROUTINE: noah36_reset
! \label{noah36_reset}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial version
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
! 
! !INTERFACE:
subroutine noah36_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use noah36_lsmMod
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
       LIS_update_timestep
  use LIS_logMod,       only : LIS_verify, LIS_logunit

!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Noah-3.6 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Noah-3.6.
!  
! The routines invoked are: 
! \begin{description}
! \item[noah36\_resetvegparms](\ref{noah36_resetvegparms}) \newline
!   initializes the vegetation-related parameters in Noah-3.6
! \item[noah36\_resettbot](\ref{noah36_resettbot}) \newline
!   initializes the bottom temperature fields
! \item[noah36\_resetsoils](\ref{noah36_resetsoils}) \newline
!   initializes the soil parameters
! \end{description}
!EOP
  implicit none
  integer                 :: i,n
  integer                 :: status


    do n=1,LIS_rc%nnest
       write(LIS_logunit,*)                        &
            'Noah-3.6 resetting'

       noah36_struc(n)%forc_count = 0 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noah36_struc(n)%noah(i)%tair = 0 
          noah36_struc(n)%noah(i)%qair = 0 
          noah36_struc(n)%noah(i)%swdown = 0 
          noah36_struc(n)%noah(i)%lwdown = 0
          noah36_struc(n)%noah(i)%uwind = 0
          noah36_struc(n)%noah(i)%vwind = 0 
          noah36_struc(n)%noah(i)%psurf = 0 
          noah36_struc(n)%noah(i)%rainf = 0 
          noah36_struc(n)%noah(i)%snowf = 0 
          noah36_struc(n)%noah(i)%rainf_c = 0 
          noah36_struc(n)%noah(i)%ch = 0 
       enddo

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
!       call LIS_update_timestep(LIS_rc, n, noah36_struc(n)%ts)


       ! Initialize min/max values to implausible values.
       noah36_struc(n)%noah(:)%tair_agl_min = 999.0


       noah36_struc(n)%z0brd_upd = 0 
       noah36_struc(n)%lai_upd = 0 
       noah36_struc(n)%embrd_upd = 0 
       noah36_struc(n)%alb_upd = 0 

       noah36_struc(n)%optStartFlag = 1 
!uninitialized variables: 
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          noah36_struc(n)%noah(i)%cqs2 = LIS_rc%udef
          noah36_struc(n)%noah(i)%qsfc = LIS_rc%udef
          noah36_struc(n)%noah(i)%sca = 0.0
       enddo
       
       
!!!       LIS_sfmodel_struc(n)%nsm_layers = noah36_struc(n)%nslay
!!!       LIS_sfmodel_struc(n)%nst_layers = noah36_struc(n)%nslay
!!!       allocate(LIS_sfmodel_struc(n)%lyrthk(noah36_struc(n)%nslay))
!!!       do i = 1,noah36_struc(n)%nslay
!!!          LIS_sfmodel_struc(n)%lyrthk(i) = noah36_struc(n)%lyrthk(i)*100.0
!!!       enddo
!!!       LIS_sfmodel_struc(n)%ts = noah36_struc(n)%ts

       call noah36_resetvegparms(LIS_rc%lsm_index)
!       call noah36_resettbot(LIS_rc%lsm_index)
!       call noah36_resetsoils(LIS_rc%lsm_index)
       
  enddo
end subroutine noah36_reset
