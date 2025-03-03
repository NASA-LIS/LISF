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
! !ROUTINE: clsmf25_reset
! \label{clsmf25_reset}
!
! !REVISION HISTORY:
! 15 Dec 2005; Sujay Kumar, Initial Code
! 23 Nov 2012: David Mocko, Additional configs for Catchment Fortuna-2.5
! 
! !INTERFACE:
subroutine clsmf25_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use clsmf25_lsmMod
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
       LIS_update_timestep
  use LIS_logMod,       only : LIS_verify, LIS_logunit

!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for CLSMF2.5 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in CLSMF2.5.
!  
! The routines invoked are: 
! \begin{description}
! \item[clsmf25\_resetvegparms](\ref{clsmf25_resetvegparms}) \newline
!   initializes the vegetation-related parameters in CLSMF2.5
! \item[clsmf25\_resettbot](\ref{clsmf25_resettbot}) \newline
!   initializes the bottom temperature fields
! \item[clsmf25\_resetsoils](\ref{clsmf25_resetsoils}) \newline
!   initializes the soil parameters
! \end{description}
!EOP
   implicit none
   integer                 :: i,n
   integer                 :: status

   do n=1,LIS_rc%nnest
      write(LIS_logunit,*)  &
           'CLSMF2.5 resetting'

      clsmf25_struc(n)%forc_count = 0 
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
         clsmf25_struc(n)%met_force(i)%tair = 0 
         clsmf25_struc(n)%met_force(i)%qair = 0 
         clsmf25_struc(n)%met_force(i)%swdown = 0 
         clsmf25_struc(n)%met_force(i)%lwdown = 0
         clsmf25_struc(n)%met_force(i)%wind = 0
         clsmf25_struc(n)%met_force(i)%psurf = 0 
         clsmf25_struc(n)%met_force(i)%rainf = 0 
         clsmf25_struc(n)%met_force(i)%snowf = 0 
         clsmf25_struc(n)%met_force(i)%rainf_c = 0 
      enddo
   enddo

end subroutine clsmf25_reset
