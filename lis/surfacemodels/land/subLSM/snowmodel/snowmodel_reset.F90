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
! !ROUTINE: snowmodel_reset
! \label{snowmodel_reset}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel
! 
! !INTERFACE:
subroutine snowmodel_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use snowmodel_lsmMod
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
       LIS_update_timestep
  use LIS_logMod,       only : LIS_verify, LIS_logunit
  use snowmodel_lsmMod
!
! !DESCRIPTION: 
! 
!  This routine is the reset point for parameters and variables
!  required for SnowModel. 
!  
! The routines invoked are: 
! \begin{description}
! \item[snowmodel\_resetvegparms](\ref{snowmodel_resetvegparms}) \newline
!   initializes the vegetation-related parameters in SnowModel
! \end{description}
!EOP
  implicit none
  integer          :: i,n


    do n=1,LIS_rc%nnest
       write(LIS_logunit,*)  &
            '[INFO] SnowModel resetting'
       write(LIS_logunit,*) "[INFO] Reset on:", LIS_rc%mo, LIS_rc%da, &
             LIS_rc%hr, snowmodel_struc(n)%forc_count

       snowmodel_struc(n)%forc_count = 0
       do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          snowmodel_struc(n)%sm(i)%tair = 0
          snowmodel_struc(n)%sm(i)%qair = 0
          snowmodel_struc(n)%sm(i)%swdown = 0
          snowmodel_struc(n)%sm(i)%lwdown = 0
          snowmodel_struc(n)%sm(i)%uwind = 0
          snowmodel_struc(n)%sm(i)%vwind = 0
          snowmodel_struc(n)%sm(i)%psurf = 0
          snowmodel_struc(n)%sm(i)%rainf = 0
!          snowmodel_struc(n)%sm(i)%snowf = 0
       enddo
 
    enddo

end subroutine snowmodel_reset
