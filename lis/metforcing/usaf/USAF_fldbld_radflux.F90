!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: USAF_fldbld_radflux
! \label{USAF_fldbld_radflux}
!
! !REVISION HISTORY:
! 14 Jun 2016  Refactor into simple caller; see USAF_fldbld_radflux_gfs
!              .........................................James Geiger/NASA
! 11 Oct 2017  Added logic to switch from GALWEM to GFS in emergency
!              .........................................Eric Kemp/GSFC
! 23 Oct 2017  Reorganized code.........................Eric Kemp/GSFC
! 28 May 2024  Generated reader for radiation fluxese...K. Arsenault/SAIC
!
! !INTERFACE:    
!subroutine USAF_fldbld_radflux(n,order,julhr,swdown,longwv)
subroutine USAF_fldbld_radflux(n,order,swdown,longwv)

! !USES: 
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_coreMod,       only : LIS_rc, LIS_masterproc
  use LIS_logMod,        only : LIS_logunit, LIS_endrun, LIS_abort,&
                                LIS_alert
  use LIS_timeMgrMod,    only : LIS_julhr_date, LIS_get_julhr

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: order
!  integer, intent(inout) :: julhr
  real, intent(out)      :: swdown(LIS_rc%lnc(n), LIS_rc%lnr(n))
  real, intent(out)      :: longwv(LIS_rc%lnc(n), LIS_rc%lnr(n))

!
! !DESCRIPTION: 
!  This routine calls the appropriate first guess forcing based on
!  the source of the first guess forcing data.
!
! The arguments and variables are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[order]
!   flag to indicate which data is to be read \newline
!   1 = end time                              \newline
!   2 = start time
!  \item[julhr]
!    julian hour used to determine names of forecast files from previous cycles
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[agrmet\_fldbld\_gfs](\ref{agrmet_fldbld_gfs}) \newline
!   read avn or nogaps data in grib format
!  \item[agrmet\_fldbld\_galwem](\ref{agrmet_fldbld_galwem}) \newline
!   read UK Unified Model (GALWEM) data in GRIB2 format
!  \end{description}
!EOP

  integer :: rc
  character(len=255) :: message(20)
  character(len=10)  :: yyyymmddhh
  integer            :: julhr
  integer            :: jultmp
  integer            :: istart
  integer            :: julend
  integer            :: ordertmp  ! Temporary placeholder
!  integer            :: fc_hr
  integer            :: ierr
  logical :: timetoReadRad
  real*8           :: time
  integer          :: doy
  real :: gmt

  ! Sanity check
  if ( agrmet_struc(n)%first_guess_source .ne. 'GALWEM' .and. &
       agrmet_struc(n)%first_guess_source .ne. "GFS") then
      write(LIS_logunit,*) '[ERR] First guess source is not correctly defined.'
      call LIS_endrun
  end if

  ! Estimate julian hour to trigger reading of radiation fields:
!  call LIS_get_julhr(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,&
!                     0,0,jultmp)
!
!  if(jultmp .gt. agrmet_struc(n)%lastRadHour) then 
!    timeToReadRad = .true. 
!  else
!    timeToReadRad = .false.
!  endif

  timeToReadRad = .true. ! EMK TEST

  if(timeToReadRad) then  
!------------------------------------------------------------------    
! Find the time to start the processing from 
!------------------------------------------------------------------    
     !call find_agrfld_starttime(LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,istart)
     call LIS_get_julhr(LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
          LIS_rc%hr, 0, 0, istart)
     
     !julend = istart+6
     julend = istart

     agrmet_struc(n)%lastRadHour = julend

!    ------------------------------------------------------------------
!    check if the current instance is a restart or continuation of a run
!    if the run is a cold start, read current and previous 6 hour data
!    else copy the current to previous and read the new current data. 
!    ------------------------------------------------------------------
     if(agrmet_struc(n)%findtime1.eq.1.and.agrmet_struc(n)%findtime2.eq.1) then
        ordertmp = 2         
!        order = 2         
!        call AGRMET_fldbld(n,order,istart)
        
!        order = 1   
!        call AGRMET_fldbld(n,order,julend)
     else 
        ordertmp = 1
!        order = 1
!        call AGRMET_fldbld(n,order,julend)
     endif

     ! Generate 10-digit Year-Month-Day-Hour for "julhr":
     call AGRMET_julhr_date10(julend, yyyymmddhh)
     write(LIS_logunit,*) " "
     write(LIS_logunit,*) &
       '[INFO] Searching for NWP radiation fields, valid ',yyyymmddhh

     ! Try fetching GALWEM, if requested
     ierr = 0
     if ( agrmet_struc(n)%first_guess_source == 'GALWEM' ) then

       julhr = julend   ! KRA -- Temporary placeholder

!       fc_hr = LIS_rc%hr  ! KRA - TEMPORARY PLACEHOLDER FOR FC_HR ...
!       write(LIS_logunit,*) '[INFO] NWP radiation fc_hr: ',fc_hr

!       call USAF_fldbld_radflux_galwem(n,julhr,fc_hr,swdown,longwv,rc)
       call USAF_fldbld_radflux_galwem(n,julhr,swdown,longwv,rc)
       ierr = rc
     end if

     ! Try fetching GFS, if requested, or if GALWEM is not available.
     if ( agrmet_struc(n)%first_guess_source == "GFS" .or. &
          ierr .ne. 0) then
       if (ierr .ne. 0) then
          write(LIS_logunit,*)'[WARN] Unable to find GALWEM data!'
          write(LIS_logunit,*)'[WARN] Rolling back to GFS...'
       end if

       call USAF_fldbld_radflux_gfs(n, julhr, swdown, longwv, rc)

       if (rc .ne. 0) then
          call AGRMET_julhr_date10(julhr, yyyymmddhh)
          if (ierr .ne. 0) then
             write(LIS_logunit,*) &
                  '[ERR] No GALWEM or GFS background found for ',yyyymmddhh
          else
             write(LIS_logunit,*) &
                  '[ERR] No GFS background found for ',yyyymmddhh
          end if
          write(LIS_logunit,*) ' ABORTING!'
          flush(LIS_logunit)
          message(:) = ''
          message(1) = '[ERR] Program:  LIS'
          message(2) = '  Routine:  USAF_fldbld_radflux.'
          message(3) = '  GALWEM and GFS GRIB data not available for '//&
             yyyymmddhh
          if(LIS_masterproc) then
             call LIS_alert( 'LIS.USAF_fldbld_radflux.', 1, &
                  message )
             call LIS_abort( message)
          endif
       end if
    end if

  endif   ! End to timeToReadRad check
 
end subroutine USAF_fldbld_radflux
