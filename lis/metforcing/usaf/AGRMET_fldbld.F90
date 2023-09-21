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
! !ROUTINE: AGRMET_fldbld
! \label{AGRMET_fldbld}
!
! !REVISION HISTORY:
! 14 Jun 2016  Refactor into simple caller; see AGRMET_fldbld_gfs
!              .........................................James Geiger/NASA
! 11 Oct 2017  Added logic to switch from GALWEM to GFS in emergency
!              .........................................Eric Kemp/GSFC
! 23 Oct 2017  Reorganized code.........................Eric Kemp/GSFC
! !INTERFACE:    
subroutine AGRMET_fldbld(n,order,julhr)
! !USES: 
  use AGRMET_forcingMod, only : agrmet_struc
  use LIS_coreMod,       only : LIS_masterproc
  use LIS_logMod,        only : LIS_logunit, LIS_endrun, LIS_abort,&
       LIS_alert
  use LIS_timeMgrMod,    only : LIS_julhr_date

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: order
  integer, intent(inout) :: julhr
!
! !DESCRIPTION: 
!  This routine calls the appropriate first guess forcing reading based on
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
  character(len=10) :: yyyymmddhh
  integer :: ierr

  ! Sanity check
  if ( agrmet_struc(n)%first_guess_source .ne. 'GALWEM' .and. &
       agrmet_struc(n)%first_guess_source .ne. "GFS") then
      write(LIS_logunit,*) '[ERR] First guess source is not correctly defined.'
      call LIS_endrun
  end if

  call AGRMET_julhr_date10(julhr, yyyymmddhh)
  write(LIS_logunit,*) &
       '[INFO] Searching for NWP fields valid ',yyyymmddhh

  ! Try fetching GALWEM, if requested
  ierr = 0
  if ( agrmet_struc(n)%first_guess_source == 'GALWEM' ) then
     call AGRMET_fldbld_galwem(n,order,julhr,rc)
     ierr = rc
  end if

  ! Try fetching GFS, if requested, or if GALWEM is not available.
  if ( agrmet_struc(n)%first_guess_source == "GFS" .or. &
       ierr .ne. 0) then
     if (ierr .ne. 0) then
        write(LIS_logunit,*)'[WARN] Unable to find GALWEM data!'
        write(LIS_logunit,*)'[WARN] Rolling back to GFS...'
     end if

     call AGRMET_fldbld_gfs(n,order,julhr,rc)
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
        message(2) = '  Routine:  AGRMET_fldbld.'
        message(3) = '  GALWEM and GFS GRIB data not available for '//&
             yyyymmddhh
        if(LIS_masterproc) then
           call LIS_alert( 'LIS.AGRMET_fldbld.', 1, &
                message )
           call LIS_abort( message)
        endif
     end if
  end if
  
end subroutine AGRMET_fldbld
