!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_galwemge
! \label{read_galwemge}
!
! !REVISION HISTORY:
! 09 May 2022: Yeosang Yoon, initial code
!
! !INTERFACE:
subroutine read_galwemge(n, m, findex, order, gribfile, rc)

! !USES:
  use LIS_coreMod,       only : LIS_rc
  use LIS_logMod
  use galwemge_forcingMod, only : galwemge_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none

!ARGUMENTS:
  integer,           intent(in)    :: n
  integer,           intent(in)    :: m       ! number of ensembles
  integer,           intent(in)    :: findex  ! forcing index
  integer,           intent(in)    :: order
  character(len=*),  intent(in)    :: gribfile

!DESCRIPTION:
!  For the given time, reads the forcing data from the
!  GALWEM-GE file, transforms into LIS forcing
!  parameters and interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    hourly instance, order=2, read the next hourly instance)
!  \item[gribfile]
!    name of the file to be read
!  \item[rc]
!    return error flag (0-fail, 1-success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[interp\_galwemge](\ref{interp_galwemge}) \newline
!    Performs spatial interpolation of GALEM-GE forecast data to the LIS grid
!  \end{description}

!EOP
  integer         :: ftn, igrib, ierr
  integer         :: center
  character*100   :: gtype
  integer         :: file_julhr
  integer         :: yr1, mo1, da1, hr1
  character*255   :: message     ( 20 )
  integer         :: iginfo      ( 40 )
  real            :: gridres_dlat, gridres_dlon
  integer         :: ifguess, jfguess
  integer         :: dataDate, dataTime
  integer         :: fc_hr

  real            :: tair(LIS_rc%lnc(n),LIS_rc%lnr(n))      !Temperature interpolated to 2 metres [K] 
  real            :: qair(LIS_rc%lnc(n),LIS_rc%lnr(n))      !Relative humidity interpolated to 2 metres[kg/kg] 
  real            :: swdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward shortwave flux at the ground [W/m^2] 
  real            :: lwdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward longwave radiation at the ground [W/m^2] 
  real            :: uwind(LIS_rc%lnc(n),LIS_rc%lnr(n))     !Instantaneous zonal wind interpolated to 10 metres [m/s]
  real            :: vwind(LIS_rc%lnc(n),LIS_rc%lnr(n))     !Instantaneous meridional wind interpolated to 10 metres[m/s]
  real            :: ps(LIS_rc%lnc(n),LIS_rc%lnr(n))        !Instantaneous Surface Pressure [Pa] 
  real            :: prectot(LIS_rc%lnc(n),LIS_rc%lnr(n))   !Accumulated total precipitation [kg/m^2] 
 
  integer         :: r, c 
  real            :: e, es
  real            :: q(LIS_rc%lnc(n),LIS_rc%lnr(n))         !Specific humidity interpolated to 2 metres[kg/kg]

  integer, intent(out) :: rc

  ! Initialize return code to "no error".  We will change it below if
  ! necessary.
  rc = 0

#if (defined USE_GRIBAPI)

   call grib_open_file(ftn,trim(gribfile),'r',ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] Failed to open - '//trim(gribfile)
      call LIS_endrun()
   end if

   ! Read in the first grib record, unpack the header and extract
   ! section 1 and section 2 information.
   call grib_new_from_file(ftn,igrib,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] failed to read - '//trim(gribfile)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'centre',center,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grib_get: ' // &
           'centre in read_galwemge'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'gridType',gtype,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'gridtype in read_galwemge'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   if(trim(gtype).ne."regular_ll") then
      write(LIS_logunit,*)'[ERR] GALWEM-GE file not on lat/lon grid!'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'Ni',iginfo(1),ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: Ni read_galwemge'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'Nj',iginfo(2),ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: Nj in read_galwemge'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'jDirectionIncrementInDegrees',gridres_dlat,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'jDirectionIncrementInDegrees ' // &
           'in read_galwemge'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   ! EMK...Added dlon
   call grib_get(igrib, 'iDirectionIncrementInDegrees', gridres_dlon, ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'iDirectionIncrementInDegrees ' // &
           'in read_galwemge'
      call grib_release(igrib, ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'dataDate',dataDate,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'dataDate in read_galwemge'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   call grib_get(igrib,'dataTime',dataTime,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'dataTime in read_galwemge'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      call LIS_endrun()
   endif

   ! Here we tentatively have a file we can use. Close it for now, and
   ! prepare to pull the appropriate variables.
   call grib_release(igrib,ierr)
   call grib_close_file(ftn)

   ifguess = iginfo(1)
   jfguess = iginfo(2)

   call fldbld_read_galwemge(n, findex, order, gribfile, ifguess, jfguess, &
           tair, qair, swdown, lwdown, uwind, vwind, ps, prectot, rc)

   call assign_processed_galwemgeforc(n, m, order, 1, tair)

!  convert relative humidity to specific humidity using Bolton (1980)
!
!   e = 6.112*exp((17.67*Td)/(Td + 243.5));
!   q = (0.622 * e)/(p - (0.378 * e));
!    where:
!       e = vapor pressure in mb;
!       Td = dew point in deg C;
!       p = surface pressure in mb;
!       q = specific humidity in kg/kg.
!
!   Td dew point temperature is
!   es = 6.112 * exp((17.67 * T)/(T + 243.5));
!   e = es * (RH/100.0);
!   Td = log(e/6.112)*243.5/(17.67-log(e/6.112));
!     where:
!       T = temperature in deg C;
!       es = saturation vapor pressure in mb;
!       e = vapor pressure in mb;
!       RH = Relative Humidity in percent;
!       Td = dew point in deg C

   do r=1,LIS_rc%lnr(n)
      do c=1,LIS_rc%lnc(n)
         tair(c,r) =tair(c,r)-273.15 !to celcius 
         es = 6.112*exp((17.67+tair(c,r))/(tair(c,r)+243.5))
         e = es*(qair(c,r)/100.0)
         q(c,r) = (0.622*e)/(ps(c,r)*0.01-0.378*e)   ! 1mb=100Pa        
      enddo
   enddo
   call assign_processed_galwemgeforc(n, m, order, 2, q)
   call assign_processed_galwemgeforc(n, m, order, 3, swdown)
   call assign_processed_galwemgeforc(n, m, order, 4, lwdown)
   call assign_processed_galwemgeforc(n, m, order, 5, uwind)
   call assign_processed_galwemgeforc(n, m, order, 6, vwind)
   call assign_processed_galwemgeforc(n, m, order, 7, ps)
   call assign_processed_galwemgeforc(n, m, order, 8, prectot)
      
#endif

end subroutine read_galwemge

subroutine fldbld_read_galwemge(n, findex, order, gribfile, ifguess, jfguess,     &
                                     tair, qair, swdown, lwdown,                &
                                     uwind, vwind, ps, prectot, rc)                            
 
! !USES:
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,        intent(in)    :: n
  integer, intent(in)           :: findex  ! Forcing index
  integer,        intent(in)    :: order
  character(len=*),  intent(in) :: gribfile
  integer,        intent(in)    :: ifguess
  integer,        intent(in)    :: jfguess
  real,           intent(out)   :: tair(LIS_rc%lnc(n),LIS_rc%lnr(n))      !Temperature interpolated to 2 metres [K]
  real,           intent(out)   :: qair(LIS_rc%lnc(n),LIS_rc%lnr(n))      !Relative humidity interpolated to 2 metres[kg/kg]
  real,           intent(out)   :: swdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward shortwave flux at the ground [W/m^2]
  real,           intent(out)   :: lwdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward longwave radiation at the ground [W/m^2]
  real,           intent(out)   :: uwind(LIS_rc%lnc(n),LIS_rc%lnr(n))     !Instantaneous zonal wind interpolated to 10 metres [m/s]
  real,           intent(out)   :: vwind(LIS_rc%lnc(n),LIS_rc%lnr(n))     !Instantaneous meridional wind interpolated to 10 metres[m/s]
  real,           intent(out)   :: ps(LIS_rc%lnc(n),LIS_rc%lnr(n))        !Instantaneous Surface Pressure [Pa]
  real,           intent(out)   :: prectot(LIS_rc%lnc(n),LIS_rc%lnr(n))   !Accumulated Total precipitation [kg/m^2]
  integer,        intent(out)   :: rc
!
! !DESCRIPTION:
!
!     To read UK Unified Model (GALWEM-GE) data in GRIB-2 format.
!
!EOP
  character*9                   :: cstat
  character*255                 :: message     ( 20 )
  character(len=7)              :: grib_msg
  character(len=7)              :: check_galwemge_message
  integer                       :: count_tair, count_qair
  integer                       :: count_swdown, count_lwdown
  integer                       :: count_uwind, count_vwind
  integer                       :: count_ps, count_prectot
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: igrib
  integer                       :: ftn
  integer                       :: kk, nvars
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, surface_val, level_val
  real, allocatable             :: dum1d   ( : )

  real, allocatable             :: fg_tair    ( : , : )
  real, allocatable             :: fg_qair    ( : , : )
  real, allocatable             :: fg_swdown  ( : , : )
  real, allocatable             :: fg_lwdown  ( : , : )
  real, allocatable             :: fg_uwind   ( : , : )
  real, allocatable             :: fg_vwind   ( : , : )
  real, allocatable             :: fg_ps      ( : , : )
  real, allocatable             :: fg_prectot ( : , : )

  logical                       :: found_inq

  rc = 1 ! Initialize as "no error"

  ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
  ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
  ! writing error messages to stdout/stderr, which may lead to runtime
  ! problems.
  inquire(file=trim(gribfile),exist=found_inq)
  if (.not. found_inq) then
     write(LIS_logunit,*) '[WARN] Cannot find file '//trim(gribfile)
     rc = 0
     return
  end if

#if (defined USE_GRIBAPI)

  ! If a problem occurs here, we can just return immediately since no
  ! memory has been allocated yet.
  call grib_open_file(ftn,trim(gribfile),'r',ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] Failed to open - '//trim(gribfile)
     rc = 0
     return
  end if

  allocate ( fg_tair    (ifguess, jfguess) )
  allocate ( fg_qair    (ifguess, jfguess) )
  allocate ( fg_swdown  (ifguess, jfguess) )
  allocate ( fg_lwdown  (ifguess, jfguess) )
  allocate ( fg_uwind   (ifguess, jfguess) )
  allocate ( fg_vwind   (ifguess, jfguess) )
  allocate ( fg_ps      (ifguess, jfguess) )
  allocate ( fg_prectot (ifguess, jfguess) )

  allocate ( dum1d   (ifguess*jfguess) )

  ! Initalization
  fg_tair = LIS_rc%udef
  fg_qair = LIS_rc%udef
  fg_swdown = LIS_rc%udef
  fg_lwdown = LIS_rc%udef
  fg_uwind = LIS_rc%udef
  fg_vwind = LIS_rc%udef
  fg_ps = LIS_rc%udef
  fg_prectot = LIS_rc%udef 
  dum1d = LIS_rc%udef

  tair = LIS_rc%udef
  qair = LIS_rc%udef
  swdown = LIS_rc%udef
  lwdown = LIS_rc%udef
  uwind = LIS_rc%udef
  vwind = LIS_rc%udef
  ps = LIS_rc%udef
  prectot = LIS_rc%udef

  ! From this point, we must deallocate memory before returning.
  ! Unfortunately this means using a GOTO statement if a problem is
  ! encountered, but such is life.
  count_tair  = 0
  count_qair   = 0
  count_swdown  = 0
  count_lwdown = 0
  count_uwind = 0
  count_vwind = 0
  count_ps = 0
  count_prectot = 0

  call grib_count_in_file(ftn,nvars,ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] in grib_count_in_file in ' // &
          'fldbld_read_galwemge'
     goto 100
  end if

  ! Tentatively loop through every field in GRIB file looking for the variables
  ! we want. The code below will exit the loop early if a problem is found *or*
  ! once all the required variables are found and read in.
  do kk=1,nvars

     call grib_new_from_file(ftn,igrib,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] failed to read - '//trim(gribfile)
        goto 100
     end if

     call grib_get(igrib,'discipline',param_disc_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grib_get: parameterNumber in ' // &
             'fldbld_read_galwemge'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterCategory in ' // &
             'fldbld_read_galwemge'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'parameterNumber',param_num_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterNumber in ' // &
             'fldbld_read_galwemge'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'typeOfFirstFixedSurface',surface_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'fldbld_read_galwemge'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'level',level_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'fldbld_read_galwemge'
        call grib_release(igrib,ierr)
        goto 100
     end if

     ! We have enough information to determine what GRIB parameter this is.
     grib_msg = check_galwemge_message(param_disc_val, &
          param_cat_val, param_num_val, surface_val, level_val)

     ! Skip this field if GRIB parameter is not required.
     if (grib_msg == 'none') then
        call grib_release(igrib,ierr)
        if (ierr .ne. 0) then
           write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
                'fldbld_read_galwemge'
           goto 100
        end if
        cycle ! Not a message we are interested in.
     end if

     call grib_get(igrib,'values',dum1d,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: values in ' // &
             'fldbld_read_galwemge'
        call grib_release(igrib,ierr)
        goto 100
     end if

     select case (grib_msg)
     case('t2')       ! 2-m temperature
        fg_tair = reshape(dum1d, (/ifguess,jfguess/))
        count_tair = count_tair + 1
     case ('q2')      ! 2-m relative humidity 
        fg_qair = reshape(dum1d, (/ifguess,jfguess/))
        count_qair = count_qair + 1
     case ('swdown')  ! downward shortwave
        fg_swdown = reshape(dum1d, (/ifguess,jfguess/))
        count_swdown = count_swdown + 1
     case ('lwdown')  ! downward longwave radiation
        fg_lwdown = reshape(dum1d, (/ifguess,jfguess/))
        count_lwdown = count_lwdown + 1
     case ('u10')     ! 10m u-wind
        fg_uwind = reshape(dum1d, (/ifguess,jfguess/))
        count_uwind = count_uwind + 1
     case('v10')      ! 10m v-wind
        fg_vwind = reshape(dum1d, (/ifguess,jfguess/))
        count_vwind = count_vwind + 1
     case('ps')       ! surface pressure
        fg_ps = reshape(dum1d, (/ifguess,jfguess/))
        count_ps = count_ps + 1
     case('prectot')  ! accumulated total precipitation
        fg_prectot = reshape(dum1d, (/ifguess,jfguess/))
        count_prectot = count_prectot + 1

     case default ! Internal error, we shouldn't be here
        write(LIS_logunit,*)'[WARN] Unknown grib_message ',grib_msg
        write(LIS_logunit,*)'Aborting...'
        flush(LIS_logunit)
        write(cstat,'(i9)',iostat=istat1) ierr
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  fldbld_read_galwemge.'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(gribfile)
        if( istat1 .eq. 0 )then
           message(5) = '  Status = ' // trim(cstat)
        endif
        call LIS_abort( message)
     end select

     ! Finished with this field
     call grib_release(igrib,ierr)
     if (ierr .ne. 0) then
        write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
             'fldbld_read_galwemge'
        goto 100
     end if

     ! Jump out of loop early if we have everything
     if ( (count_tair   .eq. 1)  .and. (count_qair    .eq. 1)   .and. &
          (count_swdown .eq. 1 ) .and. (count_lwdown  .eq. 1)   .and. &
          (count_uwind  .eq. 1)  .and. (count_vwind   .eq. 1)   .and. &
          (count_ps     .eq. 1)  .and. (count_prectot .eq. 1) ) then
        exit
     end if

  enddo ! Loop through all GRIB file fields

  ! Make sure we have all required fields.
  !if ( (count_tair   .ne. 1)  .or. (count_qair    .ne. 1)   .or. &
  !     (count_swdown .ne. 1 ) .or. (count_lwdown  .ne. 1)   .or. &
  !     (count_uwind  .ne. 1)  .or. (count_vwind   .ne. 1)   .or. &
  !     (count_ps     .ne. 1)  .or. (count_prectot .ne. 1) ) then
  !   write(LIS_logunit,*)'[WARN] Missing data from GALWEM-GE GRIB file!'
  !   !goto 100
  !end if

  ! Interpolate the fields to the LIS grid
  ! tair
  call interp_galwemge(n, findex, ifguess, jfguess, .false., fg_tair, tair)
  ! qair
  call interp_galwemge(n, findex, ifguess, jfguess, .false., fg_qair, qair)
  ! swdown
  call interp_galwemge(n, findex, ifguess, jfguess, .false., fg_swdown, swdown)
  ! lwdown
  call interp_galwemge(n, findex, ifguess, jfguess, .false., fg_lwdown, lwdown)
  ! uwind
  call interp_galwemge(n, findex, ifguess, jfguess, .false., fg_uwind, uwind)
  ! vwind
  call interp_galwemge(n, findex, ifguess, jfguess, .false., fg_vwind, vwind)
  ! ps
  call interp_galwemge(n, findex, ifguess, jfguess, .false., fg_ps, ps)
  ! prectot
  call interp_galwemge(n, findex, ifguess, jfguess, .true., fg_prectot, prectot)

  ! At this point, we have everything.  Close the file and return.
  call grib_close_file(ftn)
  rc = 1
  return

  ! Jump down here to clean up memory before returning after finding a
  ! problem.
  100 continue
  call grib_close_file(ftn)

  deallocate ( dum1d      )
  deallocate ( fg_tair    )
  deallocate ( fg_qair    )
  deallocate ( fg_swdown  )
  deallocate ( fg_lwdown  )
  deallocate ( fg_uwind   )
  deallocate ( fg_vwind   )
  deallocate ( fg_ps      )
  deallocate ( fg_prectot )
  rc = 0
#endif

end subroutine fldbld_read_galwemge

function check_galwemge_message(param_disc_val, &
     param_cat_val, param_num_val, surface_val, level_val)
! !USES:
! none

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: param_disc_val, param_cat_val, &
        param_num_val, surface_val, level_val
   character(len=7)    :: check_galwemge_message
!EOP

   if     ( param_disc_val == 0 .and. &
            param_cat_val  == 0 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 103 .and. &
            level_val == 2) then
      check_galwemge_message = 't2'      ! Temperature interpolated to 2 metres [K] 
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 1 .and. &
            param_num_val  == 1 .and. &
            surface_val    == 103 .and. &
            level_val == 2) then
      check_galwemge_message = 'q2'      ! Relative humidity interpolated to 2 metres [kg/kg]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 4 .and. &
            param_num_val  == 7 .and. &
            surface_val    == 1 ) then
      check_galwemge_message = 'swdown'  ! Downward shortwave flux at the ground [W/m^2]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 5 .and. &
            param_num_val  == 3 .and. &
            surface_val    == 1 ) then
      check_galwemge_message = 'lwdown'  ! Downward longwave radiation at the ground [W/m^2] 
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 2 .and. &
            surface_val    == 103 .and. &
            level_val == 10) then
      check_galwemge_message = 'u10'     ! Instantaneous zonal wind interpolated to 10 metres [m/s] 
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 3 .and. &
            surface_val    == 103 .and. &
            level_val == 10) then
      check_galwemge_message = 'v10'     ! Instantaneous meridional wind interpolated to 10 metres[m/s]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 3 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 1 ) then
      check_galwemge_message = 'ps'      ! Instantaneous Surface Pressure [Pa]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 1 .and. &
            param_num_val  == 8 .and. &
            surface_val    == 1 ) then
      check_galwemge_message = 'prectot' ! Accumulated total precipitation [kg/m2]
   else
      check_galwemge_message = 'none'
   endif
end function check_galwemge_message

subroutine interp_galwemge(n, findex, ifguess, jfguess, pcp_flag, input, output)

! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use galwemge_forcingMod, only : galwemge_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)   :: findex
  integer, intent(in)  :: ifguess
  integer, intent(in)  :: jfguess
  logical, intent(in)  :: pcp_flag
  real,    intent(in)  :: input  ( ifguess,jfguess )
  real,    intent(out) :: output ( LIS_rc%lnc(n),LIS_rc%lnr(n) )
!
! !DESCRIPTION:
!
! This routine interpolates the GALWEM-GE data to the LIS grid.

!EOP

  integer   :: mi, mo
  integer   :: k
  integer   :: i,j
  integer   :: iret
  integer   :: midway
  character(len=50) :: method
  real, allocatable, dimension(:,:)    :: var
  logical*1, allocatable, dimension(:) :: lb
  logical*1, allocatable, dimension(:) :: lo

  allocate(var(ifguess,jfguess))
  allocate(lb(ifguess*jfguess))
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  mi = ifguess * jfguess
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  lb = .true.
  lo = .true.
  output = LIS_rc%udef

  ! translate from (0,360) to (-180,180).
  midway = ifguess/2
  do j = 1, jfguess
     do i = 1, ifguess
        if ( (i+midway) < ifguess ) then
           var(i,j) = input((i+midway),j)
        else
           var(i,j) = input((i-midway+1),j)
        endif
     enddo
  enddo

  ! Interpolate to LIS grid
  select case( LIS_rc%met_interp(findex) )

    case( "bilinear" )
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                     &
          var,lo,output,mi,mo,                                         &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          galwemge_struc(n)%w111,galwemge_struc(n)%w121,                   &
          galwemge_struc(n)%w211,galwemge_struc(n)%w221,                   &
          galwemge_struc(n)%n111,galwemge_struc(n)%n121,                   &
          galwemge_struc(n)%n211,galwemge_struc(n)%n221,LIS_rc%udef,iret)

    case( "budget-bilinear" )
     if (pcp_flag) then
        call conserv_interp(LIS_rc%gridDesc(n,:),lb,                   &
             var,lo,output,mi,mo,                                      &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
             galwemge_struc(n)%w112,galwemge_struc(n)%w122,                &
             galwemge_struc(n)%w212,galwemge_struc(n)%w222,                &
             galwemge_struc(n)%n112,galwemge_struc(n)%n122,                &
             galwemge_struc(n)%n212,galwemge_struc(n)%n222,LIS_rc%udef,iret)
     else
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                  &
             var,lo,output,mi,mo,                                      &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
             galwemge_struc(n)%w111,galwemge_struc(n)%w121,                &
             galwemge_struc(n)%w211,galwemge_struc(n)%w221,                &
             galwemge_struc(n)%n111,galwemge_struc(n)%n121,                &
             galwemge_struc(n)%n211,galwemge_struc(n)%n221,LIS_rc%udef,iret)
     endif

     case( "neighbor" )
        call neighbor_interp(LIS_rc%gridDesc(n,:),lb,                     &
             var,lo,output,mi,mo,                                         &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                        &
             galwemge_struc(n)%n113,LIS_rc%udef,iret)

     case DEFAULT
        write(LIS_logunit,*) 'ERR: Unexpected interpolation method'
        write(LIS_logunit,*) '     in interp_galwemge_first_guess'
        write(LIS_logunit,*) '     ', trim(LIS_rc%met_interp(findex))
        call LIS_endrun()  
  end select

  deallocate(var)
  deallocate(lb)
  deallocate(lo)

end subroutine interp_galwemge

!BOP
!
! !ROUTINE: assign_processed_galwemgeforc
! \label{assign_processed_galwemgeforc}
!
! !INTERFACE:
subroutine assign_processed_galwemgeforc(n,m,order,var_index,galwemgeforc)
! !USES:
  use LIS_coreMod
  use galwemge_forcingMod, only : galwemge_struc
!
! !DESCRIPTION:
!  This routine assigns the interpolated GALWEM-GE forcing data
!  to the module data structures to be used later for
!  time interpolation
!
!EOP
  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: m 
  integer, intent(in) :: order
  integer, intent(in) :: var_index
  real,    intent(in) :: galwemgeforc(LIS_rc%lnc(n),LIS_rc%lnr(n))

  integer :: c,r

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then
           if(order.eq.1) then
              galwemge_struc(n)%metdata1(var_index,m,&
                   LIS_domain(n)%gindex(c,r)) = &
                   galwemgeforc(c,r)
           elseif(order.eq.2) then
              galwemge_struc(n)%metdata2(var_index,m,&
                   LIS_domain(n)%gindex(c,r)) = &
                   galwemgeforc(c,r)
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_galwemgeforc

