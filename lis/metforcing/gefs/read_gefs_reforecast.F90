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
! !ROUTINE: read_gefs_reforecast
!  \label{read_gefs_reforecast}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
!
! !INTERFACE:
subroutine read_gefs_reforecast(n, m, findex, order, filename, varname, ferror)

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_metforcingMod,  only : LIS_forc
  use gefs_forcingMod,    only : gefs_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)       :: findex
  integer, intent(in)       :: n
  integer, intent(in)       :: m
  integer, intent(in)       :: order
  character(len=*), intent(in) :: filename
  character*10, intent(in)  :: varname
  integer, intent(out)      :: ferror
!
! !DESCRIPTION:
!  For the given time, reads the forcing data from the 
!  GEFS file, transforms into LIS forcing 
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
!  \item[filename]
!    name of the file to be read
!  \item[ferror]
!    return error flag (0-fail, 1-success)
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gefs](\ref{interp_gefs}) \newline
!    Performs spatial interpolation of GEFS forecast data to the LIS grid
!  \end{description}

!EOP
  integer     :: ftn, rc
  logical     :: file_exists  
  integer     :: c,r,i

  ! Grib specific dimenstions:
  integer     :: igrib
  integer     :: nfcsttimes
  integer     :: numpts
  integer     :: grib_gridsize

  ! Eccodes - Keys:
  integer     :: stepRange
  ! For more Grib Eccodes keys, see web document:
  !  https://confluence.ecmwf.int/download/attachments/97363968/eccodes-keys-2018.pdf?api=v2
  !  Slides: 7 to 20, forecast keys on slide 14

  real, allocatable  :: gefs_grib_data(:)              ! Read-in data
  real        :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n)) ! Interp field
  logical     :: pcp_flag                              ! Precip flag for spatial interp
  character*50  :: shortName                           ! variable for 0.25 degree wind data

! ______________________________________________________________________________

  varfield = 0 
  ferror = 1
  numpts = gefs_struc(n)%nc*gefs_struc(n)%nr

  pcp_flag = .false.
  if( trim(varname) == "apcp_sfc" ) then
    pcp_flag = .true.
  endif

#if(defined USE_GRIBAPI)

! Check if file exists: ----------------------------------

  inquire( file=filename, exist=file_exists )
  if (file_exists) then      

     call grib_open_file(ftn,trim(filename),'r',rc)
     call LIS_verify(rc,'grib_open_file error in read_gefs_reforecast')

     ! Read in first two book-ends data values:
     if( LIS_rc%tscount(n) == 1 ) then  
       if( order == 1 ) then  ! Bookend-1
         call grib_new_from_file(ftn, igrib, rc)
         call LIS_verify(rc,'error in grib_new_from_file for data1 in read_gefs_reforecast')

       elseif( order == 2 ) then  ! Bookend-2
         do i = 1, 3   ! 6-hour increment
!         do i = 1, 2   ! 3-hour increment
           call grib_new_from_file(ftn, igrib, rc)
           call LIS_verify(rc,'error in grib_new_from_file for data2 in read_gefs_reforecast')
         enddo
       endif

       call grib_get_size(igrib,'values',grib_gridsize,rc)
       if( grib_gridsize /= numpts ) then
          write(LIS_logunit,*) &
             '[ERR] Number of values does not match expected', trim(filename)
          call LIS_endrun
       endif

       allocate( gefs_grib_data(grib_gridsize) )
       gefs_grib_data = LIS_rc%udef

       call grib_get(igrib,'values',gefs_grib_data,rc)
       call LIS_verify(rc,' grib_get error: data values in read_gefs_reforecast')
 
!       write(*,*) 'Min GEFS value: ', minval(gefs_grib_data)
!       write(*,*) 'Max GEFS value: ', maxval(gefs_grib_data)

       ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
       call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                                  numpts, gefs_grib_data )

       ! Spatially interp GEFS forcing field to LIS domain:
       call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
       deallocate( gefs_grib_data )

       call grib_release(igrib,rc)
       call LIS_verify(rc,'error in grib_release in read_gefs_reforecast')
       call grib_close_file(ftn)

    ! Loop over all other forecast times in files:
     else
       call grib_count_in_file(ftn, nfcsttimes, rc)
       call LIS_verify(rc,'grib_count_in_file error in read_gefs_reforecast')
 
       do i=1,nfcsttimes

         call grib_new_from_file(ftn, igrib, rc)
         call LIS_verify(rc,'error in grib_new_from_file in read_gefs_reforecast')

         call grib_get(igrib,'stepRange',stepRange,rc)
         call LIS_verify(rc,' grib_get error: stepRange in read_gefs_reforecast')
!         print *, 'grib records: ', i, igrib, stepRange, gefs_struc(n)%fcst_hour

         call grib_get(igrib,'shortName',shortName,rc)
         call LIS_verify(rc,'error in grib_get: shortName in read_gefs_reforecast')
         ! Check if local forecast hour exceeds max grib file forecast hour:
         if( gefs_struc(n)%fcst_hour > 192 ) then   
            write(LIS_logunit,*) &
               "[INFO] GEFS Forecast hour has exceeded the grib file's final"
            write(LIS_logunit,*) &
               '  forecast hour (record). Future 8-16 day forecast files and records'
            write(LIS_logunit,*) &
               ' will be added.  Run will end here for now ... '
            call LIS_endrun
         endif  
         ! Will extend forecast window to day-16 (hour 384) in future 

         if( gefs_struc(n)%fcst_hour == stepRange ) then

           call grib_get_size(igrib,'values',grib_gridsize,rc)
           if( grib_gridsize /= numpts ) then
              write(LIS_logunit,*) &
                '[ERR] Number of values does not match expected', trim(filename)
              call LIS_endrun
           endif

           allocate( gefs_grib_data(grib_gridsize) )
           gefs_grib_data = LIS_rc%udef

           call grib_get(igrib,'values',gefs_grib_data,rc)
           call LIS_verify(rc,' grib_get error: data values in read_gefs_reforecast')

!           write(*,*) 'Min GEFS value: ', minval(gefs_grib_data)
!           write(*,*) 'Max GEFS value: ', maxval(gefs_grib_data)

           ! Convert 1D GEFS field to 2D and shift 0-360 to -180 to 180 grid:
           call gefs_shift_longitude( gefs_struc(n)%nc, gefs_struc(n)%nr, &
                                      numpts, gefs_grib_data )
      
           ! Spatially interp GEFS forcing field to LIS domain:
           call interp_gefs(n, findex, gefs_grib_data, pcp_flag, varfield )
           deallocate( gefs_grib_data )

           call grib_release(igrib,rc)
           call LIS_verify(rc,'error in grib_release in read_gefs_reforecast')

           call grib_close_file(ftn)
           exit
        else
           ! Release grib message, if not used 
           call grib_release(igrib,rc)
           call LIS_verify(rc,'error in grib_release in read_gefs_reforecast')
           cycle
        endif
      enddo
    endif

! _______________________________________________________

    ! Assign GEFS forcing fields to LIS forcing arrays

     ! Air temp:
     if( trim(varname) == "tmp_2m" ) then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
             if(LIS_domain(n)%gindex(c,r).ne.-1) then 
               if(order.eq.1) then 
                  gefs_struc(n)%metdata1(1,m,&
                       LIS_domain(n)%gindex(c,r)) &
                       = varfield(c,r)
               elseif(order.eq.2) then 
                  gefs_struc(n)%metdata2(1,m,&
                       LIS_domain(n)%gindex(c,r)) &
                       = varfield(c,r)
               endif
             endif
          end do
       enddo
     endif

     ! Specific humidity
     if( trim(varname) == "spfh_2m" ) then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
             if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gefs_struc(n)%metdata1(2,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gefs_struc(n)%metdata2(2,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
             endif
          end do
       enddo
     endif

     ! Downward SW radiation:
     if( trim(varname) == "dswrf_sfc" ) then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gefs_struc(n)%metdata1(3,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gefs_struc(n)%metdata2(3,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
             endif
          end do
       enddo
     endif

     ! Downward LW radiation:
     if( trim(varname) == "dlwrf_sfc" ) then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gefs_struc(n)%metdata1(4,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gefs_struc(n)%metdata2(4,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
            endif
          end do
       enddo
     endif

     if( trim(varname) == "ugrd_10m" .or. shortName.eq."10u") then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gefs_struc(n)%metdata1(5,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gefs_struc(n)%metdata2(5,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
            endif
          end do
       enddo
     endif

     if( trim(varname) == "vgrd_10m" .or. shortName.eq."10v") then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then
              if(order.eq.1) then 
                 gefs_struc(n)%metdata1(6,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gefs_struc(n)%metdata2(6,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
            endif
          end do
       enddo
     endif

     ! Surface pressure:
     if( trim(varname) == "pres_sfc" ) then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gefs_struc(n)%metdata1(7,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gefs_struc(n)%metdata2(7,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
            endif
          end do
       enddo
     endif

     ! Accumulated total precipitation:
     if( trim(varname) == "apcp_sfc" ) then
       do r=1,LIS_rc%lnr(n)
          do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gefs_struc(n)%metdata1(8,m,&
                      LIS_domain(n)%gindex(c,r)) &
                      = varfield(c,r)
              elseif(order.eq.2) then 
                 gefs_struc(n)%metdata2(8,m,&
                      LIS_domain(n)%gindex(c,r))&
                      = varfield(c,r)
              endif
            endif
          end do
       enddo
     endif

! Snowfall is gefs_struc(n)%metdata2(9,m,LIS_domain(n)%gindex(c,r))
! Convective precip is gefs_struc(n)%metdata2(10,m,LIS_domain(n)%gindex(c,r))

! GEFS file not found:
  else
     write(LIS_logunit,*) &
       '[ERR] Could not find file: ',trim(filename)
     ferror = 0
  endif
#endif

end subroutine read_gefs_reforecast


!BOP
! !ROUTINE: gefs_shift_longitude
! \label{gefs_shift_longitude}
!
! !INTERFACE:
subroutine gefs_shift_longitude( nc, nr, numpts, indata1d )

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: nc, nr, numpts
  real, intent(inout)    :: indata1d(numpts) 
!
! !DESCRIPTION:
!  This subroutine converts 1D GEFS field to 2D and 
!  shifts the Reforecast data from 0-360 to -180 to 180 grid
! 
!  The arguments are:
!  \begin{description}
!  \item[nc]
!   Number of input data columns
!  \item[nr]
!   Number of input data rows
!  \item[numpts]
!   Total number of input data points
!  \item[indata1d]
!   Input data field to be shifted
!  \end{description}
!
!EOP
  integer  :: i, c, r
  real     :: indata2d(nc,nr)
  real     :: indata2d_shift(nc,nr)

  indata2d = -9999.
  indata2d_shift = -9999.

  ! Convert 1d grid to 2d:
  i = 0
  do r = 1, nr
     do c = 1, nc
        i = i + 1
        indata2d(c,r) = indata1d(i)
     enddo
  enddo

  ! Shift longitudinal grid
  do r = 1, nr
     do c = 1, (nc/2)
        indata2d_shift(c,r) = indata2d((nc/2+c),r)
     enddo
     do c = (nc/2+1), nc
        indata2d_shift(c,r) = indata2d((c-nc/2),r)
     enddo
  enddo

  ! Convert back to 1D field 
  i = 0
  do r = 1, nr
     do c = 1, nc
        i = i + 1
        indata1d(i) = indata2d_shift(c,r)
     enddo
  enddo

end subroutine gefs_shift_longitude

