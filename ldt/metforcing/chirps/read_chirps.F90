!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_chirps 
!  \label{read_chirps}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Kato; Initial Specification adopted from LDT/retberg.F90
!  15 Jul 2015: K. Arsenault Hall;  Adapted for CHIRPS
!  29 Oct 2025: J Nattala; Updated to support both CHIRPS v2.0 and v3.0;
!                          use generic chirps_struc
!   9 Mar 2026: J Nattala; Fixed udef propagation through metdata2 assignment;
!                          fixed lb initialization via where construct
!  27 Mar 2026: J Nattala; Added excl_lat_eps
!
! !INTERFACE:
subroutine read_chirps( n, findex, chirps_filename, year, mon, day, ferror )

! !USES:
  use LDT_coreMod,         only : LDT_rc, LDT_domain, LDT_localPet
  use LDT_metforcingMod,   only : LDT_forc
  use LDT_timeMgrMod,      only : LDT_tick
  use LDT_logMod,          only : LDT_logunit, LDT_endrun, LDT_verify
  use chirps_forcingMod,   only : chirps_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n       ! nest
  integer, intent(in)    :: findex
  integer, intent(in)    :: year
  integer, intent(in)    :: mon
  integer, intent(in)    :: day
  integer, intent(inout) :: ferror  ! set to non-zero if there's an error
  character(len=*), intent(in) :: chirps_filename

! !DESCRIPTION:
!  For the given time, reads the CHIRPS precipitation data, 
!   then spatially transforms the input fields to the
!   the given LDT domain.
!
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[chirps_filename]
!    CHIRPS filename to be opened and read
!  \item[ferror]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP

  integer     :: ios     ! Input/output status
  integer     :: c, i, j, index1, &
                 years_since, years_since_lp, &
                 iyr, imon, doy, total_days
  integer     :: days(12), days_lp(12)
  integer     :: iret, mo, actual_row
  real        :: lat_approx
  real,    parameter :: excl_lat_eps = 1.0e-4
  data days /31,28,31,30,31,30,31,31,30,31,30,31/
  data days_lp /31,29,31,30,31,30,31,31,30,31,30,31/
  character(3) :: cdoy

  ! Netcdf inputs:
  integer    :: nid     ! Netcdf file unit ID 
  integer    :: varid   ! Netcdf file id
  real, allocatable :: precip_nc(:,:)

  ! Input CHIRPS grid:
  real,      allocatable :: chirpsprec_1d(:)
  logical*1, allocatable :: lb(:)

  ! Output LDT grid:
  logical*1 :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real      :: lisprec_1d(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real      :: lisprec_2d(LDT_rc%lnc(n),LDT_rc%lnr(n))

  external :: conserv_interp
  external :: neighbor_interp
  external :: upscaleByAveraging
  
! ______________________________________________

  lisprec_2d = LDT_rc%udef

! --- Calculate day number since 1980-1-1 00Z ---

! Calculate the number of years since 1980:
  years_since    = 0
  years_since_lp = 0
  do iyr = 1980, (year-1)
   ! Account for leap year:
     if( (mod(iyr,4)==0 .and. mod(iyr,100).ne.0) .or. &
                 mod(iyr,400)==0 ) then
        years_since_lp = years_since_lp + 1
     else
        years_since = years_since + 1
     endif
  enddo

  ! Calculate day of year (current):
  doy = 0
  do imon = 1, (mon-1)
     ! Leap year:
     if( (mod(year,4)==0 .and. mod(year,100).ne.0) .or. &
                      mod(year,400)==0 ) then
        doy = doy + days_lp(imon)
        ! Non-leap year:
     else
        doy = doy + days(imon)
     endif
  enddo
  doy = doy + day

  ! Total day count since 1980-1-1, 0Z:
  total_days = (years_since*365) + (years_since_lp*366) + doy

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  ! --- Open NetCDF file ---
  ios = nf90_open(path=chirps_filename, &
    mode=NF90_NOWRITE, ncid=nid)
  call LDT_verify(ios, 'Error in nf90_open in read_chirps')

  !- Read in Netcdf file .... 
  allocate( precip_nc(chirps_struc(n)%nc, chirps_struc(n)%nr) )
  precip_nc = LDT_rc%udef

  ios = nf90_inq_varid(nid, "precip", varid)
  call LDT_verify(ios, 'Error in nf90_inq_varid in read_chirps')

  !- Read in subsetted domain:
  ios = nf90_get_var(nid, varid, precip_nc, &
       start=(/chirps_struc(n)%start_nc, chirps_struc(n)%start_nr, doy/), &
       count=(/chirps_struc(n)%nc, chirps_struc(n)%nr, 1/) )
  call LDT_verify(ios, 'Error in nf90_get_var in read_chirps')

  write(unit=cdoy, fmt='(i3.3)') doy
  write(LDT_logunit,*) "[INFO] Reading CHIRPS Day since 1980-1-1 (and DOY): ", &
    total_days, "("//cdoy//")"

  ios = nf90_close(nid)
  call LDT_verify(ios, 'Error in nf90_close in read_chirps')

#endif

  !- Transferring current data to 1-D array for interpolation
  write(LDT_logunit,*) "[INFO] Reprojecting the CHIRPS domain to the LIS run domain ... "

  ! --- Apply northern/southern exclusion latitudes (CHIRPS3 only) ---
  if( chirps_struc(n)%version == 3 ) then
     do j = 1, chirps_struc(n)%nr
        actual_row = chirps_struc(n)%start_nr + j - 1
        lat_approx = -60.0 + (chirps_struc(n)%yres/2.0) + &
          real(actual_row-1)*chirps_struc(n)%yres
        if( lat_approx >= chirps_struc(n)%north_excl_lat - &
          excl_lat_eps .or. lat_approx <= chirps_struc(n)%south_excl_lat &
          + excl_lat_eps ) then
           precip_nc(:,j) = LDT_rc%udef
        endif
     enddo
  endif

  c=0
  allocate( chirpsprec_1d(chirps_struc(n)%nc*chirps_struc(n)%nr) )
  allocate( lb(chirps_struc(n)%nc*chirps_struc(n)%nr) )
  chirpsprec_1d = LDT_rc%udef

  do j = 1, chirps_struc(n)%nr
     do i = 1, chirps_struc(n)%nc
        c = c + 1
        chirpsprec_1d(c) = precip_nc(i,j)
     enddo
  enddo
  deallocate( precip_nc )

  !- Interpolate if forcing and model grids are not both one deg.
  !- fix lb initialization
  lb  = .false.; lo  = .true.
  mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
  lisprec_1d = LDT_rc%udef
  where( chirpsprec_1d .ne. LDT_rc%udef .and. &
    chirpsprec_1d .ge. 0.0 ) lb = .true.

  !- Select interpolation method:
  select case( LDT_rc%met_gridtransform(findex) )

  case( "budget-bilinear" )
    call conserv_interp( LDT_rc%met_gridDesc(findex,:), &
      lb, chirpsprec_1d, lo, lisprec_1d, &
      chirps_struc(n)%mi, mo, &
      LDT_domain(n)%lat, LDT_domain(n)%lon, &
      chirps_struc(n)%w112, chirps_struc(n)%w122, &
      chirps_struc(n)%w212, chirps_struc(n)%w222, &
      chirps_struc(n)%n112, chirps_struc(n)%n122, &
      chirps_struc(n)%n212, chirps_struc(n)%n222, &
      LDT_rc%udef, iret)

  case( "neighbor" )
    call neighbor_interp( LDT_rc%met_gridDesc(findex,:), &
      lb, chirpsprec_1d, lo, lisprec_1d, &
      chirps_struc(n)%mi, mo, &
      LDT_domain(n)%lat, LDT_domain(n)%lon, &
      chirps_struc(n)%n113, LDT_rc%udef, iret)

  case( "average" )
    call upscaleByAveraging( chirps_struc(n)%mi, mo, LDT_rc%udef, &
      chirps_struc(n)%n111, lb, chirpsprec_1d, lo, lisprec_1d)

  end select

  ! Convert precip from 1D to 2D LDT output gridspace:
  lisprec_2d = LDT_rc%udef
  c = 0
  do j = 1, LDT_rc%lnr(n)
     do i = 1, LDT_rc%lnc(n)
        if( LDT_domain(n)%gindex(i,j) .ne. -1 ) then
           lisprec_2d(i,j) = lisprec_1d(i+c)
        endif
     enddo
     c = c + LDT_rc%lnc(n)
  enddo

  deallocate( lb, chirpsprec_1d )

  ! Store metdata2: fixed udef propagation:
  do j = 1, LDT_rc%lnr(n)
    do i = 1, LDT_rc%lnc(n)
      index1 = LDT_domain(n)%gindex(i,j)
      if( index1 .ne. -1 ) then
         if( lisprec_2d(i,j) >= 0. ) then
            LDT_forc(n,findex)%metdata2(1,index1) = lisprec_2d(i,j)
         else
            LDT_forc(n,findex)%metdata2(1,index1) = LDT_rc%udef
         endif
      endif
    enddo
  enddo

end subroutine read_chirps
