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
! !ROUTINE: read_chirps2 
!  \label{read_chirps2}
!
! !REVISION HISTORY:
!  26 Jan 2007: Hiroko Kato; Initial Specification adopted from LIS/retberg.F90
!  15 Jul 2015: K. Arsenault Hall;  Adapted for CHIRPS
!
! !INTERFACE:
subroutine read_chirps2( n, kk, findex, chirps_filename, year, mon, day, ferror )

! !USES:
  use LIS_coreMod,         only : LIS_rc,LIS_domain, LIS_localPet, &
                                  LIS_masterproc
  use LIS_metforcingMod,   only : LIS_forc
  use LIS_timeMgrMod,      only : LIS_tick
  use LIS_logMod,          only : LIS_logunit, LIS_endrun, LIS_verify
  use chirps2_forcingMod,  only : chirps2_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n       ! nest
  integer, intent(in)    :: kk      ! forecast member
  integer, intent(in)    :: findex
  integer, intent(in)    :: year
  integer, intent(in)    :: mon
  integer, intent(in)    :: day
  integer, intent(inout) :: ferror   ! set to non-zero if there's an error
  character(len=*), intent(in) :: chirps_filename

! !DESCRIPTION:
!  For the given time, reads the CHIRPS 2.0 precipitation data, 
!   then spatially transforms the input fields to the
!   the given LIS domain.
!
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[chirps2\_filename]
!    CHIRPS 2 filename to be opened and read
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

   integer    :: ios     ! Input/output status
   integer    :: c,i,j
   integer    :: index1

   integer    :: years_since
   integer    :: years_since_lp
   integer    :: iyr, imon
   integer    :: doy
   integer    :: total_days
   integer    :: days(12), days_lp(12)
   data days /31,28,31,30,31,30,31,31,30,31,30,31/
   data days_lp /31,29,31,30,31,30,31,31,30,31,30,31/
   character(3):: cdoy

 ! Netcdf inputs:
   integer    :: nid     ! Netcdf file unit ID 
   integer    :: varid   ! Netcdf file id
   real, allocatable :: precip_nc(:,:)
!   real       :: latitude(chirps2_struc(n)%nr)

 ! Input CHIRPS grid:
   real, allocatable :: chirpsprec_1d(:)
   logical*1, allocatable :: lb(:)

 ! Output LIS grid:
   logical*1  :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   real       :: lisprec_1d(LIS_rc%lnc(n)*LIS_rc%lnr(n))
   real       :: lisprec_2d(LIS_rc%lnc(n),LIS_rc%lnr(n))
   integer    :: iret, mo

! ______________________________________________

  lisprec_2d = LIS_rc%udef

! Calculate the day number since 1980-1-1, 00Z:

! Calculate the number of years since 1980:
  years_since = 0
  years_since_lp = 0
  do iyr = 1980, (year-1)
   ! Account for leap year:
     if( (mod(iyr,4)== 0 .and. mod(iyr,100).ne.0) .or. &
         (mod(iyr,400)==0) ) then
        years_since_lp = years_since_lp + 1
     else   ! Non-leap years
        years_since = years_since + 1
     endif
  end do

! Calculate day of year (current):
  doy = 0
  do imon = 1, (mon-1) 
   ! Leap year:
     if( (mod(year,4)== 0 .and. mod(year,100).ne.0) .or. &
         (mod(year,400)==0) ) then
        doy = doy + days_lp(imon)
   ! Non-leap year:
     else
        doy = doy + days(imon)
     endif
  enddo
  doy = doy + day

! Total day count since 1980-1-1, 0Z:
  total_days = (years_since*365)+(years_since_lp*366)+doy

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

!- Open the netcdf file:
   ios = nf90_open(path=chirps_filename,&
                   mode=NF90_NOWRITE,ncId=nid)
   call LIS_verify(ios,'Error in nf90_open in read_chirps2')
 
!- Read in Netcdf file .... 
   allocate( precip_nc(chirps2_struc(n)%nc, chirps2_struc(n)%nr) )
   precip_nc = LIS_rc%udef
   ios = nf90_inq_varid( nid, "precip", varid )
   call LIS_verify(ios,'Error in nf90_inq_varid in read_chirps2')

   ios = nf90_get_var( nid, varid, precip_nc, &
!              start=(/1,1,doy/),&
              start=(/chirps2_struc(n)%start_nc,chirps2_struc(n)%start_nr,doy/),&
              count=(/chirps2_struc(n)%nc,chirps2_struc(n)%nr,1/) )
   call LIS_verify(ios,'Error in nf90_get_var in read_chirps2')

   write(unit=cdoy, fmt='(i3.3)') doy
   write(LIS_logunit,*)"[INFO] Reading in CHIRPS-2 Day since 1980-1-1 (and DOY): "
   write(LIS_logunit,*) total_days, "("//cdoy//")"

!   latitude=LIS_rc%udef
!   ios = nf90_inq_varid( nid, "latitude", varid )
!   ios = nf90_get_var( nid, varid, latitude, &
!             start=(/1/),&
!              count=(/chirps2_struc(n)%nr/) )

! For day 121 (year: 2005)
!   write(*,*) precip_nc(4083,1115)
!   write(*,*) precip_nc(4129,862)

!- Close netCDF file.
   ios=nf90_close(nid)
   call LIS_verify(ios,'Error in nf90_close in read_chirps2')
!-
#endif

!- Transferring current data to 1-D array for interpolation
   write(LIS_logunit,*) "[INFO] Reprojecting the CHIRPS-2 domain to the LIS run domain ... "
   c=0
   allocate( chirpsprec_1d(chirps2_struc(n)%nc*chirps2_struc(n)%nr) )
   allocate( lb(chirps2_struc(n)%nc*chirps2_struc(n)%nr) )
   chirpsprec_1d = LIS_rc%udef
   do j=1,chirps2_struc(n)%nr
      do i=1,chirps2_struc(n)%nc
         c = c + 1
         chirpsprec_1d(c) = precip_nc(i,j)
      enddo
   enddo
   deallocate( precip_nc )

!- Interpolate if forcing and model grids are not both one deg.
   lb = .true.
   lo = .true.
   mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
   lisprec_1d = LIS_rc%udef

   select case( LIS_rc%met_interp(findex) )

    case( "budget-bilinear" )
      call conserv_interp( chirps2_struc(n)%gridDesc(:), &
                   lb, chirpsprec_1d, lo, lisprec_1d,&
                   chirps2_struc(n)%mi, mo, &
                   LIS_domain(n)%lat, LIS_domain(n)%lon,&
                   chirps2_struc(n)%w112, chirps2_struc(n)%w122,&
                   chirps2_struc(n)%w212, chirps2_struc(n)%w222,&
                   chirps2_struc(n)%n112, chirps2_struc(n)%n122,&
                   chirps2_struc(n)%n212, chirps2_struc(n)%n222,&
                   LIS_rc%udef, iret)

    case( "neighbor" )
      call neighbor_interp( chirps2_struc(n)%gridDesc(:), &
                    lb, chirpsprec_1d, lo, lisprec_1d,&
                    chirps2_struc(n)%mi, mo, &
                    LIS_domain(n)%lat, LIS_domain(n)%lon,&
                    chirps2_struc(n)%n113,LIS_rc%udef, iret)

     case( "average" )
        call upscaleByAveraging( chirps2_struc(n)%mi, &
               mo, LIS_rc%udef, &
               chirps2_struc(n)%n111, &
               lb, chirpsprec_1d, lo, lisprec_1d)

    end select

  ! Convert precip from the 1D to 2D LIS output gridspace:
    lisprec_2d = LIS_rc%udef
    c = 0
    do j = 1, LIS_rc%lnr(n)
       do i = 1, LIS_rc%lnc(n)
          if(LIS_domain(n)%gindex(i,j).ne.-1) then
             lisprec_2d(i,j) = lisprec_1d(i+c)
          endif
       enddo
       c = c + LIS_rc%lnc(n)
    enddo

  deallocate( lb, chirpsprec_1d )

   do j = 1,LIS_rc%lnr(n)
      do i = 1,LIS_rc%lnc(n)
         if( lisprec_2d(i,j) >= 0. ) then
            index1 = LIS_domain(n)%gindex(i,j)
            chirps2_struc(n)%metdata2(kk,1,index1) = lisprec_2d(i,j)  ! mm/day
         endif
      enddo
   enddo

end subroutine read_chirps2

