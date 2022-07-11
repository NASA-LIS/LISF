!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"

! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif

module USAF_precipBiasCorrMod

  ! DESCRIPTION:
  ! This module implements routines to calculate and apply bias ratios for
  ! select precipitation inputs to NAFPA, based on climatology.
  !
  ! REVISION HISTORY:
  ! 25 May 2022: Eric Kemp/SSAI: First version

  ! Defaults
  implicit none
  private

  ! Public type for calculating bias ratios for precipitation
  type USAF_PrecipBiasRatio_t
     private
     integer :: ncols
     integer :: nrows
     real, allocatable :: bias_ratio(:,:)
   contains
     ! Object methods
     procedure :: new => USAF_precipBiasRatio_new
     procedure :: delete => USAF_precipBiasRatio_delete
     procedure :: update_back_ratio => USAF_precipBiasRatio_update_back_ratio
  end type USAF_PrecipBiasRatio_t

contains

  ! Constructor method
  subroutine USAF_precipBiasRatio_new(this, ncols, nrows)
    implicit none
    class(USAF_PrecipBiasRatio_t), intent(inout) :: this
    integer, intent(in) :: ncols
    integer, intent(in) :: nrows
    this%ncols = ncols
    this%nrows = nrows
    allocate(this%bias_ratio(ncols, nrows))
    this%bias_ratio = 1
  end subroutine USAF_precipBiasRatio_new

  ! Destructor method
  subroutine USAF_precipBiasRatio_delete(this)
    implicit none
    class(USAF_precipBiasRatio_t), intent(inout) :: this
    this%ncols = 0
    this%nrows = 0
    deallocate(this%bias_ratio)
  end subroutine USAF_precipBiasRatio_delete

  ! Handle calculation of bias_ratio
  subroutine USAF_precipBiasRatio_update_back_ratio(this, first_guess_source, &
       imonth, n)

    ! Imports
    use AGRMET_forcingMod, only: agrmet_struc
    use LIS_logMod, only: LIS_logunit, LIS_endrun

    ! Defaults
    implicit none

    ! Arguments
    class(USAF_precipBiasRatio_t), intent(inout) :: this
    character(*), intent(in) :: first_guess_source
    integer, intent(in) :: imonth
    integer, intent(in) :: n

    ! Locals
    integer :: ncols, nrows
    real, allocatable :: chelsa_climo(:,:)
    real, allocatable :: back_climo(:,:)
    real, allocatable :: chelsa_udef, back_udef
    integer :: c, r

    ! Get new month of chelsa climo data
    nrows = this%nrows
    ncols = this%ncols
    allocate(chelsa_climo(ncols, nrows))
    chelsa_climo = 0
    call fetch_climo(agrmet_struc(n)%chelsa_climo_file, imonth, ncols, nrows, &
         chelsa_climo, chelsa_udef)

    ! Get new month of background climo data
    allocate(back_climo(ncols, nrows))
    back_climo = 0
    if (trim(first_guess_source) .eq. "GFS") then
       call fetch_climo(agrmet_struc(n)%gfs_climo_file, imonth, ncols, nrows, &
            back_climo, back_udef)
    else if (trim(first_guess_source) .eq. "GALWEM") then
       call fetch_climo(agrmet_struc(n)%galwem_climo_file, imonth, &
            ncols, nrows, back_climo, back_udef)
    else
       write(LIS_logunit,*)'[ERR] Invalid selection for NAFPA background field'
       write(LIS_logunit,*)'[ERR] Expected GFS or GALWEM'
       write(LIS_logunit,*)'[ERR] Found ', trim(first_guess_source)
       call LIS_endrun()
    end if

    ! Calculate new bias ratio
    do r = 1, nrows
       do c = 1, ncols
          this%bias_ratio(c,r) = 1.0
          if (chelsa_climo(c,r) .ne. chelsa_udef) then
             if (back_climo(c,r) .ne. back_udef) then
                this%bias_ratio(c,r) = (chelsa_climo(c,r) + 1.e-3) / &
                     (back_climo(c,r) + 1.e-3)
             end if
          end if
       end do
    end do

    ! Clean up
    deallocate(chelsa_climo)
    deallocate(back_climo)
  end subroutine USAF_precipBiasRatio_update_back_ratio

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  ! Fetch desired month of remapped precip climatology
  subroutine fetch_climo(filename, imonth, ncols, nrows, climo_field, udef)

    ! Imports
    use LIS_logMod, only: LIS_logunit, LIS_endrun
    use netcdf

    ! Defaults
    implicit none

    ! Arguments
    character(*), intent(in) :: filename
    integer, intent(in) :: imonth
    integer, intent(in) :: ncols
    integer, intent(in) :: nrows
    real, intent(out) :: climo_field(ncols, nrows)
    real, intent(out) :: udef

    ! Locals
    logical :: found_inq
    integer :: ncid, iret, varid, ndims
    integer :: dimids(3), dims(3), start(3), count(3), stride(3)
    integer :: i

    ! Ensure file exists.
    inquire(file=trim(filename), exist=found_inq)
    if (.not. found_inq) then
       write(LIS_logunit,*)'[ERR] Cannot find ', trim(filename)
       write(LIS_logunit,*)'[ERR] Stopping...'
       call LIS_endrun()
    end if
    write(LIS_logunit,*)'[INFO] Reading ', trim(filename)

    ! Get netCDF file ID.
    iret = nf90_open(path=trim(filename), mode=NF90_NOWRITE, ncid=ncid)
    if (iret .ne. NF90_NOERR) then
       write(LIS_logunit,*)'[ERR] nf90_open failed for ', trim(filename)
       call LIS_endrun()
    end if

    ! Get netCDF variable ID
    iret = nf90_inq_varid(ncid, 'PPTCLIM', varid)
    if (iret .ne. NF90_NOERR) then
       write(LIS_logunit,*)'[ERR] nf90_inq_varid failed for PPTCLIM'
       call LIS_endrun()
    end if

    ! Get number of variable dimensions
    iret = nf90_inquire_variable(ncid, varid, ndims=ndims)
    if (iret .ne. NF90_NOERR) then
       write(LIS_logunit,*)'[ERR] nf90_inquire_variable failed for PPTCLIM'
       call LIS_endrun()
    end if
    if (ndims .ne. 3) then
       write(LIS_logunit,*)'[ERR] Unexpected shape for PPTCLIM variable'
       write(LIS_logunit,*)'[ERR] Expected 3, found ', ndims
       call LIS_endrun()
    end if

    ! Get netCDF dimension IDs
    iret = nf90_inquire_variable(ncid, varid, dimids=dimids)
    if (iret .ne. NF90_NOERR) then
       write(LIS_logunit,*)'[ERR] nf90_inquire_variable failed for PPTCLIM'
       call LIS_endrun()
    end if

    ! Get dimensions
    do i = 1, 3
       iret = nf90_inquire_dimension(ncid, dimids(i), len=dims(i))
       if (iret .ne. NF90_NOERR) then
          write(LIS_logunit,*)'[ERR] nf90_inquire_dimension failed for PPTCLIM'
          call LIS_endrun()
       end if
    end do

    ! Check dimensions
    if (dims(1) .ne. ncols .or. &
         dims(2) .ne. nrows .or. &
         dims(3) .ne. 12) then
       write(LIS_logunit,*)'[ERR] Shape mismatch for PPTCLIM'
       write(LIS_logunit,*)'[ERR] Expected ', ncols, nrows, 12
       write(LIS_logunit,*)'[ERR] Found ', dims(1), dims(2), dims(3)
       call LIS_endrun()
    end if

    ! Fetch one month slice of PPTCLIM
    start(1) = 1
    start(2) = 1
    start(3) = imonth
    count(1) = ncols
    count(2) = nrows
    count(3) = 1
    iret = nf90_get_var(ncid, varid, climo_field, start=start, count=count)
    if (iret .ne. NF90_NOERR) then
       write(LIS_logunit,*)'[ERR] nf90_get_var failed for PPTCLIM'
       call LIS_endrun()
    end if

    ! Get missing flag
    iret = nf90_get_att(ncid, varid, 'missing_value', udef)
    if (iret .ne. NF90_NOERR) then
       write(LIS_logunit,*)'[ERR] nf90_get_att failed for PPTCLIM'
       call LIS_endrun()
    end if

    ! Close the file
    iret = nf90_close(ncid)
    if (iret .ne. NF90_NOERR) then
       write(LIS_logunit,*)'[ERR] nf90_close failed!'
       call LIS_endrun()
    end if

  end subroutine fetch_climo

#else

  ! Dummy version if netCDF is missing
  subroutine fetch_climo(filename, imonth, climo_field, udef)
    use LIS_logMod, only: LIS_logunit, LIS_endrun
    implicit none
    character(*), intent(in) :: filename
    integer, intent(in) :: imonth
    real, intent(inout) :: climo_field(:,:)
    real, intent(out) :: udef
    write(LIS_logunit,*)'[ERR] LIS was compiled without netCDF support!'
    write(LIS_logunit,*)'[ERR] Cannot process precipitation climatology!'
    write(LIS_logunit,*)'[ERR] Recompile LIS with netCDF support and try again'
    call LDT_endrun()
  end subroutine fetch_climo
#endif

end module USAF_precipBiasCorrMod


