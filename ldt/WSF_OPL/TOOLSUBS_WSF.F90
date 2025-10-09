!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: TOOLSUBS_WSF
!
! DESCRIPTION: Module for reading WSF NetCDF data
!
!-------------------------------------------------------------------------

#include "LDT_misc.h"

MODULE TOOLSUBS_WSF
    USE LDT_logMod, only: LDT_logunit, LDT_endrun
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    USE netcdf
#endif
    IMPLICIT NONE

CONTAINS

    SUBROUTINE get_wsf_data(filename, &
        tb_lowres, lat, lon, land_frac_low, &
        quality_flag, earth_inc_angle, &
        nscans, nfovs, nchans, ierr)
    
    ! Arguments
    character(*), intent(in) :: filename
    real*4, allocatable, intent(out) :: tb_lowres(:,:,:)   ! (nchans, nscans, nfovs)
    real*4, allocatable, intent(out) :: lat(:,:)           ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: lon(:,:)           ! (nscans, nfovs)
    real*4, allocatable, intent(out) :: land_frac_low(:,:) ! (nscans, nfovs)
    integer*1, allocatable, intent(out) :: quality_flag(:,:,:) ! (nbands, nscans, nfovs)
    real*4, allocatable, intent(out) :: earth_inc_angle(:,:,:) ! (nbands, nscans, nfovs)
    integer, intent(out) :: nscans, nfovs, nchans
    integer, intent(out) :: ierr
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    ! Local variables
    integer :: ncid, varid
    integer :: nscanr_dimid, nfovr_dimid, nchan_dimid, nband_dimid
    integer :: nbands
    logical :: file_exists
    
    ierr = 0
    nscans = 0
    nfovs = 0
    nchans = 0
    
    ! Check if file exists
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
        write(LDT_logunit,*)'[ERR] Cannot find file ', trim(filename)
        ierr = 1
        return
    end if
    
    ! Open the NetCDF file
    ierr = nf90_open(trim(filename), NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot open file ', trim(filename)
        write(LDT_logunit,*)nf90_strerror(ierr)
        ierr = 1
        return
    end if
    
    write(LDT_logunit,*)'[INFO] Successfully opened ', trim(filename)
    
    ! Get dimensions
    ierr = nf90_inq_dimid(ncid, 'nScanR', nscanr_dimid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find dimension nScanR'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inquire_dimension(ncid, nscanr_dimid, len=nscans)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot get nScanR dimension'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inq_dimid(ncid, 'nFOVR', nfovr_dimid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find dimension nFOVR'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inquire_dimension(ncid, nfovr_dimid, len=nfovs)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot get nFOVR dimension'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inq_dimid(ncid, 'nChan', nchan_dimid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find dimension nChan'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inquire_dimension(ncid, nchan_dimid, len=nchans)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot get nChan dimension'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inq_dimid(ncid, 'nBand', nband_dimid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find dimension nBand'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    ierr = nf90_inquire_dimension(ncid, nband_dimid, len=nbands)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot get nBand dimension'
        call nf90_close(ncid)
        ierr = 1
        return
    end if
    
    write(LDT_logunit,*)'[INFO] Dimensions: nScanR=', nscans, &
        ' nFOVR=', nfovs, ' nChan=', nchans, ' nBand=', nbands
    
    ! Allocate arrays
    allocate(tb_lowres(nchans, nscans, nfovs))
    allocate(lat(nscans, nfovs))
    allocate(lon(nscans, nfovs))
    allocate(land_frac_low(nscans, nfovs))
    allocate(quality_flag(nbands, nscans, nfovs))
    allocate(earth_inc_angle(nbands, nscans, nfovs))
    
    ! Read Latitude
    ierr = nf90_inq_varid(ncid, 'Latitude', varid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find variable Latitude'
        call cleanup_and_exit()
        return
    end if
    ierr = nf90_get_var(ncid, varid, lat)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot read Latitude'
        call cleanup_and_exit()
        return
    end if
    write(LDT_logunit,*)'[INFO] Read Latitude: range [', &
        minval(lat), ',', maxval(lat), ']'
    
    ! Read Longitude
    ierr = nf90_inq_varid(ncid, 'Longitude', varid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find variable Longitude'
        call cleanup_and_exit()
        return
    end if
    ierr = nf90_get_var(ncid, varid, lon)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot read Longitude'
        call cleanup_and_exit()
        return
    end if
    write(LDT_logunit,*)'[INFO] Read Longitude: range [', &
        minval(lon), ',', maxval(lon), ']'
    
    ! Read TbLowRes
    ierr = nf90_inq_varid(ncid, 'TbLowRes', varid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find variable TbLowRes'
        call cleanup_and_exit()
        return
    end if
    ierr = nf90_get_var(ncid, varid, tb_lowres)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot read TbLowRes'
        call cleanup_and_exit()
        return
    end if
    write(LDT_logunit,*)'[INFO] Read TbLowRes'
    
    ! Read LandFractionLowRes
    ierr = nf90_inq_varid(ncid, 'LandFractionLowRes', varid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot find variable LandFractionLowRes'
        call cleanup_and_exit()
        return
    end if
    ierr = nf90_get_var(ncid, varid, land_frac_low)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Cannot read LandFractionLowRes'
        call cleanup_and_exit()
        return
    end if
    write(LDT_logunit,*)'[INFO] Read LandFractionLowRes'
    
    ! Read QualityFlag
    ierr = nf90_inq_varid(ncid, 'QualityFlag', varid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[WARN] Cannot find variable QualityFlag'
        quality_flag = 0
    else
        ierr = nf90_get_var(ncid, varid, quality_flag)
        if (ierr /= NF90_NOERR) then
            write(LDT_logunit,*)'[WARN] Cannot read QualityFlag'
            quality_flag = 0
        else
            write(LDT_logunit,*)'[INFO] Read QualityFlag'
        end if
    end if
    
    ! Read EarthIncidenceAngle
    ierr = nf90_inq_varid(ncid, 'EarthIncidenceAngle', varid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[WARN] Cannot find variable EarthIncidenceAngle'
        earth_inc_angle = 52.0  ! Default value
    else
        ierr = nf90_get_var(ncid, varid, earth_inc_angle)
        if (ierr /= NF90_NOERR) then
            write(LDT_logunit,*)'[WARN] Cannot read EarthIncidenceAngle'
            earth_inc_angle = 52.0
        else
            write(LDT_logunit,*)'[INFO] Read EarthIncidenceAngle'
        end if
    end if
    
    ! Close the file
    ierr = nf90_close(ncid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[WARN] Error closing file'
    end if
    
    write(LDT_logunit,*)'[INFO] Successfully read all WSF data'
    ierr = 0
    return
    
CONTAINS
    subroutine cleanup_and_exit()
        if (allocated(tb_lowres)) deallocate(tb_lowres)
        if (allocated(lat)) deallocate(lat)
        if (allocated(lon)) deallocate(lon)
        if (allocated(land_frac_low)) deallocate(land_frac_low)
        if (allocated(quality_flag)) deallocate(quality_flag)
        if (allocated(earth_inc_angle)) deallocate(earth_inc_angle)
        call nf90_close(ncid)
        ierr = 1
    end subroutine cleanup_and_exit
    
#else
    ! Dummy version if LDT was compiled w/o NetCDF support
    write(LDT_logunit,*) '[ERR] get_wsf_data called without NetCDF support!'
    write(LDT_logunit,*) '[ERR] Recompile LDT with NetCDF support and try again!'
    call LDT_endrun()
    ierr = 1
#endif
    
    END SUBROUTINE get_wsf_data

END MODULE TOOLSUBS_WSF