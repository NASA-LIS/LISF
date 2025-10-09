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
! SUBROUTINE: WSF_ARFS_RESAMPLE
!
! DESCRIPTION: Resample WSF TB to Air Force Grid
!
!-------------------------------------------------------------------------

subroutine WSF_ARFS_RESAMPLE(WSF_filename, output_dir)

    USE VARIABLES
    USE DATADOMAIN
    USE FUNCTIONS
    USE TOOLSUBS_WSF
    USE invdist_wsf2arfs
    USE LDT_logMod
    USE netcdf
    
    IMPLICIT NONE
    
    ! Arguments
    character(len=*), intent(in) :: WSF_filename
    character(len=*), intent(in) :: output_dir
    
    ! Local variables
    integer :: ierr, nscans, nfovs, nchans
    real*4, allocatable :: tb_lowres(:,:,:)
    real*4, allocatable :: lat_in(:,:), lon_in(:,:)
    real*4, allocatable :: land_frac_low(:,:)
    integer*1, allocatable :: quality_flag_in(:,:,:)
    real*4, allocatable :: earth_inc_angle(:,:,:)
    
    real*8, allocatable :: ARFS_LAT(:), ARFS_LON(:)
    real*4, allocatable :: ARFS_TB(:,:,:)
    real*4, allocatable :: ARFS_LAND_FRAC(:,:)
    integer*1, allocatable :: ARFS_QUALITY_FLAG(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_COUNT(:,:)
    
    character(len=255) :: output_filename
    character(len=100) :: basename
    integer :: filename_start_pos, ncid, varid
    integer :: lon_dimid, lat_dimid, chan_dimid
    integer :: lon_varid, lat_varid
    integer :: tb_varid, land_varid, qf_varid, count_varid
    integer :: ichan
    real, allocatable :: lats(:), lons(:)
    character(len=20) :: chan_name
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Starting WSF ARFS Resampling'
    write(LDT_logunit,*)'[INFO] Input file: ', trim(WSF_filename)
    write(LDT_logunit,*)'[INFO] Using LOW resolution (30 km) data'
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Read WSF data
    write(LDT_logunit,*)'[INFO] Reading WSF NetCDF file...'
    call get_wsf_data(WSF_filename, &
        tb_lowres, lat_in, lon_in, land_frac_low, &
        quality_flag_in, earth_inc_angle, &
        nscans, nfovs, nchans, ierr)
    
    if (ierr /= 0) then
        write(LDT_logunit,*)'[ERR] Failed to read WSF data'
        return
    end if
    
    ! Setup ARFS grid
    write(LDT_logunit,*)'[INFO] Setting up ARFS grid...'
    call ARFS_GEO
    allocate(ARFS_LAT(arfs_nrow_lat))
    allocate(ARFS_LON(arfs_mcol_lon))
    ARFS_LAT = LAT(arfs_geo_lat_lo, arfs_geo_lat_up, -arfs_lat_space)
    ARFS_LON = LON(arfs_geo_lon_lf, arfs_geo_lon_rt, arfs_lon_space)
    
    write(LDT_logunit,*)'[INFO] ARFS grid dimensions: ', &
        arfs_mcol_lon, ' x ', arfs_nrow_lat
    
    ! Allocate output arrays
    allocate(ARFS_TB(nchans, arfs_mcol_lon, arfs_nrow_lat))
    allocate(ARFS_LAND_FRAC(arfs_mcol_lon, arfs_nrow_lat))
    allocate(ARFS_QUALITY_FLAG(arfs_mcol_lon, arfs_nrow_lat))
    allocate(ARFS_SAMPLE_COUNT(arfs_mcol_lon, arfs_nrow_lat))
    
    ! Perform inverse distance resampling
    write(LDT_logunit,*)'[INFO] Performing inverse distance resampling...'
    call WSF2ARFS_INVDIS(tb_lowres, lat_in, lon_in, &
        quality_flag_in, land_frac_low, &
        nscans, nfovs, nchans, &
        ARFS_LAT, ARFS_LON, &
        ARFS_TB, ARFS_QUALITY_FLAG, ARFS_LAND_FRAC, &
        ARFS_SAMPLE_COUNT)
    
    ! Generate output filename
    filename_start_pos = index(WSF_filename, '/', back=.true.) + 1
    basename = WSF_filename(filename_start_pos:)
    output_filename = trim(output_dir) // '/WSF_resampled_LOW_' // trim(basename)
    
    write(LDT_logunit,*)'[INFO] Writing output to: ', trim(output_filename)
    
    ! Create NetCDF output file
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    ierr = nf90_create(trim(output_filename), NF90_NETCDF4, ncid)
    if (ierr /= NF90_NOERR) then
        write(LDT_logunit,*)'[ERR] Failed to create output file'
        write(LDT_logunit,*)nf90_strerror(ierr)
        return
    end if
    
    ! Define dimensions
    ierr = nf90_def_dim(ncid, 'lon', arfs_mcol_lon, lon_dimid)
    ierr = nf90_def_dim(ncid, 'lat', arfs_nrow_lat, lat_dimid)
    ierr = nf90_def_dim(ncid, 'channel', nchans, chan_dimid)
    
    ! Define coordinate variables
    ierr = nf90_def_var(ncid, 'lon', NF90_FLOAT, lon_dimid, lon_varid)
    ierr = nf90_put_att(ncid, lon_varid, 'units', 'degrees_east')
    ierr = nf90_put_att(ncid, lon_varid, 'standard_name', 'longitude')
    ierr = nf90_put_att(ncid, lon_varid, 'long_name', 'Longitude')
    
    ierr = nf90_def_var(ncid, 'lat', NF90_FLOAT, lat_dimid, lat_varid)
    ierr = nf90_put_att(ncid, lat_varid, 'units', 'degrees_north')
    ierr = nf90_put_att(ncid, lat_varid, 'standard_name', 'latitude')
    ierr = nf90_put_att(ncid, lat_varid, 'long_name', 'Latitude')
    
    ! Define data variables
    ierr = nf90_def_var(ncid, 'brightness_temperature', NF90_FLOAT, &
        [lon_dimid, lat_dimid, chan_dimid], tb_varid)
    ierr = nf90_put_att(ncid, tb_varid, 'units', 'K')
    ierr = nf90_put_att(ncid, tb_varid, 'long_name', &
        'Brightness Temperature - Resampled to ARFS Grid')
    ierr = nf90_put_att(ncid, tb_varid, '_FillValue', 0.0)
    
    ierr = nf90_def_var(ncid, 'land_fraction', NF90_FLOAT, &
        [lon_dimid, lat_dimid], land_varid)
    ierr = nf90_put_att(ncid, land_varid, 'units', 'fraction')
    ierr = nf90_put_att(ncid, land_varid, 'long_name', 'Land Fraction')
    ierr = nf90_put_att(ncid, land_varid, '_FillValue', 0.0)
    
    ierr = nf90_def_var(ncid, 'quality_flag', NF90_BYTE, &
        [lon_dimid, lat_dimid], qf_varid)
    ierr = nf90_put_att(ncid, qf_varid, 'long_name', 'Quality Flag')
    
    ierr = nf90_def_var(ncid, 'sample_count', NF90_INT, &
        [lon_dimid, lat_dimid], count_varid)
    ierr = nf90_put_att(ncid, count_varid, 'long_name', &
        'Number of input samples per output pixel')
    
    ! Global attributes
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.10')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'title', &
        'WSF Data Resampled to ARFS Grid')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'institution', &
        'NASA GSFC Hydrological Sciences Laboratory')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'source', &
        'Land Data Toolkit (LDT)')
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'source_file', &
        trim(WSF_filename))
    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'resolution', &
        'Low resolution 30 km')
    
    ! End define mode
    ierr = nf90_enddef(ncid)
    
    ! Write coordinate data
    allocate(lats(arfs_nrow_lat))
    allocate(lons(arfs_mcol_lon))
    lats = real(ARFS_LAT)
    lons = real(ARFS_LON)
    
    ierr = nf90_put_var(ncid, lat_varid, lats)
    ierr = nf90_put_var(ncid, lon_varid, lons)
    
    ! Write data variables
    ierr = nf90_put_var(ncid, tb_varid, ARFS_TB)
    ierr = nf90_put_var(ncid, land_varid, ARFS_LAND_FRAC)
    ierr = nf90_put_var(ncid, qf_varid, ARFS_QUALITY_FLAG)
    ierr = nf90_put_var(ncid, count_varid, ARFS_SAMPLE_COUNT)
    
    ! Close file
    ierr = nf90_close(ncid)
    
    deallocate(lats, lons)
    
    write(LDT_logunit,*)'[INFO] Successfully wrote output file'
#else
    write(LDT_logunit,*)'[ERR] LDT not compiled with NetCDF support'
#endif
    
    ! Cleanup
    deallocate(tb_lowres)
    deallocate(lat_in, lon_in)
    deallocate(land_frac_low)
    deallocate(quality_flag_in, earth_inc_angle)
    deallocate(ARFS_LAT, ARFS_LON)
    deallocate(ARFS_TB, ARFS_LAND_FRAC)
    deallocate(ARFS_QUALITY_FLAG, ARFS_SAMPLE_COUNT)
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] WSF ARFS Resampling Complete'
    write(LDT_logunit,*)'[INFO] ========================================='

end subroutine WSF_ARFS_RESAMPLE