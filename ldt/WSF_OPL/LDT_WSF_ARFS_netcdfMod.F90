!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LDT_WSF_ARFS_netcdfMod
!
! DESCRIPTION: NetCDF writer for WSF hourly stitched data on ARFS grid
!              FIXED: Proper time format with hour, quality flag documentation
!
!-------------------------------------------------------------------------
#include "LDT_misc.h"

module LDT_WSF_ARFS_netcdfMod

  implicit none
  private
  
  public :: LDT_WSF_ARFS_write_netcdf_hourly

contains

  subroutine LDT_WSF_ARFS_write_netcdf_hourly(nc, nr, &
      tb_10h, tb_10v, tb_18h, tb_18v, &
      tb_23h, tb_23v, tb_36h, tb_36v, &
      tb_89h, tb_89v, land_frac, &
      quality_flag, sample_v, sample_h, &
      lat, lon, output_fname, &
      yyyymmdd, hour_str, n_files)
      
    use LDT_coreMod
    use LDT_logMod
    use netcdf
    use ESMF
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: nc, nr  ! Grid dimensions
    real*4, intent(in) :: tb_10h(nc,nr), tb_10v(nc,nr)
    real*4, intent(in) :: tb_18h(nc,nr), tb_18v(nc,nr)
    real*4, intent(in) :: tb_23h(nc,nr), tb_23v(nc,nr)
    real*4, intent(in) :: tb_36h(nc,nr), tb_36v(nc,nr)
    real*4, intent(in) :: tb_89h(nc,nr), tb_89v(nc,nr)
    real*4, intent(in) :: land_frac(nc,nr)
    integer*1, intent(in) :: quality_flag(nc,nr)
    integer*4, intent(in) :: sample_v(nc,nr), sample_h(nc,nr)
    real*8, intent(in) :: lat(nr), lon(nc)
    character(len=*), intent(in) :: output_fname
    character(len=8), intent(in) :: yyyymmdd
    character(len=2), intent(in) :: hour_str
    integer, intent(in) :: n_files
    
    ! Fill value constant
    real*4, parameter :: FILL_VALUE = -9999.0
    
    ! Local variables
    integer :: ncid, dimids(2)
    integer :: lat_dimid, lon_dimid, time_dimid
    integer :: lat_varid, lon_varid, time_varid, hour_varid
    integer :: tb_10h_varid, tb_10v_varid
    integer :: tb_18h_varid, tb_18v_varid
    integer :: tb_23h_varid, tb_23v_varid
    integer :: tb_36h_varid, tb_36v_varid
    integer :: tb_89h_varid, tb_89v_varid
    integer :: land_frac_varid, qf_varid
    integer :: sample_v_varid, sample_h_varid
    integer :: iret
    real, allocatable :: lats(:), lons(:)
    real*8 :: time_seconds
    integer :: year, month, day, hour, minute
    type(ESMF_Time) :: file_time, reference_time
    type(ESMF_TimeInterval) :: time_diff
    integer :: rc_time
    character(len=50) :: time_units_str, time_long_name
    character(len=100) :: history_str
    character(len=200) :: qf_desc_str
    integer :: hour_value
    
    ! Only master process writes the file
    if (LDT_masterproc) then
       write(LDT_logunit,*)'[INFO] Creating NetCDF file: ', trim(output_fname)
       write(LDT_logunit,*)'[INFO] Using fill value: ', FILL_VALUE
       
       ! Create the output file
#if(defined USE_NETCDF4)
       iret = nf90_create(path=trim(output_fname), &
                         cmode=NF90_NETCDF4, ncid=ncid)
#else
       iret = nf90_create(path=trim(output_fname), &
                         cmode=NF90_CLOBBER, ncid=ncid)
#endif
       call LDT_verify(iret, '[ERR] nf90_create failed')
       
       ! Define dimensions
       call LDT_verify(nf90_def_dim(ncid, 'lon', nc, lon_dimid), &
            '[ERR] nf90_def_dim failed for lon')
       call LDT_verify(nf90_def_dim(ncid, 'lat', nr, lat_dimid), &
            '[ERR] nf90_def_dim failed for lat')
       call LDT_verify(nf90_def_dim(ncid, 'time', 1, time_dimid), &
            '[ERR] nf90_def_dim failed for time')
       
       ! =====================================================================
       ! COORDINATE VARIABLES
       ! =====================================================================
       
       ! Longitude
       call LDT_verify(nf90_def_var(ncid, 'lon', NF90_FLOAT, &
            lon_dimid, lon_varid), &
            '[ERR] nf90_def_var failed for lon')
       call LDT_verify(nf90_put_att(ncid, lon_varid, 'long_name', &
            'longitude'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, lon_varid, 'standard_name', &
            'longitude'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, lon_varid, 'units', &
            'degrees_east'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, lon_varid, 'axis', &
            'X'), '[ERR] nf90_put_att failed')
       
       ! Latitude
       call LDT_verify(nf90_def_var(ncid, 'lat', NF90_FLOAT, &
            lat_dimid, lat_varid), &
            '[ERR] nf90_def_var failed for lat')
       call LDT_verify(nf90_put_att(ncid, lat_varid, 'long_name', &
            'latitude'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, lat_varid, 'standard_name', &
            'latitude'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, lat_varid, 'units', &
            'degrees_north'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, lat_varid, 'axis', &
            'Y'), '[ERR] nf90_put_att failed')
       
       ! Time (seconds since epoch)
       call LDT_verify(nf90_def_var(ncid, 'time', NF90_DOUBLE, &
            time_dimid, time_varid), &
            '[ERR] nf90_def_var failed for time')
       
       ! Build time description with date and hour
       write(time_long_name, '(A,A8,A,A2,A)') 'time for ', yyyymmdd, ' hour ', hour_str, ':00'
       call LDT_verify(nf90_put_att(ncid, time_varid, 'long_name', &
            trim(time_long_name)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, time_varid, 'standard_name', &
            'time'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, time_varid, 'units', &
            'seconds since 1970-01-01 00:00:00'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, time_varid, 'calendar', &
            'gregorian'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, time_varid, 'axis', &
            'T'), '[ERR] nf90_put_att failed')
       
       ! Build readable time string
       write(time_units_str, '(A,A2,A)') yyyymmdd, hour_str, ':00:00'
       call LDT_verify(nf90_put_att(ncid, time_varid, 'time_coverage_start', &
            trim(time_units_str)), '[ERR] nf90_put_att failed')
       
       ! Hour variable (for easy reference)
       call LDT_verify(nf90_def_var(ncid, 'hour', NF90_INT, &
            time_dimid, hour_varid), &
            '[ERR] nf90_def_var failed for hour')
       call LDT_verify(nf90_put_att(ncid, hour_varid, 'long_name', &
            'hour of day'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, hour_varid, 'units', &
            'hours'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, hour_varid, 'valid_range', &
            (/0, 23/)), '[ERR] nf90_put_att failed')
       
       ! =====================================================================
       ! DATA VARIABLES
       ! =====================================================================
       
       dimids = (/lon_dimid, lat_dimid/)
       
       ! TB_10H
       call LDT_verify(nf90_def_var(ncid, 'TB_10H', NF90_FLOAT, &
            dimids, tb_10h_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_10h_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_10h_varid, 'long_name', &
            'Brightness Temperature at 10.65 GHz H-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_10h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_10h_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_10h_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_10V
       call LDT_verify(nf90_def_var(ncid, 'TB_10V', NF90_FLOAT, &
            dimids, tb_10v_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_10v_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_10v_varid, 'long_name', &
            'Brightness Temperature at 10.65 GHz V-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_10v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_10v_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_10v_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_18H
       call LDT_verify(nf90_def_var(ncid, 'TB_18H', NF90_FLOAT, &
            dimids, tb_18h_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_18h_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_18h_varid, 'long_name', &
            'Brightness Temperature at 18.7 GHz H-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_18h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_18h_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_18h_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_18V
       call LDT_verify(nf90_def_var(ncid, 'TB_18V', NF90_FLOAT, &
            dimids, tb_18v_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_18v_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_18v_varid, 'long_name', &
            'Brightness Temperature at 18.7 GHz V-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_18v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_18v_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_18v_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_23H
       call LDT_verify(nf90_def_var(ncid, 'TB_23H', NF90_FLOAT, &
            dimids, tb_23h_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_23h_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_23h_varid, 'long_name', &
            'Brightness Temperature at 23.8 GHz H-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_23h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_23h_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_23h_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_23V
       call LDT_verify(nf90_def_var(ncid, 'TB_23V', NF90_FLOAT, &
            dimids, tb_23v_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_23v_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_23v_varid, 'long_name', &
            'Brightness Temperature at 23.8 GHz V-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_23v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_23v_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_23v_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_36H
       call LDT_verify(nf90_def_var(ncid, 'TB_36H', NF90_FLOAT, &
            dimids, tb_36h_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_36h_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_36h_varid, 'long_name', &
            'Brightness Temperature at 36.5 GHz H-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_36h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_36h_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_36h_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_36V
       call LDT_verify(nf90_def_var(ncid, 'TB_36V', NF90_FLOAT, &
            dimids, tb_36v_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_36v_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_36v_varid, 'long_name', &
            'Brightness Temperature at 36.5 GHz V-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_36v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_36v_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_36v_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_89H
       call LDT_verify(nf90_def_var(ncid, 'TB_89H', NF90_FLOAT, &
            dimids, tb_89h_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_89h_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_89h_varid, 'long_name', &
            'Brightness Temperature at 89.0 GHz H-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_89h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_89h_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_89h_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! TB_89V
       call LDT_verify(nf90_def_var(ncid, 'TB_89V', NF90_FLOAT, &
            dimids, tb_89v_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, tb_89v_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, tb_89v_varid, 'long_name', &
            'Brightness Temperature at 89.0 GHz V-pol'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_89v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_89v_varid, 'valid_range', &
            (/100.0, 350.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, tb_89v_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! Land fraction
       call LDT_verify(nf90_def_var(ncid, 'LAND_FRAC', NF90_FLOAT, &
            dimids, land_frac_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_def_var_fill(ncid, land_frac_varid, 0, FILL_VALUE), &
            '[ERR] nf90_def_var_fill failed')
       call LDT_verify(nf90_put_att(ncid, land_frac_varid, 'long_name', &
            'Land Fraction'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, land_frac_varid, 'units', &
            '1'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, land_frac_varid, 'valid_range', &
            (/0.0, 1.0/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, land_frac_varid, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       ! =====================================================================
       ! QUALITY FLAG with proper documentation
       ! =====================================================================
       
       call LDT_verify(nf90_def_var(ncid, 'QUALITY_FLAG', NF90_BYTE, &
            dimids, qf_varid), '[ERR] nf90_def_var failed')
       
      call LDT_verify(nf90_put_att(ncid, qf_varid, 'long_name', &
            'Quality Flag (8-bit packed)'), '[ERR] nf90_put_att failed')
       
       ! Detailed description for 8-bit quality flag
       qf_desc_str = 'Bit-packed quality indicators: ' // &
                     'Bit 0 = Ocean (1=ocean/water, 0=land); ' // &
                     'Bit 1 = Precipitation (1=precip detected); ' // &
                     'Bit 2 = Snow (1=snow detected); ' // &
                     'Bit 3 = Sensor quality 10GHz (1=bad); ' // &
                     'Bit 4 = Sensor quality 18GHz (1=bad); ' // &
                     'Bit 5 = Sensor quality 23GHz (1=bad); ' // &
                     'Bit 6 = Sensor quality 36GHz (1=bad); ' // &
                     'Bit 7 = Sensor quality 89GHz (1=bad)'
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'description', &
            trim(qf_desc_str)), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'flag_masks', &
            (/1, 2, 4, 8, 16, 32, 64, 128/)), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'flag_meanings', &
            'ocean precipitation snow sensor_10ghz sensor_18ghz ' // &
            'sensor_23ghz sensor_36ghz sensor_89ghz'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'valid_range', &
            (/0, 255/)), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_put_att(ncid, qf_varid, '_FillValue', &
            int(-1,kind=1)), '[ERR] nf90_put_att failed')
       
       ! =====================================================================
       ! SAMPLE COUNT VARIABLES
       ! =====================================================================
       
       call LDT_verify(nf90_def_var(ncid, 'SAMPLE_V', NF90_INT, &
            dimids, sample_v_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_put_att(ncid, sample_v_varid, 'long_name', &
            'Number of V-pol samples averaged'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, sample_v_varid, 'description', &
            'Number of observations averaged for V-pol channels in hourly stitching'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, sample_v_varid, 'units', &
            'count'), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_def_var(ncid, 'SAMPLE_H', NF90_INT, &
            dimids, sample_h_varid), '[ERR] nf90_def_var failed')
       call LDT_verify(nf90_put_att(ncid, sample_h_varid, 'long_name', &
            'Number of H-pol samples averaged'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, sample_h_varid, 'description', &
            'Number of observations averaged for H-pol channels in hourly stitching'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, sample_h_varid, 'units', &
            'count'), '[ERR] nf90_put_att failed')
       
       ! =====================================================================
       ! GLOBAL ATTRIBUTES
       ! =====================================================================
       
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'title', &
            'WSF Brightness Temperature Data Resampled to ARFS Grid (Hourly Stitched)'), &
            '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'institution', &
            'NASA GSFC'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'source', &
            'WSF L1R SDR Data'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'sensor', &
            'WSF Microwave Radiometer'), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'grid_resolution', &
            '0.09375 x 0.140625 degrees'), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'grid_type', &
            'ARFS (Air Force) Grid'), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'processing_date', &
            yyyymmdd), '[ERR] nf90_put_att failed')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'processing_hour', &
            hour_str), '[ERR] nf90_put_att failed')
       
       ! Build time coverage string
       write(time_units_str, '(A,A2,A)') yyyymmdd, hour_str, ':00:00'
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'time_coverage', &
            trim(time_units_str)), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'number_of_input_files', &
            n_files), '[ERR] nf90_put_att failed')
       
       write(history_str, '(A,I3,A)') 'Stitched from ', n_files, &
            ' WSF files with mean averaging for overlapping regions and snow/precip filtering'
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'history', &
            trim(history_str)), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, '_FillValue', &
            FILL_VALUE), '[ERR] nf90_put_att failed')
       
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', &
            'CF-1.8'), '[ERR] nf90_put_att failed')
       
       ! End definition mode
       call LDT_verify(nf90_enddef(ncid), '[ERR] nf90_enddef failed')
       
       ! =====================================================================
       ! WRITE DATA
       ! =====================================================================
       
       ! Write coordinate data
       allocate(lats(nr))
       allocate(lons(nc))
       lats = real(lat)
       lons = real(lon)
       
       call LDT_verify(nf90_put_var(ncid, lat_varid, lats), &
            '[ERR] nf90_put_var failed for lat')
       call LDT_verify(nf90_put_var(ncid, lon_varid, lons), &
            '[ERR] nf90_put_var failed for lon')
       
       ! Calculate time value (seconds since epoch)
       read(yyyymmdd(1:4), '(I4)') year
       read(yyyymmdd(5:6), '(I2)') month
       read(yyyymmdd(7:8), '(I2)') day
       read(hour_str, '(I2)') hour
       minute = 0
       
       call ESMF_TimeSet(reference_time, yy=1970, mm=1, dd=1, &
                        h=0, m=0, s=0, rc=rc_time)
       call ESMF_TimeSet(file_time, yy=year, mm=month, dd=day, &
                        h=hour, m=minute, s=0, rc=rc_time)
       
       time_diff = file_time - reference_time
       call ESMF_TimeIntervalGet(time_diff, s_r8=time_seconds, rc=rc_time)
       
       call LDT_verify(nf90_put_var(ncid, time_varid, time_seconds), &
            '[ERR] nf90_put_var failed for time')
       
       ! Write hour value
       hour_value = hour
       call LDT_verify(nf90_put_var(ncid, hour_varid, hour_value), &
            '[ERR] nf90_put_var failed for hour')
       
       write(LDT_logunit,*)'[INFO] Writing data for: ', trim(time_units_str)
       write(LDT_logunit,*)'[INFO] Time (seconds since epoch): ', time_seconds
       write(LDT_logunit,*)'[INFO] Hour: ', hour_value
       
       ! Write TB data
       call LDT_verify(nf90_put_var(ncid, tb_10h_varid, tb_10h), &
            '[ERR] nf90_put_var failed for TB_10H')
       call LDT_verify(nf90_put_var(ncid, tb_10v_varid, tb_10v), &
            '[ERR] nf90_put_var failed for TB_10V')
       call LDT_verify(nf90_put_var(ncid, tb_18h_varid, tb_18h), &
            '[ERR] nf90_put_var failed for TB_18H')
       call LDT_verify(nf90_put_var(ncid, tb_18v_varid, tb_18v), &
            '[ERR] nf90_put_var failed for TB_18V')
       call LDT_verify(nf90_put_var(ncid, tb_23h_varid, tb_23h), &
            '[ERR] nf90_put_var failed for TB_23H')
       call LDT_verify(nf90_put_var(ncid, tb_23v_varid, tb_23v), &
            '[ERR] nf90_put_var failed for TB_23V')
       call LDT_verify(nf90_put_var(ncid, tb_36h_varid, tb_36h), &
            '[ERR] nf90_put_var failed for TB_36H')
       call LDT_verify(nf90_put_var(ncid, tb_36v_varid, tb_36v), &
            '[ERR] nf90_put_var failed for TB_36V')
       call LDT_verify(nf90_put_var(ncid, tb_89h_varid, tb_89h), &
            '[ERR] nf90_put_var failed for TB_89H')
       call LDT_verify(nf90_put_var(ncid, tb_89v_varid, tb_89v), &
            '[ERR] nf90_put_var failed for TB_89V')
       
       ! Write auxiliary data
       call LDT_verify(nf90_put_var(ncid, land_frac_varid, land_frac), &
            '[ERR] nf90_put_var failed for LAND_FRAC')
       call LDT_verify(nf90_put_var(ncid, qf_varid, quality_flag), &
            '[ERR] nf90_put_var failed for QUALITY_FLAG')
       call LDT_verify(nf90_put_var(ncid, sample_v_varid, sample_v), &
            '[ERR] nf90_put_var failed for SAMPLE_V')
       call LDT_verify(nf90_put_var(ncid, sample_h_varid, sample_h), &
            '[ERR] nf90_put_var failed for SAMPLE_H')
       
       ! Close the file
       call LDT_verify(nf90_close(ncid), '[ERR] nf90_close failed')
       
       deallocate(lats)
       deallocate(lons)
       
       write(LDT_logunit,*)'[INFO] Successfully wrote hourly stitched file: ', &
                          trim(output_fname)
       write(LDT_logunit,*)'[INFO] Time coverage: ', trim(time_units_str)
    endif
    
  end subroutine LDT_WSF_ARFS_write_netcdf_hourly

end module LDT_WSF_ARFS_netcdfMod