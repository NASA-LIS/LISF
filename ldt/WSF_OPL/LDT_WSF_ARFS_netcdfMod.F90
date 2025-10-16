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
! MODULE: LDT_WSF_ARFS_netcdfMod
!
! REVISION HISTORY:
!  16 Oct 2025: Ehsan Jalilvand. Initial implementation for WSF.
!
! DESCRIPTION: Outputs WSF resampled brightness temperatures in netCDF format.
!
!------------------------------------------------------------------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_WSF_ARFS_netcdfMod

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: LDT_WSF_ARFS_write_netcdf

contains

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ! Subroutine for writing WSF ARFS resampled data to netCDF
  subroutine LDT_WSF_ARFS_write_netcdf(nc, nr, &
       tb_10h, tb_10v, tb_18h, tb_18v, tb_23h, tb_23v, &
       tb_36h, tb_36v, tb_89h, tb_89v, land_frac, quality_flag, &
       sample_v, sample_h, lats, lons, output_fname)

    ! Imports
    use LDT_coreMod, only: LDT_rc, LDT_masterproc
    use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
#if ( defined SPMD)
    use mpi
#endif
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    ! Arguments
    integer, intent(in) :: nc, nr
    real*4, intent(in) :: tb_10h(nc,nr), tb_10v(nc,nr)
    real*4, intent(in) :: tb_18h(nc,nr), tb_18v(nc,nr)
    real*4, intent(in) :: tb_23h(nc,nr), tb_23v(nc,nr)
    real*4, intent(in) :: tb_36h(nc,nr), tb_36v(nc,nr)
    real*4, intent(in) :: tb_89h(nc,nr), tb_89v(nc,nr)
    real*4, intent(in) :: land_frac(nc,nr)
    integer*1, intent(in) :: quality_flag(nc,nr)
    integer*4, intent(in) :: sample_v(nc,nr), sample_h(nc,nr)
    real*8, intent(in) :: lats(nr), lons(nc)
    character(*), intent(in) :: output_fname

    ! Locals
    integer :: shuffle, deflate, deflate_level
    integer :: iret
    integer :: ncid
    integer :: lon_dimid, lat_dimid, time_dimid
    integer :: lon_varid, lat_varid, time_varid
    integer :: tb_10h_varid, tb_10v_varid, tb_18h_varid, tb_18v_varid
    integer :: tb_23h_varid, tb_23v_varid, tb_36h_varid, tb_36v_varid
    integer :: tb_89h_varid, tb_89v_varid
    integer :: lwf_varid, qf_varid, sv_varid, sh_varid
    real*4, allocatable :: lats_float(:), lons_float(:)
    integer :: ierr
    integer :: i, j

    ! Only the master process handles the file output
    if (LDT_masterproc) then
       write(LDT_logunit,*)'[INFO] ========================================'
       write(LDT_logunit,*)'[INFO] Creating NETCDF file:'
       write(LDT_logunit,*)'[INFO] ', trim(output_fname)
       write(LDT_logunit,*)'[INFO] ========================================'

       ! Copy netCDF compression settings
       shuffle = NETCDF_shuffle
       deflate = NETCDF_deflate
       deflate_level = NETCDF_deflate_level

       ! Create the output file
#if(defined USE_NETCDF3)
         iret=nf90_create(path=trim(output_fname),&
              cmode=NF90_CLOBBER, ncid=ncid)
         call LDT_verify(iret,&
              '[ERR] nf90_create failed')
#endif
#if(defined USE_NETCDF4)
         iret=nf90_create(path=trim(output_fname),&
              cmode=NF90_NETCDF4, ncid=ncid)
         call LDT_verify(iret, &
              '[ERR] nf90_create failed')
#endif

         ! Define dimensions
         call LDT_verify(nf90_def_dim(ncid, 'lon', nc, lon_dimid), &
              '[ERR] nf90_def_dim failed for lon')
         call LDT_verify(nf90_def_dim(ncid, 'lat', nr, lat_dimid), &
              '[ERR] nf90_def_dim failed for lat')
         call LDT_verify(nf90_def_dim(ncid, 'time', 1, time_dimid), &
              '[ERR] nf90_def_dim failed for time')

         ! Define coordinate variables
         call LDT_verify(nf90_def_var(ncid, 'lon', NF90_FLOAT, &
              lon_dimid, lon_varid), &
              '[ERR] nf90_def_var failed for lon')
         call LDT_verify(nf90_put_att(ncid, lon_varid, 'long_name', &
              'longitude'), &
              '[ERR] nf90_put_att failed for lon long_name')
         call LDT_verify(nf90_put_att(ncid, lon_varid, 'units', &
              'degrees_east'), &
              '[ERR] nf90_put_att failed for lon units')

         call LDT_verify(nf90_def_var(ncid, 'lat', NF90_FLOAT, &
              lat_dimid, lat_varid), &
              '[ERR] nf90_def_var failed for lat')
         call LDT_verify(nf90_put_att(ncid, lat_varid, 'long_name', &
              'latitude'), &
              '[ERR] nf90_put_att failed for lat long_name')
         call LDT_verify(nf90_put_att(ncid, lat_varid, 'units', &
              'degrees_north'), &
              '[ERR] nf90_put_att failed for lat units')

         call LDT_verify(nf90_def_var(ncid, 'time', NF90_DOUBLE, &
              time_dimid, time_varid), &
              '[ERR] nf90_def_var failed for time')
         call LDT_verify(nf90_put_att(ncid, time_varid, 'long_name', &
              'time'), &
              '[ERR] nf90_put_att failed for time long_name')
         call LDT_verify(nf90_put_att(ncid, time_varid, 'units', &
              'seconds since 1970-01-01 00:00:00'), &
              '[ERR] nf90_put_att failed for time units')

         ! Define TB data variables with fill values
         ! TB_10H
         call LDT_verify(nf90_def_var(ncid, 'TB_10H', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_10h_varid), &
              '[ERR] nf90_def_var failed for TB_10H')
         call LDT_verify(nf90_def_var_fill(ncid, tb_10h_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_10H')
         call LDT_verify(nf90_put_att(ncid, tb_10h_varid, 'long_name', &
              'brightness temperature 10.85 GHz H-pol'), &
              '[ERR] nf90_put_att failed for TB_10H')
         call LDT_verify(nf90_put_att(ncid, tb_10h_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_10H units')

         ! TB_10V
         call LDT_verify(nf90_def_var(ncid, 'TB_10V', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_10v_varid), &
              '[ERR] nf90_def_var failed for TB_10V')
         call LDT_verify(nf90_def_var_fill(ncid, tb_10v_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_10V')
         call LDT_verify(nf90_put_att(ncid, tb_10v_varid, 'long_name', &
              'brightness temperature 10.85 GHz V-pol'), &
              '[ERR] nf90_put_att failed for TB_10V')
         call LDT_verify(nf90_put_att(ncid, tb_10v_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_10V units')

         ! TB_18H
         call LDT_verify(nf90_def_var(ncid, 'TB_18H', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_18h_varid), &
              '[ERR] nf90_def_var failed for TB_18H')
         call LDT_verify(nf90_def_var_fill(ncid, tb_18h_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_18H')
         call LDT_verify(nf90_put_att(ncid, tb_18h_varid, 'long_name', &
              'brightness temperature 18.85 GHz H-pol'), &
              '[ERR] nf90_put_att failed for TB_18H')
         call LDT_verify(nf90_put_att(ncid, tb_18h_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_18H units')

         ! TB_18V
         call LDT_verify(nf90_def_var(ncid, 'TB_18V', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_18v_varid), &
              '[ERR] nf90_def_var failed for TB_18V')
         call LDT_verify(nf90_def_var_fill(ncid, tb_18v_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_18V')
         call LDT_verify(nf90_put_att(ncid, tb_18v_varid, 'long_name', &
              'brightness temperature 18.85 GHz V-pol'), &
              '[ERR] nf90_put_att failed for TB_18V')
         call LDT_verify(nf90_put_att(ncid, tb_18v_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_18V units')

         ! TB_23H
         call LDT_verify(nf90_def_var(ncid, 'TB_23H', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_23h_varid), &
              '[ERR] nf90_def_var failed for TB_23H')
         call LDT_verify(nf90_def_var_fill(ncid, tb_23h_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_23H')
         call LDT_verify(nf90_put_att(ncid, tb_23h_varid, 'long_name', &
              'brightness temperature 23.80 GHz H-pol'), &
              '[ERR] nf90_put_att failed for TB_23H')
         call LDT_verify(nf90_put_att(ncid, tb_23h_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_23H units')

         ! TB_23V
         call LDT_verify(nf90_def_var(ncid, 'TB_23V', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_23v_varid), &
              '[ERR] nf90_def_var failed for TB_23V')
         call LDT_verify(nf90_def_var_fill(ncid, tb_23v_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_23V')
         call LDT_verify(nf90_put_att(ncid, tb_23v_varid, 'long_name', &
              'brightness temperature 23.80 GHz V-pol'), &
              '[ERR] nf90_put_att failed for TB_23V')
         call LDT_verify(nf90_put_att(ncid, tb_23v_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_23V units')

         ! TB_36H
         call LDT_verify(nf90_def_var(ncid, 'TB_36H', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_36h_varid), &
              '[ERR] nf90_def_var failed for TB_36H')
         call LDT_verify(nf90_def_var_fill(ncid, tb_36h_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_36H')
         call LDT_verify(nf90_put_att(ncid, tb_36h_varid, 'long_name', &
              'brightness temperature 36.75 GHz H-pol'), &
              '[ERR] nf90_put_att failed for TB_36H')
         call LDT_verify(nf90_put_att(ncid, tb_36h_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_36H units')

         ! TB_36V
         call LDT_verify(nf90_def_var(ncid, 'TB_36V', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_36v_varid), &
              '[ERR] nf90_def_var failed for TB_36V')
         call LDT_verify(nf90_def_var_fill(ncid, tb_36v_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_36V')
         call LDT_verify(nf90_put_att(ncid, tb_36v_varid, 'long_name', &
              'brightness temperature 36.75 GHz V-pol'), &
              '[ERR] nf90_put_att failed for TB_36V')
         call LDT_verify(nf90_put_att(ncid, tb_36v_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_36V units')

         ! TB_89H
         call LDT_verify(nf90_def_var(ncid, 'TB_89H', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_89h_varid), &
              '[ERR] nf90_def_var failed for TB_89H')
         call LDT_verify(nf90_def_var_fill(ncid, tb_89h_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_89H')
         call LDT_verify(nf90_put_att(ncid, tb_89h_varid, 'long_name', &
              'brightness temperature 89.00 GHz H-pol'), &
              '[ERR] nf90_put_att failed for TB_89H')
         call LDT_verify(nf90_put_att(ncid, tb_89h_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_89H units')

         ! TB_89V
         call LDT_verify(nf90_def_var(ncid, 'TB_89V', NF90_FLOAT, &
              [lon_dimid, lat_dimid], tb_89v_varid), &
              '[ERR] nf90_def_var failed for TB_89V')
         call LDT_verify(nf90_def_var_fill(ncid, tb_89v_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for TB_89V')
         call LDT_verify(nf90_put_att(ncid, tb_89v_varid, 'long_name', &
              'brightness temperature 89.00 GHz V-pol'), &
              '[ERR] nf90_put_att failed for TB_89V')
         call LDT_verify(nf90_put_att(ncid, tb_89v_varid, 'units', 'K'), &
              '[ERR] nf90_put_att failed for TB_89V units')

         ! Define auxiliary variables
         call LDT_verify(nf90_def_var(ncid, 'LAND_FRACTION', NF90_FLOAT, &
              [lon_dimid, lat_dimid], lwf_varid), &
              '[ERR] nf90_def_var failed for LAND_FRACTION')
         call LDT_verify(nf90_def_var_fill(ncid, lwf_varid, 0, -9999.0), &
              '[ERR] nf90_def_var_fill failed for LAND_FRACTION')
         call LDT_verify(nf90_put_att(ncid, lwf_varid, 'long_name', &
              'land fraction'), &
              '[ERR] nf90_put_att failed for LAND_FRACTION')

         call LDT_verify(nf90_def_var(ncid, 'QUALITY_FLAG', NF90_BYTE, &
              [lon_dimid, lat_dimid], qf_varid), &
              '[ERR] nf90_def_var failed for QUALITY_FLAG')
         call LDT_verify(nf90_put_att(ncid, qf_varid, 'long_name', &
              'quality flag'), &
              '[ERR] nf90_put_att failed for QUALITY_FLAG')
         call LDT_verify(nf90_put_att(ncid, qf_varid, 'flag_meanings', &
              'bit0=ocean bit1=precipitation bit2=snow bit3=sensor_quality'), &
              '[ERR] nf90_put_att failed for QUALITY_FLAG meanings')

         call LDT_verify(nf90_def_var(ncid, 'SAMPLE_V', NF90_INT, &
              [lon_dimid, lat_dimid], sv_varid), &
              '[ERR] nf90_def_var failed for SAMPLE_V')
         call LDT_verify(nf90_put_att(ncid, sv_varid, 'long_name', &
              'number of V-pol samples'), &
              '[ERR] nf90_put_att failed for SAMPLE_V')

         call LDT_verify(nf90_def_var(ncid, 'SAMPLE_H', NF90_INT, &
              [lon_dimid, lat_dimid], sh_varid), &
              '[ERR] nf90_def_var failed for SAMPLE_H')
         call LDT_verify(nf90_put_att(ncid, sh_varid, 'long_name', &
              'number of H-pol samples'), &
              '[ERR] nf90_put_att failed for SAMPLE_H')

         ! Add global attributes
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'title', &
              'WSF Low Resolution Resampled Brightness Temperatures'), &
              '[ERR] nf90_put_att failed for title')
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'source', &
              'WSF Microwave Sounder'), &
              '[ERR] nf90_put_att failed for source')
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'grid_resolution', &
              '0.140625 deg lon x 0.09375 deg lat'), &
              '[ERR] nf90_put_att failed for grid_resolution')

         ! End define mode
         call LDT_verify(nf90_enddef(ncid), &
              '[ERR] nf90_enddef failed')

         ! Write coordinate data
         allocate(lons_float(nc))
         allocate(lats_float(nr))

         do i = 1, nc
            lons_float(i) = real(lons(i), 4)
         end do
         do j = 1, nr
            lats_float(j) = real(lats(j), 4)
         end do

         call LDT_verify(nf90_put_var(ncid, lon_varid, lons_float), &
              '[ERR] nf90_put_var failed for lon')
         call LDT_verify(nf90_put_var(ncid, lat_varid, lats_float), &
              '[ERR] nf90_put_var failed for lat')
         call LDT_verify(nf90_put_var(ncid, time_varid, 0.0d0), &
              '[ERR] nf90_put_var failed for time')

         deallocate(lons_float, lats_float)

         ! Write TB data
         write(LDT_logunit,*)'[INFO] Writing brightness temperature data...'
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
         write(LDT_logunit,*)'[INFO] Writing auxiliary data...'
         call LDT_verify(nf90_put_var(ncid, lwf_varid, land_frac), &
              '[ERR] nf90_put_var failed for LAND_FRACTION')
         call LDT_verify(nf90_put_var(ncid, qf_varid, quality_flag), &
              '[ERR] nf90_put_var failed for QUALITY_FLAG')
         call LDT_verify(nf90_put_var(ncid, sv_varid, sample_v), &
              '[ERR] nf90_put_var failed for SAMPLE_V')
         call LDT_verify(nf90_put_var(ncid, sh_varid, sample_h), &
              '[ERR] nf90_put_var failed for SAMPLE_H')

         ! Close the file
         call LDT_verify(nf90_close(ncid), &
              '[ERR] nf90_close failed')

         write(LDT_logunit,*)'[INFO] ========================================'
         write(LDT_logunit,*)'[INFO] ✓ Successfully wrote NetCDF output file'
         write(LDT_logunit,*)'[INFO] ', trim(output_fname)
         write(LDT_logunit,*)'[INFO] ========================================'
      endif

      ierr = 0
#if (defined SPMD)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
    end subroutine LDT_WSF_ARFS_write_netcdf
#else
  ! Dummy version
  subroutine LDT_WSF_ARFS_write_netcdf(nc, nr, &
       tb_10h, tb_10v, tb_18h, tb_18v, tb_23h, tb_23v, &
       tb_36h, tb_36v, tb_89h, tb_89v, land_frac, quality_flag, &
       sample_v, sample_h, lats, lons, output_fname)
    use LDT_logMod, only: LDT_logunit, LDT_endrun
    implicit none
    integer, intent(in) :: nc, nr
    real*4, intent(in) :: tb_10h(nc,nr), tb_10v(nc,nr)
    real*4, intent(in) :: tb_18h(nc,nr), tb_18v(nc,nr)
    real*4, intent(in) :: tb_23h(nc,nr), tb_23v(nc,nr)
    real*4, intent(in) :: tb_36h(nc,nr), tb_36v(nc,nr)
    real*4, intent(in) :: tb_89h(nc,nr), tb_89v(nc,nr)
    real*4, intent(in) :: land_frac(nc,nr)
    integer*1, intent(in) :: quality_flag(nc,nr)
    integer*4, intent(in) :: sample_v(nc,nr), sample_h(nc,nr)
    real*8, intent(in) :: lats(nr), lons(nc)
    character(*), intent(in) :: output_fname
    write(LDT_logunit,*)'[ERR] LDT not compiled with NETCDF support!'
    write(LDT_logunit,*)'Cannot write WSF ARFS resampled data in NETCDF format!'
    write(LDT_logunit,*)'Recompile with NETCDF support and try again!'
    call LDT_endrun()
  end subroutine LDT_WSF_ARFS_write_netcdf
#endif

end module LDT_WSF_ARFS_netcdfMod