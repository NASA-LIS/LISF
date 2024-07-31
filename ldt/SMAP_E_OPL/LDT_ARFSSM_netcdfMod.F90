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
! MODULE: LDT_ARFSSM_netcdfMod
!
! REVISION HISTORY:
!  10 Feb 2023: Eric Kemp.  Initial implementation.
!
! DESCRIPTION: Outputs SMAP based soil moisture retrieval in netCDF format.
!
!------------------------------------------------------------------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_ARFSSM_netcdfMod

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: LDT_ARFSSM_write_netcdf

contains

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ! Subroutine for writing ARFS SM retrieval to netCDF
  subroutine LDT_ARFSSM_write_netcdf(nc, nr, arfs_sm, retrieval_fname, &
       yyyymmdd, hhmmss)

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
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real*4, intent(in) :: arfs_sm(nc,nr)
    character(*), intent(in) :: retrieval_fname
    character(8), intent(in) :: yyyymmdd
    character(6), intent(in) :: hhmmss

    ! Locals
    integer :: shuffle, deflate, deflate_level
    integer :: iret
    integer :: ncid
    integer :: dim_ids(3)
    real :: dlat, dlon
    real :: swlat, swlon
    real :: nelat, nelon
    integer :: lon_varid, lat_varid, time_varid, sm_varid
    character(120) :: time_units
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    integer :: values(8)
    real, allocatable :: lats(:)
    real, allocatable :: lons(:)
    integer :: ierr
    integer :: i, j

    ! Only the master process handles the file output
    if (LDT_masterproc) then
       write(LDT_logunit,*)'[INFO] Creating NETCDF file ', &
            trim(retrieval_fname)

       ! Copy netCDF compression settings
       shuffle = NETCDF_shuffle
       deflate = NETCDF_deflate
       deflate_level = NETCDF_deflate_level

       ! Create the output file
#if(defined USE_NETCDF3)
         iret=nf90_create(path=trim(retrieval_fname),&
              cmode=NF90_CLOBBER, ncid=ncid)
         call LDT_verify(iret,&
              '[ERR] nf90_create failed')
#endif
#if(defined USE_NETCDF4)
         iret=nf90_create(path=trim(retrieval_fname),&
              cmode=NF90_NETCDF4, ncid=ncid)
         call LDT_verify(iret, &
              '[ERR] nf90_create failed')
#endif

         ! Write out dimensions headers
         call LDT_verify(nf90_def_dim(ncid, 'time', 1, dim_ids(3)), &
              '[ERR] nf90_def_dim failed')
         call LDT_verify(nf90_def_dim(ncid, 'lat', nr, dim_ids(2)), &
              '[ERR] nf90_def_dim failed')
         call LDT_verify(nf90_def_dim(ncid, 'lon', nc, dim_ids(1)), &
              '[ERR] nf90_def_dim failed')

         ! Map projection
         !FIXME: Allow map projections other than lat/lon
         select case("latlon")
         case ("latlon")
            dlon = LDT_rc%gridDesc(1,9)
            dlat = LDT_rc%gridDesc(1,10)
            swlat = LDT_rc%gridDesc(1,4)
            swlon = LDT_rc%gridDesc(1,5)
            nelat = LDT_rc%gridDesc(1,7)
            nelon = LDT_rc%gridDesc(1,8)

            call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL,&
                 "MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL"), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL,&
                 "SOUTH_WEST_CORNER_LAT", swlat), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL,&
                 "SOUTH_WEST_CORNER_LON", swlon), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, &
                 "NORTH_EAST_CORNER_LAT", nelat), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, &
                 "NORTH_EAST_CORNER_LON", nelon), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, &
                 "DX", dlon),&
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, &
                 "DY", dlat), &
                 '[ERR] nf90_put_att failed')

         case default
            write(LDT_logunit,*) &
                 '[ERR] Only latlon map projection supported for SMAP_E_OPL'
            call LDT_endrun()
         end select

         ! Include the water points
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, &
              "INC_WATER_PTS", "true"), &
              '[ERR] nf90_put_att failed')

         ! Construct the longitudes
         ! FIXME:  Add support for other map projections
         call LDT_verify(nf90_def_var(ncid, "lon", NF90_FLOAT, dim_ids(1), &
              lon_varid),'[ERR] nf90_def_var failed')
#if(defined USE_NETCDF4)
         call LDT_verify(nf90_def_var_deflate(ncid, &
              lon_varid, shuffle, deflate, deflate_level),&
              '[ERR] nf90_def_var_deflate')
#endif
         call LDT_verify(nf90_put_att(ncid, lon_varid, &
              "units", "degrees_east"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, lon_varid, &
              "long_name", "longitude"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, lon_varid, &
              "standard_name", "longitude"), &
              '[ERR] nf90_put_att failed')

         ! Construct the latitudes
         ! FIXME:  Add support for other map projections
         call LDT_verify(nf90_def_var(ncid, "lat", NF90_FLOAT, dim_ids(2), &
              lat_varid), '[ERR] nf90_def_var failed')
#if(defined USE_NETCDF4)
         call LDT_verify(nf90_def_var_deflate(ncid, &
              lat_varid, shuffle, deflate, deflate_level),&
              '[ERR] nf90_def_var_deflate')
         call LDT_verify(nf90_put_att(ncid, lat_varid, &
              "units", "degrees_north"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, lat_varid, &
              "long_name", "latitude"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, lat_varid, &
              "standard_name", "latitude"), &
              '[ERR] nf90_put_att failed')
#endif

         ! Define the time array.  The valid time will be written as an
         ! attribute.
         call LDT_verify(nf90_def_var(ncid, 'time', NF90_DOUBLE, &
              dimids=dim_ids(3), varid=time_varid), &
              '[ERR] nf90_def_var failed')
         write(time_units,'(A)') &
              "seconds since "//yyyymmdd(1:4)//"-" &
              //yyyymmdd(5:6)//"-" &
              //yyyymmdd(7:8)//" " &
              //hhmmss(1:2)//":"//hhmmss(3:4)//":"//hhmmss(5:6)
         call LDT_verify(nf90_put_att(ncid, time_varid, &
              "units", trim(time_units)), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, time_varid, &
              "long_name", "time"), &
              '[ERR] nf90_put_att failed')

         ! Define the soil moisture retrieval
         call LDT_verify(nf90_def_var(ncid, "arfs_sm", NF90_FLOAT, &
              dimids=dim_ids, &
              varid=sm_varid), '[ERR] nf90_def_var failed')
#if (defined USE_NETCDF4)
         call LDT_verify(nf90_def_var_deflate(ncid, &
              sm_varid, shuffle, deflate, deflate_level), &
              '[ERR] nf90_def_var_deflate failed')
#endif
         call LDT_verify(nf90_put_att(ncid, sm_varid, &
              "units", "m^3/m^3"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, sm_varid, &
              "long_name","soil moisture"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, sm_varid, &
              '_FillValue', -9999.), &
              '[ERR] nf90_put_att failed')

         ! Miscellaneous header information
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", &
              "CF-1.10"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, "title", &
              "LDT SMAP_E_OPL retrieval"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, "institution", &
              "NASA GSFC Hydrological Sciences Laboratory"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, "source", &
              "Land Data Toolkit (LDT)"), &
              '[ERR] nf90_put_att failed')
         call date_and_time(date, time, zone, values)
#ifndef LDT_SKIP_HISTORY
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, "history", &
              "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
              date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)), &
              '[ERR] nf90_put_att failed')
#endif
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, "references", &
              "Arsenault_etal_GMD_2018, Kumar_etal_EMS_2006"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, "comment", &
              "website: http://lis.gsfc.nasa.gov/"), &
              '[ERR] nf90_put_att failed')

         ! We are ready to write the actual data.  This requires taking NETCDF
         ! out of define mode.
         call LDT_verify(nf90_enddef(ncid), &
              '[ERR] ncf90_enddef failed')

         ! Write the latitude data
         allocate(lats(nr))
         do j = 1, nr
            lats(j) = swlat + (j-1)*(dlat)
         end do
         call LDT_verify(nf90_put_var(ncid, lat_varid, &
              lats(:), (/1/), (/nr/)), &
              '[ERR] nf90_put_var failed for lats')
         deallocate(lats)

         ! Write the longitude data
         allocate(lons(nc))
         do i = 1, nc
            lons(i) = swlon + (i-1)*(dlon)
         end do
         call LDT_verify(nf90_put_var(ncid, lon_varid, lons(:),&
              (/1/), (/nc/)), &
              '[ERR] nf90_put_var failed for lon')
         deallocate(lons)

         ! Write the time data
         call LDT_verify(nf90_put_var(ncid, time_varid, 0.0), &
              '[ERR] nf90_put_var failed for time')

         ! Write the ARFS SM field
         call LDT_verify(nf90_put_var(ncid, sm_varid, &
              arfs_sm(:,:), &
              (/1,1/), (/nc,nr/)), &
              '[ERR] nf90_put_var failed for sm')

         ! Close the file and clean up
         call LDT_verify(nf90_close(ncid), &
              '[ERR] nf90_close failed!')
      endif

      ierr = 0
#if (defined SPMD)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif
    end subroutine LDT_ARFSSM_write_netcdf
#else
  ! Dummy version
  subroutine LDT_ARFSSM_write_netcdf(nc, nr, arfs_sm, retrieval_fname, &
       yyyymmdd, hhmmss)
    use LDT_logMod, only: LDT_logunit, LDT_endrun
    implicit none
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real*4, intent(in) :: arfs_sm(nc,nr)
    character(*), intent(in) :: retrieval_fname
    character(8), intent(in) :: yyyymmdd
    character(6), intent(in) :: hhmmss
    write(LDT_logunit,*)'[ERR] LDT not compiled with NETCDF support!'
    write(LDT_logunit,*)'Cannot write ARFS SM retrieval in NETCDF format!'
    write(LDT_logunit,*)'Recompile with NETCDF support and try again!'
    call LDT_endrun()
  end subroutine LDT_ARFSSM_write_netcdf
#endif

end module LDT_ARFSSM_netcdfMod

