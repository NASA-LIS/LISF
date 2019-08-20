!-----------------------Begin NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LDTSI_netcdfMod
! 
! REVISION HISTORY:
! 01 Mar 2019  Eric Kemp  First version.
!
! DESCRIPTION:
! Source code for writing Air Force snow depth analysis variables to 
! NETCDF.
!-------------------------------------------------------------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDTSI_netcdfMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: LDTSI_write_netcdf
   public :: LDTSI_read_netcdf
   public :: LDTSI_read_netcdf_12z

contains

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   ! Subroutine for writing LDTSI analysis to netCDF
   subroutine LDTSI_write_netcdf(date10)

      ! Imports
      use LDT_coreMod, only: LDT_rc, LDT_masterproc
      use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
      use LDT_ldtsiMod, only: ldtsi_settings
#if ( defined SPMD )
      use mpi
#endif
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif
      use LDTSI_arraysMod, only: LDTSI_arrays
      use LDTSI_paramsMod

      ! Arguments
      character*10, intent(in) :: date10

      ! Local variables
      character*255 :: outfilename
      integer :: iret
      integer :: ncid
      character*8 :: date
      character*10 :: time
      character*5 :: zone
      integer :: values(8)
      integer :: dim_ids(3)
      integer :: snoanl_varid, snoage_varid, icecon_varid, icemask_varid, &
           iceage_varid
      integer :: lon_varid, lat_varid, time_varid
      integer               :: shuffle, deflate, deflate_level
      real, allocatable :: lats(:)
      real, allocatable :: lons(:)
      real, allocatable :: icecon_tmp(:,:)
      integer :: i,j
      integer :: ierr
      character*4 :: cyyyy
      character*2 :: cmm,cdd,chh
      character*120 :: time_units
      integer :: nc,nr
      real :: dlat, dlon
      real :: swlat, swlon
      real :: nelat, nelon

      ! Only the master process handles the file output
      if (LDT_masterproc) then

         nc = LDT_rc%lnc(1)
         nr = LDT_rc%lnr(1)

         ! FIXME:  Set this in ldt.config
         outfilename = "ldtsi_"//date10//".nc"

         write(LDT_logunit,*)'[INFO] Creating NETCDF file ',trim(outfilename)

         ! Copy netcdf compression settings
         shuffle = NETCDF_shuffle
         deflate = NETCDF_deflate
         deflate_level = NETCDF_deflate_level

         ! Create the output file
#if(defined USE_NETCDF3) 
         iret=nf90_create(path=trim(outfilename),&
              cmode=nf90_clobber, ncid=ncid)
         call LDT_verify(iret,&
              '[ERR] nf90_create failed')
#endif
#if(defined USE_NETCDF4) 
         iret=nf90_create(path=trim(outfilename),&
              cmode=nf90_netcdf4, ncid=ncid)
         call LDT_verify(iret, &
              '[ERR] nf90_create failed')
#endif    
         
         ! Write out dimensions headers
         call LDT_verify(nf90_def_dim(ncid,'time',1,dim_ids(3)), &
              '[ERR] nf90_def_dim failed')
         call LDT_verify(nf90_def_dim(ncid,'lat',nr,dim_ids(2)), &
              '[ERR] nf90_def_dim failed')
         call LDT_verify(nf90_def_dim(ncid,'lon',nc,dim_ids(1)), &
              '[ERR] nf90_def_dim failed')
         
         ! Map projection
         !FIXME:  Allow map projections other than lat/lon
         select case("latlon")
         case ("latlon")

            dlon = LDT_rc%gridDesc(1,9)
            dlat = LDT_rc%gridDesc(1,10)
            swlat = LDT_rc%gridDesc(1,4)
            swlon = LDT_rc%gridDesc(1,5)
            nelat = LDT_rc%gridDesc(1,7)
            nelon = LDT_rc%gridDesc(1,8)

            call LDT_verify(nf90_put_att(ncid,nf90_global,&
                 "MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL"), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid,nf90_global,&
                 "SOUTH_WEST_CORNER_LAT", swlat), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid,nf90_global,&
                 "SOUTH_WEST_CORNER_LON", swlon), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid,nf90_global, &
                 "NORTH_EAST_CORNER_LAT", nelat), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid,nf90_global, &
                 "NORTH_EAST_CORNER_LON", nelon), &
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid,nf90_global, &
                 "DX",dlon),&
                 '[ERR] nf90_put_att failed')
            call LDT_verify(nf90_put_att(ncid,nf90_global, &
                 "DY",dlat), &
                 '[ERR] nf90_put_att failed')
            
         case default
            write(LDT_logunit,*) &
                 '[ERR] Only latlon map projection supported for LDTSI'
            call LDT_endrun()
         end select
         
         ! Include the water points
         call LDT_verify(nf90_put_att(ncid,nf90_global, &
              "INC_WATER_PTS","true"), &
              '[ERR] nf90_put_att failed')

         ! Construct the longitudes
         ! FIXME:  Add support for other map projections
         call LDT_verify(nf90_def_var(ncid,"lon",nf90_float,dim_ids(1), &
              lon_varid),'[ERR] nf90_def_var failed')
#if(defined USE_NETCDF4) 
         call LDT_verify(nf90_def_var_deflate(ncid,&
              lon_varid, shuffle, deflate, deflate_level),&
              '[ERR] nf90_def_var_deflate')
#endif
         call LDT_verify(nf90_put_att(ncid,lon_varid, &
              "units","degrees_east"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,lon_varid, &
              "long_name","longitude"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,lon_varid, &
              "standard_name","longitude"),&
              '[ERR] nf90_put_att failed')

         ! Construct the latitudes
         ! FIXME:  Add support for other map projections
         call LDT_verify(nf90_def_var(ncid,"lat",nf90_float,dim_ids(2), &
              lat_varid),'[ERR] nf90_def_var failed')
#if(defined USE_NETCDF4) 
         call LDT_verify(nf90_def_var_deflate(ncid,&
              lat_varid, shuffle, deflate, deflate_level),&
              '[ERR] nf90_def_var_deflate')
#endif
         call LDT_verify(nf90_put_att(ncid,lat_varid, &
              "units","degrees_north"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,lat_varid, &
              "long_name","latitude"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,lat_varid, &
              "standard_name","latitude"),&
              '[ERR] nf90_put_att failed')

         ! Define the time array.  The valid time will be written as an
         ! attribute.
         call LDT_verify(nf90_def_var(ncid,'time',nf90_double,&
              dimids=dim_ids(3), varid=time_varid), &
              '[ERR] nf90_def_var failed')
         cyyyy = ldtsi_settings%date10(1:4)
         cmm = ldtsi_settings%date10(5:6)
         cdd = ldtsi_settings%date10(7:8)
         chh = ldtsi_settings%date10(9:10)
         write(time_units,'(A)') &
              "seconds since "//trim(cyyyy)//"-"//trim(cmm)//"-"//trim(cdd)//&
              " "//trim(chh)//":00:00"
         call LDT_verify(nf90_put_att(ncid,time_varid, &
              "units",trim(time_units)),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,time_varid, &
              "long_name","time"),&             
              '[ERR] nf90_put_att failed')
         
         ! Define the snow depth analysis
         call LDT_verify(nf90_def_var(ncid,"snoanl",nf90_float, &
              dimids=dim_ids, &
              varid=snoanl_varid),'[ERR] nf90_def_var failed')
#if (defined USE_NETCDF4) 
         call LDT_verify(nf90_def_var_deflate(ncid,&
              snoanl_varid, shuffle, deflate, deflate_level), &
              '[ERR] nf90_def_var_deflate failed')
#endif
         call LDT_verify(nf90_put_att(ncid,snoanl_varid, &
              "units","m"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,snoanl_varid, &
              "long_name","depth of surface snow over land"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,snoanl_varid, &
              '_FillValue',-1.), &
              '[ERR] nf90_put_att failed for SNOANL')

         ! Define the snow age analysis
         call LDT_verify(nf90_def_var(ncid,"snoage",nf90_float, &
              dimids=dim_ids, &
              varid=snoage_varid),'[ERR] nf90_def_var failed')
#if (defined USE_NETCDF4) 
         call LDT_verify(nf90_def_var_deflate(ncid,&
              snoage_varid, shuffle, deflate, deflate_level), &
              '[ERR] nf90_def_var_deflate failed')
#endif
         call LDT_verify(nf90_put_att(ncid,snoage_varid, &
              "units","days"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,snoage_varid, &
              "long_name","age of surface snow over land"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,snoage_varid, &
              "standard_name","age_of_surface_snow"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,snoage_varid, &
              '_FillValue',-1.), &
              '[ERR] nf90_put_att failed for SNOAGE')

         ! Define the ice concentration analysis
         call LDT_verify(nf90_def_var(ncid,"icecon",nf90_float, &
              dimids=dim_ids, &
              varid=icecon_varid),'[ERR] nf90_def_var failed')
#if (defined USE_NETCDF4) 
         call LDT_verify(nf90_def_var_deflate(ncid,&
              icecon_varid, shuffle, deflate, deflate_level), &
              '[ERR] nf90_def_var_deflate failed')
#endif
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              "units","1"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              "long_name","concentration of sea ice (0-1)"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              "standard_name","sea_ice_area_fraction"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              "scale_factor",100.),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              "add_offset",0.),&
              '[ERR] nf90_put_att failed')         
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              'missing_value',-1.), &
              '[ERR] nf90_put_att failed for ICECON')
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              '_FillValue',-1.), &
              '[ERR] nf90_put_att failed for ICECON')
         call LDT_verify(nf90_put_att(ncid,icecon_varid, &
              'valid_range',(/0.,100./)), &
              '[ERR] nf90_put_att failed for ICECON')

         ! Define the ice mask analysis
         call LDT_verify(nf90_def_var(ncid,"icemask",nf90_float, &
              dimids=dim_ids, &
              varid=icemask_varid),'[ERR] nf90_def_var failed')
#if (defined USE_NETCDF4) 
         call LDT_verify(nf90_def_var_deflate(ncid,&
              icemask_varid, shuffle, deflate, deflate_level), &
              '[ERR] nf90_def_var_deflate failed')
#endif
         call LDT_verify(nf90_put_att(ncid,icemask_varid, &
              "units","1"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,icemask_varid, &
              "long_name","sea ice mask (1=ice,0=no ice)"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,icemask_varid, &
              '_FillValue',-1.), &
              '[ERR] nf90_put_att failed for ICEMASK')

         ! Define the ice age analysis
         call LDT_verify(nf90_def_var(ncid,"iceage",nf90_float, &
              dimids=dim_ids, &
              varid=iceage_varid),'[ERR] nf90_def_var failed')      
#if (defined USE_NETCDF4) 
         call LDT_verify(nf90_def_var_deflate(ncid,&
              iceage_varid, shuffle, deflate, deflate_level), &
              '[ERR] nf90_def_var_deflate failed')
#endif
         call LDT_verify(nf90_put_att(ncid,iceage_varid, &
              "units","days"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,iceage_varid, &
              "long_name","age of sea ice"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,iceage_varid, &
              "standard_name","age_of_sea_ice"),&
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,iceage_varid, &
              '_FillValue',-1.), &
              '[ERR] nf90_put_att failed for ICEAGE')

         ! Miscellaneous header information
         call LDT_verify(nf90_put_att(ncid,nf90_global,"Conventions", &
              "CF-1.7"), &
              '[ERR] nf90_put_att failed')         
         call LDT_verify(nf90_put_att(ncid,nf90_global,"title", &
              "LDT LDTSI analysis"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,nf90_global,"institution", &
              "NASA GSFC Hydrological Sciences Laboratory"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,nf90_global,"source", &
              "Land Data Toolkit (LDT)"), &
              '[ERR] nf90_put_att failed')
         call date_and_time(date,time,zone,values)
#ifndef LDT_SKIP_HISTORY
         call LDT_verify(nf90_put_att(ncid,nf90_global,"history", &
              "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
              date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)), &
              '[ERR] nf90_put_att failed')
#endif
         call LDT_verify(nf90_put_att(ncid,nf90_global,"references", &
              "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"), &
              '[ERR] nf90_put_att failed')
         call LDT_verify(nf90_put_att(ncid,nf90_global,"comment", &
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
         call LDT_verify(nf90_put_var(ncid,lat_varid, &
              lats(:),(/1/),(/nr/)), &
              '[ERR] nf90_put_var failed for lats')
         deallocate(lats)

         ! Write the longitude data
         allocate(lons(nc))
         do i = 1, nc
            lons(i) = swlon + (i-1)*(dlon)
         end do
         call LDT_verify(nf90_put_var(ncid,lon_varid,lons(:),&
              (/1/),(/nc/)), &
              '[ERR] nf90_put_var failed for lon')
         deallocate(lons)

         ! Write the time data
         call LDT_verify(nf90_put_var(ncid,time_varid,0.0), &
              '[ERR] nf90_put_var failed for time')
         
         ! Write the LDTSI fields
         call LDT_verify(nf90_put_var(ncid,snoanl_varid,&
              LDTSI_arrays%snoanl(:,:), &
              (/1,1/),(/nc,nr/)), &
              '[ERR] nf90_put_var failed for snoanl')
         call LDT_verify(nf90_put_var(ncid,snoage_varid,&
              LDTSI_arrays%snoage(:,:), &
              (/1,1/),(/nc,nr/)), &
              '[ERR] nf90_put_var failed')

         ! Special handling of ice concentration.  Since this is a standard
         ! variable with range 0-1, we compress the original values (0-100),
         ! but taking care not to change the missing points.
         allocate(icecon_tmp(nc,nr))
         do j = 1,nr
            do i = 1,nc
               if (LDTSI_arrays%icecon(i,j) < 0) then
                  icecon_tmp(i,j) = -1
               else
                  icecon_tmp(i,j) = 0.01*(LDTSI_arrays%icecon(i,j))
               endif
            end do
         end do
         call LDT_verify(nf90_put_var(ncid,icecon_varid,icecon_tmp(:,:), &
              (/1,1/),(/nc,nr/)), &
              '[ERR] nf90_put_var failed for icecon')
         deallocate(icecon_tmp)

         ! The rest of the LDTSI fields.
         call LDT_verify(nf90_put_var(ncid,icemask_varid, &
              LDTSI_arrays%icemask(:,:), &
              (/1,1/),(/nc,nr/)), &
              '[ERR] nf90_put_var failed for icemask')
         call LDT_verify(nf90_put_var(ncid,iceage_varid, &
              LDTSI_arrays%iceage(:,:), &
              (/1,1/),(/nc,nr/)), &
              '[ERR] nf90_put_var failed for iceage')
                  
         ! Close the file and clean up
         call LDT_verify(nf90_close(ncid), &
              '[ERR] nf90_close failed!')
      end if

      ierr = 0
#if (defined SPMD)
      call mpi_barrier(mpi_comm_world,ierr)
#endif

      write(LDT_logunit,*) &
           '[INFO] Finished writing LDT LDTSI parameters to netcdf file'

   end subroutine LDTSI_write_netcdf
      
#else
   ! Dummy version
   subroutine LDTSI_write_netcdf(date10)
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      implicit none
      character*10, intent(in) :: date10
      write(LDT_logunit,*)'[ERR] LDT not compiled with NETCDF support!'
      write(LDT_logunit,*)'Cannot write out LDTSI data'
      write(LDT_logunit,*)'Recompile with NETCDF support and try again!'
      call LDT_endrun()
   end subroutine LDTSI_write_netcdf
#endif

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   ! Subroutine for reading LDTSI analysis to netCDF
   subroutine LDTSI_read_netcdf(date10,ierr)

      ! Imports
      use LDT_coreMod, only: LDT_rc
      use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif
      use LDTSI_arraysMod, only: LDTSI_arrays

      ! Defaults
      implicit none

      ! Arguments
      character*10, intent(in) :: date10
      integer, intent(out) :: ierr

      ! Local variables
      character*255 :: infilename
      logical :: file_exists
      integer :: ncid, dim_ids(3)
      integer :: snoanl_varid, snoage_varid, icecon_varid, icemask_varid, &
           iceage_varid
      integer :: nc,nr,ntime,nlat,nlon
      integer :: c,r
      real, allocatable :: tmp(:,:,:)

      ierr = 1 ! Update below

      nc = LDT_rc%lnc(1)
      nr = LDT_rc%lnr(1)

      ! See if file exists
      infilename = "ldtsi_"//date10//".nc"
      inquire(file=trim(infilename), exist=file_exists)
      if (.not. file_exists) return

      write(LDT_logunit,*)'[INFO] Reading prior NETCDF file ',trim(infilename)

      ! Open the file for reading
      call LDT_verify(nf90_open(path=trim(infilename), &
           mode=NF90_NOWRITE, &
           ncid=ncid), &
           '[ERR] Error in nf90_open for '//trim(infilename))

      ! Get the dimension IDs
      call LDT_verify(nf90_inq_dimid(ncid=ncid,&
           name='time',&
           dimid=dim_ids(3)), &
           '[ERR] Error in nf90_inq_dimid for dimension time')
      call LDT_verify(nf90_inq_dimid(ncid=ncid,&
           name='lat',&
           dimid=dim_ids(2)), &
           '[ERR] Error in nf90_inq_dimid for dimension lat')
      call LDT_verify(nf90_inq_dimid(ncid=ncid,&
           name='lon',&
           dimid=dim_ids(1)), &
           '[ERR] Error in nf90_inq_dimid for dimension lon')

      ! Get the actual dimension sizes
      call LDT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(3), &
           len=ntime), &
           '[ERR] Error in nf90_inquire_dimension for dimension time')
      call LDT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(2), &
           len=nlat), &
           '[ERR] Error in nf90_inquire_dimension for dimension lat')
      call LDT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(1), &
           len=nlon), &
           '[ERR] Error in nf90_inquire_dimension for dimension lon')
      
      ! Sanity checks
      if (ntime .ne. 1 .or. &
           nlat .ne. nr .or. &
           nlon .ne. nc) then
         write(LDT_logunit, *) &
              '[ERR] Dimension mismatch between LDT and LDTSI input'
         write(LDT_logunit,*) &
              '[ERR] Expected time = 1, lat = ',nr,' lon = ',nc
         write(LDT_logunit,*) &
              '[ERR] Found time = ',ntime,' lat = ',nlat,' lon = ',nlon
         call LDT_endrun()
      end if

      ! Fetch the LDTSI analysis variable IDs
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='snoanl', &
           varid=snoanl_varid), &
           '[ERR] Error in nf90_inq_varid for snoanl')
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='snoage', &
           varid=snoage_varid), &
           '[ERR] Error in nf90_inq_varid for snoage')
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='icecon', &
           varid=icecon_varid), &
           '[ERR] Error in nf90_inq_varid for icecon')
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='icemask', &
           varid=icemask_varid), &
           '[ERR] Error in nf90_inq_varid for icemask')
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='iceage', &
           varid=iceage_varid), &
           '[ERR] Error in nf90_inq_varid for iceage')

      ! Read the LDTSI variables
      allocate(tmp(nc,nr,1)) ! Need 3D array
      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=snoanl_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoanl')
      LDTSI_arrays%olddep(:,:) = tmp(:,:,1)

      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=snoage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoage')
      !LDTSI_arrays%snoage(:,:) = tmp(:,:,1)
      do r = 1, nr
         do c = 1, nc
            LDTSI_arrays%snoage(c,r) = nint(tmp(c,r,1))
         end do ! c
      end do ! r

      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=icecon_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icecon')
      ! A unit transform is needed here
      do r = 1, nr
         do c = 1, nc
            if (tmp(c,r,1) < 0) then
               LDTSI_arrays%oldcon(c,r) = -1
            else
               LDTSI_arrays%oldcon(c,r) = nint(100*tmp(c,r,1))
            end if
         end do
      end do

      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=icemask_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for icemask')
      !LDTSI_arrays%oldmask(:,:) = tmp(:,:,1)
      do r = 1, nr
         do c = 1, nc
            LDTSI_arrays%oldmask(c,r) = nint(tmp(c,r,1))
         end do ! c
      end do ! r
 
      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=iceage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for iceage')
      !LDTSI_arrays%iceage(:,:) = tmp(:,:,1)
      do r = 1, nr
         do c = 1, nc
            LDTSI_arrays%iceage(c,r) = nint(tmp(c,r,1))
         end do ! c
      end do ! r

      deallocate(tmp)

      ! Close the file
      call LDT_verify(nf90_close(ncid), &
           '[ERR] Error in nf90_close for '//trim(infilename))

      ! Normal exit
      ierr = 0
   end subroutine LDTSI_read_netcdf

#else
   ! Dummy version
   subroutine LDTSI_read_netcdf(date10, ierr)
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      implicit none
      character*10, intent(in) :: date10
      integer, intent(out) :: ierr
      ierr = 1
      write(LDT_logunit,*)'[ERR] LDT not compiled with NETCDF support!'
      write(LDT_logunit,*)'Cannot read out LDTSI data'
      write(LDT_logunit,*)'Recompile with NETCDF support and try again!'
      call LDT_endrun()
   end subroutine LDTSI_read_netcdf

#endif

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   ! Subroutine for reading 12Z LDTSI analysis to netCDF
   ! Only snoage and iceage are read.
   subroutine LDTSI_read_netcdf_12z(date10,ierr)

      ! Imports
      use LDT_coreMod, only: LDT_rc
      use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
      use netcdf
#endif
      use LDTSI_arraysMod, only: LDTSI_arrays

      ! Defaults
      implicit none

      ! Arguments
      character*10, intent(in) :: date10
      integer, intent(out) :: ierr

      ! Local variables
      character*255 :: infilename
      logical :: file_exists
      integer :: ncid, dim_ids(3)
      integer :: snoage_varid, iceage_varid
      integer :: nc,nr,ntime,nlat,nlon
      integer :: c,r
      real, allocatable :: tmp(:,:,:)

      ierr = 1 ! Update below

      nc = LDT_rc%lnc(1)
      nr = LDT_rc%lnr(1)

      ! See if file exists
      infilename = "ldtsi_"//date10//".nc"
      inquire(file=trim(infilename), exist=file_exists)
      if (.not. file_exists) return

      write(LDT_logunit,*)'[INFO] Reading prior 12Z NETCDF file ', &
           trim(infilename)

      ! Open the file for reading
      call LDT_verify(nf90_open(path=trim(infilename), &
           mode=NF90_NOWRITE, &
           ncid=ncid), &
           '[ERR] Error in nf90_open for '//trim(infilename))

      ! Get the dimension IDs
      call LDT_verify(nf90_inq_dimid(ncid=ncid,&
           name='time',&
           dimid=dim_ids(3)), &
           '[ERR] Error in nf90_inq_dimid for dimension time')
      call LDT_verify(nf90_inq_dimid(ncid=ncid,&
           name='lat',&
           dimid=dim_ids(2)), &
           '[ERR] Error in nf90_inq_dimid for dimension lat')
      call LDT_verify(nf90_inq_dimid(ncid=ncid,&
           name='lon',&
           dimid=dim_ids(1)), &
           '[ERR] Error in nf90_inq_dimid for dimension lon')

      ! Get the actual dimension sizes
      call LDT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(3), &
           len=ntime), &
           '[ERR] Error in nf90_inquire_dimension for dimension time')
      call LDT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(2), &
           len=nlat), &
           '[ERR] Error in nf90_inquire_dimension for dimension lat')
      call LDT_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(1), &
           len=nlon), &
           '[ERR] Error in nf90_inquire_dimension for dimension lon')
      
      ! Sanity checks
      if (ntime .ne. 1 .or. &
           nlat .ne. nr .or. &
           nlon .ne. nc) then
         write(LDT_logunit, *) &
              '[ERR] Dimension mismatch between LDT and LDTSI input'
         write(LDT_logunit,*) &
              '[ERR] Expected time = 1, lat = ',nr,' lon = ',nc
         write(LDT_logunit,*) &
              '[ERR] Found time = ',ntime,' lat = ',nlat,' lon = ',nlon
         call LDT_endrun()
      end if

      ! Fetch the LDTSI analysis variable IDs
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='snoage', &
           varid=snoage_varid), &
           '[ERR] Error in nf90_inq_varid for snoage')
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='iceage', &
           varid=iceage_varid), &
           '[ERR] Error in nf90_inq_varid for iceage')

      ! Read the LDTSI variables
      allocate(tmp(nc,nr,1)) ! Need 3D array
      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=snoage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoage')
      !LDTSI_arrays%snoage12z(:,:) = tmp(:,:,1)
      do r = 1, nr
         do c = 1, nc
            LDTSI_arrays%snoage12z(c,r) = nint(tmp(c,r,1))            
         end do ! c
      end do ! r

      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=iceage_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for iceage')
      !LDTSI_arrays%iceage12z(:,:) = tmp(:,:,1)
      do r = 1, nr
         do c = 1, nc
            LDTSI_arrays%iceage12z(c,r) = nint(tmp(c,r,1))
         end do ! c
      end do ! r

      deallocate(tmp)

      ! Close the file
      call LDT_verify(nf90_close(ncid), &
           '[ERR] Error in nf90_close for '//trim(infilename))

      ! Normal exit
      ierr = 0
   end subroutine LDTSI_read_netcdf_12z

#else
   ! Dummy version
   subroutine LDTSI_read_netcdf_12z(date10,ierr)
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      implicit none
      character*10, intent(in) :: date10
      integer, intent(out) :: ierr
      ierr = 1
      write(LDT_logunit,*)'[ERR] LDT not compiled with NETCDF support!'
      write(LDT_logunit,*)'Cannot read out LDTSI data'
      write(LDT_logunit,*)'Recompile with NETCDF support and try again!'
      call LDT_endrun()
   end subroutine LDTSI_read_netcdf_12z

#endif

end module LDTSI_netcdfMod
