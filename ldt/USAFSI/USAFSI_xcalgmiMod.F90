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
! MODULE: USAFSI_xcalgmiMod
!
! REVISION HISTORY:
!  16 Nov 2020: Yonghwan Kwon; Initial Implementation
!  28 Jan 2021: Yeosang Yoon; Fix bug in calculate_sea_ice_concentration
!                             subroutine
!  08 Feb 2021: Yeosang Yoon; Remove unused variable 
!
! DESCRIPTION:
! Source code for the retrieval of snow depth and sea ice concentration from 
! XCAL GMI brightness temperature observations
!-------------------------------------------------------------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module USAFSI_xcalgmiMod

   ! Defaults
   implicit none
   private
   
   ! Public routines
   public :: USAFSI_proc_xcalgmi

contains

   ! Public routine for processing XCAL GMI data
   subroutine USAFSI_proc_xcalgmi(date10, gmi_in, gmi, option)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
      use USAFSI_utilMod
      
      ! Defaults
      implicit none

      ! Arguments
      character(len=10),  intent(in)  :: date10
      character(len=255), intent(in)  :: gmi_in
      character(len=255), intent(in)  :: gmi
      integer,            intent(in)  :: option

      ! Local variables
      integer                         :: eof, i, j, n, nArr, nFile, x, y
      character(len=255)              :: filename, nc_filename
      !integer,           dimension(:), allocatable :: surflag0, surflag   !YY
      character(len=10), dimension(:), allocatable :: date0, date10_arr
      real,              dimension(:), allocatable :: lat0, lon0, tb10h0, tb10v0, tb19h0, tb19v0, &
                                                      tb23v0, tb37h0, tb37v0, tb89h0, tb89v0
      real,              dimension(:), allocatable :: tb10h, tb10v, tb19h, tb19v, &
                                                      tb23v, tb37h, tb37v, tb89h, tb89v, &
                                                      snowdepth, ct
      real,              dimension(:), allocatable :: lat, lon
      integer,           dimension(:), allocatable :: qcflag, qcflag0
      integer                                      :: hemi
      integer                                      :: satid
      real                                         :: icoord, jcoord
      INTEGER, PARAMETER                           :: POSE    = 0   ! LONGITUDE ORIENTATION FLAG; 0 = POSITIVE WEST
 
      call search_files (date10, gmi_in)

      open(unit=90, file='./gmi_filelist.txt', form='formatted', &
           status='old', action='read')
 
      nFile = 0
      do
         read(90,'(A)',iostat=eof) filename
         if (eof/=0) exit
         
         write (LDT_logunit,*) '[INFO] get gmi obs:', trim(filename)
         
         call read_xcalgmi_attributes(filename, date0, lat0, lon0, tb10h0, tb10v0, &
                  tb19h0, tb19v0, tb23v0, tb37h0, tb37v0, tb89h0, tb89v0, qcflag0)

         ! resize the array
         if(nFile==0) then
            call move_alloc(date0, date10_arr)
            call move_alloc(lat0, lat)
            call move_alloc(lon0, lon)
            call move_alloc(tb10h0, tb10h)
            call move_alloc(tb10v0, tb10v)
            call move_alloc(tb19h0, tb19h)
            call move_alloc(tb19v0, tb19v)
            call move_alloc(tb23v0, tb23v)
            call move_alloc(tb37h0, tb37h)
            call move_alloc(tb37v0, tb37v)
            call move_alloc(tb89h0, tb89h)
            call move_alloc(tb89v0, tb89v)
            call move_alloc(qcflag0, qcflag)
         else
            call resize_array_str(date10_arr, date0)
            call resize_array(lat, lat0)
            call resize_array(lon, lon0)
            call resize_array(tb10h, tb10h0)
            call resize_array(tb10v, tb10v0)
            call resize_array(tb19h, tb19h0)
            call resize_array(tb19v, tb19v0)
            call resize_array(tb23v, tb23v0)
            call resize_array(tb37h, tb37h0)
            call resize_array(tb37v, tb37v0)
            call resize_array(tb89h, tb89h0)
            call resize_array(tb89v, tb89v0)
            call resize_array_int(qcflag, qcflag0)
         endif

         nFile = nFile + 1         

      end do
      close(90)

      ! calculate snow depth
      nArr=size(tb19h)
      if(allocated(snowdepth)) deallocate(snowdepth)
      if(allocated(ct))        deallocate(ct)
      allocate(snowdepth(nArr))
      allocate(ct(nArr))

      write (LDT_logunit,*) '[INFO] calculate snow depth from gmi'
      !call calculate_snowdepth(nArr, surflag, tb10h, tb10v, tb19h, tb19v, &
      !     tb23v, tb37h, tb37v, tb89h, tb89v, snowdepth, option, lat, lon, &
      !     qcflag)
      ! Yeosang Yoon
      call calculate_snowdepth(nArr, tb10h, tb10v, tb19h, tb19v, &
           tb23v, tb37h, tb37v, tb89h, tb89v, snowdepth, option, lat, lon, &
           qcflag)

      write (LDT_logunit,*) '[INFO] calculate sea ice concentration from gmi'
      !call calculate_sea_ice_concentration(nArr, surflag, lat, lon, tb19h, &
      !     tb19v, tb37h, tb37v, ct, qcflag)
      ! Yeosang Yoon
      call calculate_sea_ice_concentration(nArr, lat, lon, tb19h, &
           tb19v, tb37h, tb37v, ct, qcflag)

      ! write file (kept the name of "ssmis" to minimize code changes in other parts)
      filename = trim(gmi)//'ssmis_snoice_nh.06hr.'//date10//'.txt'
      write (LDT_logunit,*) &
           '[INFO] Writing GMI data to ', trim(filename)
      open(unit=10, file=filename,status='unknown', action='write')
      filename = trim(gmi)//'ssmis_snoice_sh.06hr.'//date10//'.txt'
      write (LDT_logunit,*) &
           '[INFO] Writing GMI data to ', trim(filename)
      open(unit=20, file=filename,status='unknown', action='write')

      satid = 13  !satid does not work for GMI, but need a space for satid in txt file 
                  !for further procecure in USAFSI, which was origianlly developed for SSMIS
      do i=1,size(lat)
         if (snowdepth(i) < 0) then      ! remove unrealistic values
            snowdepth(i)=-0.1
         end if

         ! lat, lon*100, snow depth(mm)*10
         !if (surflag(i)==0 .and. snowdepth(i)>=0) then      !  land points only
         if (ct(i) >=0 .or. snowdepth(i) >=0) then
            if (lat(i) > 0) then
               hemi=1
               call LLTOPS (POSE, lat(i), lon(i), 16, hemi, icoord, jcoord)
               write(10, '(A10, I3, I6, I7, 2I5, 2I6)') date10_arr(i), &
                    satid, nint(lat(i)*100), nint(lon(i)*100), &
                    nint(icoord), nint(jcoord), nint(ct(i)), &
                    nint(snowdepth(i)*10)
            else
               hemi=2
               call LLTOPS (POSE, lat(i), lon(i), 16, hemi, icoord, jcoord)
               write(20, '(A10, I3, I6, I7, 2I5, 2I6)') date10_arr(i), &
                    satid, nint(lat(i)*100), nint(lon(i)*100), &
                    nint(icoord), nint(jcoord), nint(ct(i)), &
                    nint(snowdepth(i)*10)
            end if ! if (lat(i) > 0) then
         end if
      end do ! do i=1,size(lat)
      close(10)
      close(20)

      ! write netCDF file for USAFSI
      nc_filename = trim(gmi)//'gmi_snoice_0p25deg.'//date10//'.nc'
      write (LDT_logunit,*) &
           '[INFO] Writing GMI data to ', trim(nc_filename)
      call write_netcdf(nc_filename, nArr, lat, lon, snowdepth, ct)

      if(allocated(lat0)) deallocate(lat0)
      if(allocated(lon0)) deallocate(lon0)
      if(allocated(tb10h0)) deallocate(tb10h0)
      if(allocated(tb10v0)) deallocate(tb10v0)
      if(allocated(tb19h0)) deallocate(tb19h0)
      if(allocated(tb19v0)) deallocate(tb19v0)
      if(allocated(tb23v0)) deallocate(tb23v0)
      if(allocated(tb37h0)) deallocate(tb37h0)
      if(allocated(tb37v0)) deallocate(tb37v0)
      if(allocated(tb89h0)) deallocate(tb89h0)
      if(allocated(tb89v0)) deallocate(tb89v0)

      if(allocated(lat)) deallocate(lat)
      if(allocated(lon)) deallocate(lon)
      if(allocated(tb10h)) deallocate(tb10h)
      if(allocated(tb10v)) deallocate(tb10v)
      if(allocated(tb19h)) deallocate(tb19h)
      if(allocated(tb19v)) deallocate(tb19v)
      if(allocated(tb23v)) deallocate(tb23v)
      if(allocated(tb37h)) deallocate(tb37h)
      if(allocated(tb37v)) deallocate(tb37v)
      if(allocated(tb89h)) deallocate(tb89h)
      if(allocated(tb89v)) deallocate(tb89v)
      
   end subroutine USAFSI_proc_xcalgmi

   ! *** Remaining routines are private ***

   subroutine read_xcalgmi_attributes(filename, date10_arr, lat, lon, tb10h, tb10v, &
                   tb19h, tb19v, tb23v, tb37h, tb37v, tb89h, tb89v, qcflag_1d)

#if (defined USE_HDF5)
      use hdf5
#endif

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify

      !Defaults
      implicit none

      !Arguments
      character(len=255), intent(in) :: filename
      character(len=10), allocatable, intent(inout) :: date10_arr(:)
      real,         allocatable, intent(inout) :: lat(:), lon(:)
      real,         allocatable                :: lat_raw(:,:), lon_raw(:,:)
      real,         allocatable, intent(inout) :: tb10h(:), tb10v(:), tb19h(:), tb19v(:), &
                                                  tb23v(:), tb37h(:), tb37v(:), tb89h(:), &
                                                  tb89v(:)
      integer,      allocatable, intent(inout) :: qcflag_1d(:)
      real,         allocatable                :: tb_field(:,:,:)
#if (defined USE_HDF5)
      integer(hsize_t), dimension(2)           :: dims_lat, maxdims_lat
      integer(hsize_t), dimension(2)           :: dims_lon, maxdims_lon
      integer(hsize_t), dimension(3)           :: dims_tb, maxdims_tb
      integer(hsize_t), dimension(2)           :: dims_qc, maxdims_qc
      integer(hsize_t), dimension(1)           :: dims_yr, dims_mon, dims_day, dims_hr
      integer(hsize_t), dimension(2)           :: maxdims_yr, maxdims_mon, maxdims_day, maxdims_hr
      integer(hid_t)                           :: file_id, gr_id, grst_id
      integer(hid_t)                           :: dspace_lat_id, dspace_lon_id
      integer(hid_t)                           :: dspace_tb_id
      integer(hid_t)                           :: dspace_yr_id, dspace_mon_id
      integer(hid_t)                           :: dspace_day_id, dspace_hr_id 
      integer(hid_t)                           :: tb_field_id, lat_id, lon_id
      integer(hid_t)                           :: yr_id, mon_id, day_id, hr_id
      integer(hid_t)                           :: qc_id, dspace_qc_id
      integer                                  :: status
      integer                                  :: dim_1d
      integer                                  :: sta_arr, end_arr
      integer                                  :: idims01, inscan
      integer(2), allocatable                  :: yyyy(:)
      integer(1), allocatable                  :: mm(:), dd(:), hh(:)
      integer(1), allocatable                  :: qcflag(:,:)
      character(len=10)                        :: date_st
      character(len=10), allocatable           :: date_st_arr(:)

      !---------
      ! #filespec of the XCAL GMI
      ! - (nchan1) number of channels in Swath 1 (S1) = 9
      ! - (nfreq1) number of frequencies in S1 = 5
      ! - (npix1) number of pixels in S1 = 221
      ! - Swath S1 has channels 1-9: 10V 10H 19V 19H 23V 37V 37H 89V 89H
      !---------

      call h5open_f(status)
      call LDT_verify(status, 'Error opening HDF fortran interface')

      call h5fopen_f(trim(filename),H5F_ACC_RDONLY_F, file_id, status)
      call LDT_verify(status, 'Error opening XCAL GMI file ')

      call h5gopen_f(file_id,"S1",gr_id, status)
      call LDT_verify(status, 'Error opening S1 group in XCAL GMI file')

      call h5gopen_f(file_id,"/S1/ScanTime",grst_id, status)
      call LDT_verify(status, 'Error opening S1/ScanTime group in XCAL GMI file')

      call h5dopen_f(gr_id,"Tc",tb_field_id, status)   !Tc dimension: nchan1 x npix1 x nscan
      call LDT_verify(status, 'Error opening Tb field in XCAL GMI file')

      call h5dopen_f(gr_id,"Latitude",lat_id, status)  !dimension: npix1 x nscan
      call LDT_verify(status, 'Error opening Latitude field in XCAL GMI file')

      call h5dopen_f(gr_id,"Longitude",lon_id, status) !dimension: npix1 x nscan
      call LDT_verify(status, 'Error opening Longitude field in XCAL GMI file')

      call h5dopen_f(grst_id,"Year",yr_id, status) !dimension: nscan
      call LDT_verify(status, 'Error opening Year field in XCAL GMI file')

      call h5dopen_f(grst_id,"Month",mon_id, status) !dimension: nscan
      call LDT_verify(status, 'Error opening Month field in XCAL GMI file')

      call h5dopen_f(grst_id,"DayOfMonth",day_id, status) !dimension: nscan
      call LDT_verify(status, 'Error opening DayOfMonth field in XCAL GMI file')

      call h5dopen_f(grst_id,"Hour",hr_id, status) !dimension: nscan
      call LDT_verify(status, 'Error opening Hour field in XCAL GMI file')

      call h5dopen_f(gr_id,"Quality",qc_id, status) !dimension: npix1 x nscan
      call LDT_verify(status, 'Error opening Quality field in XCAL GMI file')

      call h5dget_space_f(lat_id, dspace_lat_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI lat')

      call h5sget_simple_extent_dims_f(dspace_lat_id, dims_lat, maxdims_lat, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: Latitude')
      endif

      call h5dget_space_f(lon_id, dspace_lon_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI lon')

      call h5sget_simple_extent_dims_f(dspace_lon_id, dims_lon, maxdims_lon, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: Longitude')
      endif

      call h5dget_space_f(yr_id, dspace_yr_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI Year')

      call h5sget_simple_extent_dims_f(dspace_yr_id, dims_yr, maxdims_yr, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: Year')
      endif

      call h5dget_space_f(mon_id, dspace_mon_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI Month')

      call h5sget_simple_extent_dims_f(dspace_mon_id, dims_mon, maxdims_mon, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: Month')
      endif

      call h5dget_space_f(day_id, dspace_day_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI Day')

      call h5sget_simple_extent_dims_f(dspace_day_id, dims_day, maxdims_day, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: DayOfMonth')
      endif

      call h5dget_space_f(hr_id, dspace_hr_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI Hour')

      call h5sget_simple_extent_dims_f(dspace_hr_id, dims_hr, maxdims_hr, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: Hour')
      endif

      call h5dget_space_f(tb_field_id, dspace_tb_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI tb')

      call h5sget_simple_extent_dims_f(dspace_tb_id, dims_tb, maxdims_tb, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: tb')
      endif

      call h5dget_space_f(qc_id, dspace_qc_id, status)
      call LDT_verify(status, 'Error in h5dget_space_f: read GMI Quality')

      call h5sget_simple_extent_dims_f(dspace_qc_id, dims_qc, maxdims_qc, status)
      if(status.eq.-1) then
         call LDT_verify(status, 'Error in h5sget_simple_extent_dims_f: Quality')
      endif

      if(allocated(lat_raw)) deallocate(lat_raw)
      if(allocated(lon_raw)) deallocate(lon_raw)
      if(allocated(tb_field)) deallocate(tb_field)
      if(allocated(qcflag)) deallocate(qcflag)
      if(allocated(yyyy)) deallocate(yyyy)
      if(allocated(mm)) deallocate(mm)
      if(allocated(dd)) deallocate(dd)
      if(allocated(hh)) deallocate(hh)
      allocate(lat_raw(dims_lat(1),dims_lat(2)))
      allocate(lon_raw(dims_lon(1),dims_lon(2)))
      allocate(tb_field(dims_tb(1),dims_tb(2),dims_tb(3)))
      allocate(qcflag(dims_qc(1),dims_qc(2)))
      allocate(yyyy(dims_yr(1)))
      allocate(mm(dims_mon(1)))
      allocate(dd(dims_day(1)))
      allocate(hh(dims_hr(1)))

      call h5dread_f(lat_id, H5T_IEEE_F32LE, lat_raw, dims_lat, status)
      call LDT_verify(status, 'Error extracting Latitude field from XCAL GMI file')

      call h5dread_f(lon_id, H5T_IEEE_F32LE, lon_raw, dims_lon, status)
      call LDT_verify(status, 'Error extracting Longitude field from XCAL GMI file')

      call h5dread_f(tb_field_id, H5T_IEEE_F32LE, tb_field, dims_tb, status)
      call LDT_verify(status, 'Error extracting TB field from XCAL GMI file')

      call h5dread_f(qc_id, H5T_STD_I8LE, qcflag, dims_qc, status)
      call LDT_verify(status, 'Error extracting Quality field from XCAL GMI file')

      call h5dread_f(yr_id, H5T_STD_I16LE, yyyy, dims_yr, status)
      call LDT_verify(status, 'Error extracting Year field from XCAL GMI file')

      call h5dread_f(mon_id, H5T_STD_I8LE, mm, dims_mon, status)
      call LDT_verify(status, 'Error extracting Month field from XCAL GMI file')

      call h5dread_f(day_id, H5T_STD_I8LE, dd, dims_day, status)
      call LDT_verify(status, 'Error extracting DayOfMonth field from XCAL GMI file')

      call h5dread_f(hr_id, H5T_STD_I8LE, hh, dims_hr, status)
      call LDT_verify(status, 'Error extracting Hour field from XCAL GMI file')

      if(allocated(date_st_arr)) deallocate(date_st_arr)
      if(allocated(date10_arr)) deallocate(date10_arr)
      if(allocated(lat)) deallocate(lat)
      if(allocated(lon)) deallocate(lon)
      if(allocated(tb10h)) deallocate(tb10h)
      if(allocated(tb10v)) deallocate(tb10v)
      if(allocated(tb19h)) deallocate(tb19h)
      if(allocated(tb19v)) deallocate(tb19v)
      if(allocated(tb23v)) deallocate(tb23v)
      if(allocated(tb37h)) deallocate(tb37h)
      if(allocated(tb37v)) deallocate(tb37v)
      if(allocated(tb89h)) deallocate(tb89h)
      if(allocated(tb89v)) deallocate(tb89v)
      if(allocated(qcflag_1d)) deallocate(qcflag_1d)

      dim_1d = dims_lat(1) * dims_lat(2)  ! npix (221) * nscan (e.g., 160)
      
      allocate(date_st_arr(dims_lat(2)))
      allocate(date10_arr(dim_1d))
      allocate(lat(dim_1d))
      allocate(lon(dim_1d))
      allocate(tb10h(dim_1d))
      allocate(tb10v(dim_1d))
      allocate(tb19h(dim_1d))
      allocate(tb19v(dim_1d))
      allocate(tb23v(dim_1d))
      allocate(tb37h(dim_1d))
      allocate(tb37v(dim_1d))
      allocate(tb89h(dim_1d))
      allocate(tb89v(dim_1d))
      allocate(qcflag_1d(dim_1d))

      do inscan = 1, dims_lat(2)
         write(date_st(1:4),'(I4)') yyyy(inscan)
         write(date_st(5:6),'(I0.2)') mm(inscan)
         write(date_st(7:8),'(I0.2)') dd(inscan)
         write(date_st(9:10),'(I0.2)') hh(inscan)

         date_st_arr(inscan) = date_st
      end do
         
      sta_arr = 1
      end_arr = dims_lat(2)   !nscan (e.g., 160)
      do idims01 = 1, dims_lat(1)  !221
         date10_arr(sta_arr:end_arr) = date_st_arr(:)
         lat(sta_arr:end_arr) = lat_raw(idims01,:)
         lon(sta_arr:end_arr) = lon_raw(idims01,:)
         tb10v(sta_arr:end_arr) = tb_field(1,idims01,:)
         tb10h(sta_arr:end_arr) = tb_field(2,idims01,:) 
         tb19v(sta_arr:end_arr) = tb_field(3,idims01,:)
         tb19h(sta_arr:end_arr) = tb_field(4,idims01,:)
         tb23v(sta_arr:end_arr) = tb_field(5,idims01,:)
         tb37v(sta_arr:end_arr) = tb_field(6,idims01,:)
         tb37h(sta_arr:end_arr) = tb_field(7,idims01,:)
         tb89v(sta_arr:end_arr) = tb_field(8,idims01,:)
         tb89h(sta_arr:end_arr) = tb_field(9,idims01,:)
         qcflag_1d(sta_arr:end_arr) = qcflag(idims01,:)

         sta_arr = sta_arr + dims_lat(2)
         end_arr = end_arr + dims_lat(2)
      enddo

      call h5dclose_f(lat_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5dclose_f(lon_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5dclose_f(tb_field_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5dclose_f(qc_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5dclose_f(yr_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5dclose_f(mon_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5dclose_f(day_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5dclose_f(hr_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5gclose_f(gr_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5gclose_f(grst_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5fclose_f(file_id,status)
      call LDT_verify(status,'Error in H5DCLOSE call')

      call h5close_f(status)
      call LDT_verify(status,'Error in H5CLOSE call') 
#endif
   end subroutine read_xcalgmi_attributes


   !subroutine calculate_snowdepth(n, surflag, tb10h, tb10v, tb19h, tb19v, &
   !               tb23v, tb37h, tb37v, tb89h, tb89v, snowdepth, option, lat, lon, &
   !               qcflag) 
   subroutine calculate_snowdepth(n, tb10h, tb10v, tb19h, tb19v, &
                  tb23v, tb37h, tb37v, tb89h, tb89v, snowdepth, option, lat, lon, &
                  qcflag)
      ! option == 1: Hollinger, 1991,    SD=4445.0-17.95TB_37V
      ! option == 2: Markus (Chang et al., 1987), SD=1.58(TB_19H-TB_37H)
      ! option == 3; Foster et al., 1997
      ! option == 4; Kelly, 2009 

      ! Imports
      use LDT_usafsiMod, only: usafsi_settings
 
      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in)   :: n
      !integer, intent(in)   :: surflag(:)   !YY
      real,    intent(in)   :: tb10h(:)
      real,    intent(in)   :: tb10v(:)
      real,    intent(in)   :: tb19h(:)
      real,    intent(in)   :: tb19v(:)
      real,    intent(in)   :: tb23v(:)
      real,    intent(in)   :: tb37h(:)
      real,    intent(in)   :: tb37v(:)
      real,    intent(in)   :: tb89h(:)
      real,    intent(in)   :: tb89v(:)
      real,    intent(out)  :: snowdepth(1:n)
      integer, intent(in)   :: option
      real,    intent(in)   :: lat(:)
      real,    intent(in)   :: lon(:)
      integer, intent(in)   :: qcflag(:)

      integer, parameter                  :: nc=1440
      integer, parameter                  :: nr=720

      ! Local variables
      integer                             :: i
      character(len=255)                  :: ff_filename
      character(len=255)                  :: fd_filename
      real                                :: pd19,pd37,pd89,tt,si89, &
           scat,sc37,sc89,scx
      logical                             :: flag
      real                                :: lon_grid(nc), &
           lat_grid(nr), ratio
      real                                :: ff(nc,nr)
      real                                :: fd(nc,nr)
      integer                             :: plat, plon

      real                                :: Tphy
      real                                :: SD_f, SD_o

      ! set 0.25 deg grid, for forest fraction
      if (option==3 .or. option==4) then
         do i=1, nc
            lon_grid(i)=-179.875+0.25*(i-1)
         end do
         do i=1, nr
            lat_grid(i)=-89.875+0.25*(i-1)
         end do

         ff_filename=trim(usafsi_settings%ff_file)
         ! read forest fraction    
         call read_forestfraction(ff_filename,ff)

         if (option==4) then
            fd_filename=trim(usafsi_settings%fd_file)
            ! read forest density
            call read_forestdensity(fd_filename,fd)
         endif
      endif

      do i=1,n
         ! quality check
         pd19=tb19v(i)-tb19h(i)
         pd37=tb37v(i)-tb37h(i)
         pd89=tb89v(i)-tb89h(i)
         tt=169+0.5*tb89v(i)
         sc89=tb23v(i)-tb89v(i)-3.0
         sc37=tb19v(i)-tb37v(i)-3.0
         scx=tb37v(i)-tb89v(i)-1.0
         scat=sc89
         if (sc37>scat) then
            scat=sc37
         end if
         if (scx>scat) then
            scat=scx
         end if

         ! NG, 2002
         if (scat>0) then ! scattering materials
            flag=.true.
            if (tb23v(i)>=260.0 .and. pd89>3.0 .and. scat<=3.0) then  ! cold rain
               flag=.false.
            end if
            if (tb23v(i)>=264.0 .or. tb23v(i)>=tt) then  ! other rain events
               flag=.false.
            end if
            if (pd19>=18.0 .and. sc37<=10.0 .and. scx<=10.0) then ! cold desert
               flag=.false.
            end if
            if (pd19>=8.0 .and. sc89<=6.0 .and. sc37<=2.0) then ! frozen ground
               flag=.false.
            end if
            if ((tb23v(i)<=216.0) .or. ((tb23v(i)<=235.0) .and. (pd19>=23.0))) then ! gracier
               flag=.false.
            end if
            if (tb19v(i)<-100 .or. tb19h(i)<-100 .or. tb37v(i)<-100 .or. tb37h(i)<-100) then  !check NaN value (-9999)
               flag=.false.
            end if
         else
            flag=.false.
         end if
         
         ! qcflag
         if (qcflag(i) > 0) then
            flag=.false.
         end if

         !Get SurfaceFlag (0:Land, 1:Reserved, 2:Near coast, 3:Ice, 4:Possilbe ice, 5:Ocean, 6:Coast,
         !                  7-14: Reserved, 15: Missing value)
         if (flag .eqv. .true.) then
         !if ((flag .eqv. .true.) .and. (surflag(i)==0 .or. surflag(i)==2)) then
            if (option==1) then
               snowdepth(i)=4445.0-17.95*tb37v(i)
            else if (option==2) then
               snowdepth(i)=1.58*(tb19h(i)-tb37h(i))*10.0  ! unit cm -> mm
               !if (snowdepth(i) <= 25) then   ! derived snow depth less than 2.5 cm will be assigned as no snow.
               !  snowdepth(i) = -0.1
               !end if
            else if (option==3) then
               call checkgrid(lat_grid,lon_grid,lat(i),lon(i),plat,plon)

               ! check forest fraction
               if (ff(plon,plat)<0) then
                  ff(plon,plat)=0
               else if (ff(plon,plat)>0.9) then
                  ff(plon,plat)=0.9
               end if
               snowdepth(i)=0.78*(1/(1-ff(plon,plat)))*(tb19h(i)-tb37h(i))*10.0  ! unit cm -> mm
               !if (snowdepth(i) <= 25) then   ! derived snow depth less than 2.5 cm will be assigned as no snow.
               !  snowdepth(i) = -0.1
               !end if
               !print *, snowdepth(i)
            else !option == 4 (Kelly, 2009)
               call checkgrid(lat_grid,lon_grid,lat(i),lon(i),plat,plon)
 
               ! check forest fraction
               if (ff(plon,plat)<0) then
                  ff(plon,plat)=0
               else if (ff(plon,plat)>0.9) then
                  ff(plon,plat)=0.9
               end if

               ! check forest density
               fd(plon,plat) = fd(plon,plat)/100 !% --> fraction
               if (ff(plon,plat)==0) then
                  fd(plon,plat) = 0
               elseif (ff(plon,plat)>0) then
                  if (fd(plon,plat)<0) then
                     fd(plon,plat) = 0.001
                  elseif (fd(plon,plat)>1) then
                     fd(plon,plat) = 1
                  endif
               end if

               ! Calculate land surface physical temperature 
               Tphy = 58.08 - 0.39*tb19v(i) + 1.21*tb23v(i) - 0.37*tb37h(i) + 0.36*tb89v(i)  !K

               ! Test for moderate to deep snow presence
               ! Threshold are checked to ensure cold snow conditions are potentially present 
               !    in the 36 GHz Tb

               if (tb37h(i) < 245 .and. tb37v(i) < 255) then
                  ! snow is possible
                  ! shallow or medium depth of snow is retrieved
         
                  ! A minimum polarization difference of 3 is applied to ensure that
                  ! the denominator is not too small (Kelly, 2009)
                  if (pd19 < 3) then
                     pd19 = 3
                  end if
                  if (pd37 < 3) then
                     pd37 = 3
                  end if

                     !forest component
                  SD_f = (1/log10(pd37)) * (tb19v(i)-tb37v(i))/(1-fd(plon,plat)*0.6)     !cm
                     
                  if (SD_f < 0) then
                     SD_f = 0
                  end if                     

                     !non-forested component
                  SD_o = ((1/log10(pd37)) * (tb10v(i)-tb37v(i))) + &
                           ((1/log10(pd19)) * (tb10v(i)-tb19v(i)))      !cm

                  if (SD_o < 0) then
                     SD_o = 0
                  end if

                  snowdepth(i) = (ff(plon,plat)*SD_f) + ((1-ff(plon,plat))*SD_o)    !cm
                  snowdepth(i) = snowdepth(i) * 10.0    !cm -> mm
               else
                  ! test for shallow snow
                  if (tb10v(i)-tb37v(i)>0 .or. tb10h(i)-tb37h(i)>0) then
                     ! medium to deep snow is assumed to be present

                     ! A minimum polarization difference of 3 is applied to ensure that
                     ! the denominator is not too small (Kelly, 2009)
                     if (pd19 < 3) then
                        pd19 = 3
                     end if
                     if (pd37 < 3) then
                        pd37 = 3
                     end if

                        !forest component
                     SD_f = (1/log10(pd37)) * (tb19v(i)-tb37v(i))/(1-fd(plon,plat)*0.6)     !cm
                       
                     if (SD_f < 0) then
                        SD_f = 0
                     end if

                        !non-forested component
                     SD_o = ((1/log10(pd37)) * (tb10v(i)-tb37v(i))) + &
                              ((1/log10(pd19)) * (tb10v(i)-tb19v(i)))      !cm

                     if (SD_o < 0) then
                        SD_o = 0
                     end if

                     snowdepth(i) = (ff(plon,plat)*SD_f) + ((1-ff(plon,plat))*SD_o)    !cm
                     snowdepth(i) = snowdepth(i) * 10.0    !cm -> mm
                  else
                     !
                     ! snow presence is possible, but it is likely to be shallow snow
                     if (tb89v(i)<255 .and. tb89h(i)<265 .and. tb23v(i)-tb89v(i)>0 .and. Tphy<267) then
                        ! "tb23h-tb89h>0" has not been used for GMI due to the absence of tb23h 
                        ! snow depth is assumed to be 5.0 cm
                        snowdepth(i) = 50 ! mm
                     else
                        snowdepth(i) = -0.1
                     end if
                  end if                  
               end if

            end if ! if (option==1) then
         else
            snowdepth(i) = -0.1
         end if
      end do

   end subroutine calculate_snowdepth


   !subroutine calculate_sea_ice_concentration(n, surflag, lat, lon, tb19h, &
   !     tb19v, tb37h, tb37v, ct, qcflag)
   subroutine calculate_sea_ice_concentration(n, lat, lon, tb19h, &
        tb19v, tb37h, tb37v, ct, qcflag)

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in)   :: n
      !integer, intent(in)   :: surflag(:) !YY
      real,    intent(in)   :: lat(:)
      real,    intent(in)   :: lon(:)
      real,    intent(in)   :: tb19h(:)
      real,    intent(in)   :: tb19v(:)
      real,    intent(in)   :: tb37h(:)
      real,    intent(in)   :: tb37v(:)
      real,    intent(out)  :: ct(:)
      integer, intent(in)   :: qcflag(:)

      ! Local variables
      integer                                      :: i
      real                                         :: pr, gr, d, cf, cm
      logical                                      :: flag

      !flag=.true.

      !TODO: need mask only for ocean tile
      do i=1,n
         flag=.true.  !Yeosang Yoon

         ! qcflag
         if (qcflag(i) > 0) then
            flag=.false.
         end if

         if (tb19v(i)<-100 .or. tb19h(i)<-100 .or. tb37v(i)<-100 .or. tb37h(i)<-100) then  !check NaN value (-9999)
            flag=.false.
         end if

         if (lat(i)>44.5 .or. lat(i)<-52) then
            if (flag .eqv. .true.) then
               pr=(tb19v(i)-tb19h(i))/(tb19v(i)+tb19h(i))   ! polarization ratio
               gr=(tb37v(i)-tb19v(i))/(tb37v(i)+tb19v(i))   ! spectral gradient ratio
               d=2035.3+9244.6*pr-5665.8*gr-12875.1*pr*gr
               cf=(3290.2-20761.2*pr+23934.0*gr+47985.4*pr*gr)/d   ! first-year ice concentration
               cm=(-790.9+13825.3*pr-33155.8*gr-47771.9*pr*gr)/d   ! multiyear ice concentration
               ct(i)=(cf+cm)*100                                   ! sea ice concentration [0 100]
               if (ct(i) > 100) then
                  ct(i)=100
               end if
            else
               ct(i)=-1
            end if !flag
         else
            ct(i)=-1
         end if ! lat
      end do
      
   end subroutine calculate_sea_ice_concentration


   subroutine resize_array(arr1, arr2)

      ! Defaults
      implicit none

      ! Arguments
      real, allocatable, intent(inout) :: arr1(:)
      real, allocatable, intent(in)    :: arr2(:)

      ! Local variables
      real, dimension(:), allocatable  :: temp
      integer :: nArr

      call move_alloc(arr1, temp)

      nArr=size(temp)+size(arr2)
      allocate(arr1(nArr))

      arr1(1:size(temp))=temp
      arr1(size(temp)+1:nArr)=arr2

      deallocate(temp)

   end subroutine resize_array


   subroutine resize_array_int(arr1, arr2)

      ! Defaults
      implicit none

      ! Arguments
      integer, allocatable, intent(inout) :: arr1(:)
      integer, allocatable, intent(in)    :: arr2(:)

      ! Local variables
      integer, dimension(:), allocatable                :: temp
      integer :: nArr

      call move_alloc(arr1, temp)

      nArr=size(temp)+size(arr2)
      allocate(arr1(nArr))

      arr1(1:size(temp))=temp
      arr1(size(temp)+1:nArr)=arr2

      deallocate(temp)

   end subroutine resize_array_int

   
   subroutine resize_array_str(arr1, arr2)

      ! Defaults
      implicit none

      ! Arguments
      character(len=10), allocatable, intent(inout) :: arr1(:)
      character(len=10), allocatable, intent(in)    :: arr2(:)

      ! Local variables
      character(len=10), dimension(:), allocatable                :: temp
      integer :: nArr

      call move_alloc(arr1, temp)

      nArr=size(temp)+size(arr2)
      allocate(arr1(nArr))

      arr1(1:size(temp))=temp
      arr1(size(temp)+1:nArr)=arr2

      deallocate(temp)

   end subroutine resize_array_str


#if(defined USE_NETCDF4) 
   subroutine write_netcdf(nc_filename, n, lat, lon, snowdepth, ct)
   
      ! Imports
      use LDT_logMod, only: LDT_verify
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in)       :: nc_filename
      integer, intent(in) :: n
      real,    intent(in) :: lat(1:n)
      real,    intent(in) :: lon(1:n)
      real,    intent(in) :: snowdepth(1:n)
      real,    intent(in) :: ct(1:n)

      ! Local constants
      integer, parameter                             :: nc=1440 ! fixed 1/4 deg
      integer, parameter                             :: nr=720
      real, parameter                                :: fillval = -1.0

      ! Local variables
      real, dimension(1:nc)                          :: lon_grid
      real, dimension(1:nr)                          :: lat_grid
      real, dimension(1:nc,1:nr)                     :: gmidep, ndep, gmicon, &
           ncon
      integer                                        :: plat, plon
      integer                                        :: iret, ncid, &
           gmidep_varid, gmicon_varid, lat_varid, lon_varid, dim_ids(2)
      integer :: i,x,y

      ! set 0.25 deg grid
      do i=1, nc
         lon_grid(i)=-179.875+0.25*(i-1)
      end do
      do i=1, nr
         lat_grid(i)=-89.875+0.25*(i-1)
      end do

      ! initialization
      gmidep(1:nc,1:nr)=0.0
      ndep(1:nc,1:nr)=0.0
      gmicon(1:nc,1:nr)=0.0
      ncon(1:nc,1:nr)=0.0

      do i=1,n

         ! find nearest grid location
         call checkgrid(lat_grid,lon_grid,lat(i),lon(i),plat,plon)

         ! for snow depth
         if (snowdepth(i)>=0) then
            gmidep(plon,plat)=gmidep(plon,plat)+snowdepth(i)
            ndep(plon,plat)=ndep(plon,plat)+1
         end if

         ! for ice concentration
         if (ct(i)>=0) then
            gmicon(plon,plat)=gmicon(plon,plat)+ct(i)
            ncon(plon,plat)=ncon(plon,plat)+1
         end if

      end do ! i=1,size(lat)

      ! Averaging snow depth & ice concentration
      do x=1,nc
         do y=1,nr
            if (ndep(x,y)>=1) then
               gmidep(x,y)=gmidep(x,y)/ndep(x,y)
            else
               gmidep(x,y)=fillval   ! apply null vaule
            end if

            if (ncon(x,y)>=1) then
               gmicon(x,y)=gmicon(x,y)/ncon(x,y)
            else
               gmicon(x,y)=fillval   ! apply null vaule
            end if
         end do ! do y
      end do ! do x

      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
      ! overwrite this file, if it already exists.
      iret=nf90_create(path=trim(nc_filename),cmode=nf90_netcdf4, ncid=ncid)
      call LDT_verify(iret, '[ERR] nf90_create failed')

      ! Define the dimensions. NetCDF will hand back an ID for each.
      call LDT_verify(nf90_def_dim(ncid,'lat',nr,dim_ids(2)), &
           '[ERR] nf90_def_dim failed')
      call LDT_verify(nf90_def_dim(ncid,'lon',nc,dim_ids(1)), &
           '[ERR] nf90_def_dim failed')

      ! Add map projection
      ! Hard code: 0.25 deg and global lat/lon
      call LDT_verify(nf90_put_att(ncid,nf90_global,&
           "MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL"), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global,&
           "SOUTH_WEST_CORNER_LAT", -89.875), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global,&
           "SOUTH_WEST_CORNER_LON", -179.875), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global, &
           "NORTH_EAST_CORNER_LAT", 89.875), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global, &
           "NORTH_EAST_CORNER_LON", 179.875), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global, &
           "DX",0.25),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global, &
           "DY",0.25), &
           '[ERR] nf90_put_att failed')

      ! Define the variable.
      ! latitudes
      call LDT_verify(nf90_def_var(ncid,"lat",nf90_float,dim_ids(2), &
           lat_varid),'[ERR] nf90_def_var failed')
      call LDT_verify(nf90_put_att(ncid,lat_varid, &
           "units","degrees_north"), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,lat_varid, &
           "long_name","latitude"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,lat_varid, &
           "standard_name","latitude"),&
           '[ERR] nf90_put_att failed')

      ! longitudes
      call LDT_verify(nf90_def_var(ncid,"lon",nf90_float,dim_ids(1), &
           lon_varid),'[ERR] nf90_def_var failed')
      call LDT_verify(nf90_put_att(ncid,lon_varid, &
           "units","degrees_east"), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,lon_varid, &
           "long_name","longitude"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,lon_varid, &
           "standard_name","longitude"),&
           '[ERR] nf90_put_att failed')

      ! gmi snow depth
      call LDT_verify(nf90_def_var(ncid,"gmidep",nf90_float, &
           dimids=dim_ids, varid=gmidep_varid),'[ERR] nf90_def_var failed')
      call LDT_verify(nf90_put_att(ncid,gmidep_varid, &
           "units","mm"),'[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,gmidep_varid, &
           "long_name","gmi snow depth of surface snow over land"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,gmidep_varid, &
           "standard_name","gmi_snow_depth_of_surface_snow"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,gmidep_varid, &
           '_FillValue',fillval), &
           '[ERR] nf90_put_att failed for SMIDEP')

      ! sea ice concentration analysis
      call LDT_verify(nf90_def_var(ncid,"gmicon",nf90_float, &
           dimids=dim_ids, varid=gmicon_varid),'[ERR] nf90_def_var failed')
      call LDT_verify(nf90_put_att(ncid,gmicon_varid, &
           "units","%"),'[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,gmicon_varid, &
           "long_name","concentration of sea ice (0-100)"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,gmicon_varid, &
           "standard_name","sea_ice_area_fraction"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,gmicon_varid, &
           'missing_value',fillval), &
           '[ERR] nf90_put_att failed for SMICON')
      call LDT_verify(nf90_put_att(ncid,gmicon_varid, &
           '_FillValue',fillval), &
           '[ERR] nf90_put_att failed for SMICON')
      call LDT_verify(nf90_put_att(ncid,gmicon_varid, &
           'valid_range',(/0.,100./)), &
           '[ERR] nf90_put_att failed for SMICON')

      ! Miscellaneous header information
      call LDT_verify(nf90_put_att(ncid,nf90_global,"title", &
           "LDT GMI snow depth analysis"), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global,"institution", &
           "NASA GSFC Hydrological Sciences Laboratory"), &
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,nf90_global,"source", &
           "Land Data Toolkit (LDT)"), &
           '[ERR] nf90_put_att failed')

      ! End define mode. This tells netCDF we are done defining metadata.
      call LDT_verify(nf90_enddef(ncid),'[ERR] ncf90_enddef failed')

      ! Write the lat/lon data
      call LDT_verify(nf90_put_var(ncid,lat_varid,lat_grid(:),&
           (/1/),(/nr/)),'[ERR] nf90_put_var failed for lat')
      call LDT_verify(nf90_put_var(ncid,lon_varid,lon_grid(:),&
           (/1/),(/nc/)), '[ERR] nf90_put_var failed for lon')

      ! Write the GMI snow depth/ice concentration fields
      call LDT_verify(nf90_put_var(ncid,gmidep_varid,&
           gmidep(:,:), (/1,1/),(/nc,nr/)), &
           '[ERR] nf90_put_var failed for gmidep')
      call LDT_verify(nf90_put_var(ncid,gmicon_varid,gmicon(:,:), &
           (/1,1/),(/nc,nr/)), '[ERR] nf90_put_var failed for icecon')

      ! Close the file. This frees up any internal netCDF resources
      ! associated with the file, and flushes any buffers.
      call LDT_verify(nf90_close(ncid),'[ERR] nf90_close failed!')

   end subroutine write_netcdf

#else

   ! Dummy version
   subroutine write_netcdf(nc_filename, n, lat,lon, snowdepth, ct)
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      implicit none
      character(len=*), intent(in)       :: nc_filename
      integer, intent(in) :: n
      real,    intent(in) :: lat(1:n)
      real,    intent(in) :: lon(1:n)
      real,    intent(in) :: snowdepth(1:n)
      real,    intent(in) :: ct(1:n)
      write(LDT_logunit,*) &
           "[ERR] Recompiled LDT with netCDF4 support!"
      write(LDT_logunit,*) &
           "[ERR] Stopping in write_netcdf..."
      call LDT_endrun()
   end subroutine write_netcdf
#endif


   subroutine checkgrid(lat_grid, lon_grid, lat, lon, plat, plon)

      ! Defaults
      implicit none

      ! Local variables
      integer, parameter                    :: nc=1440
      integer, parameter                    :: nr=720

      ! Arguments
      real, intent(in)  :: lon_grid(1:nc)
      real, intent(in)  :: lat_grid(1:nr)
      real, intent(in)  :: lat
      real, intent(in)  :: lon
      integer, intent(out) :: plat
      integer, intent(out) :: plon

      ! Local variables
      integer :: x,y

      do x=1, nc
         if (lon >= lon_grid(x)-0.125 .and. lon <= lon_grid(x)+0.125) then
            plon=x
            exit
         end if
      end do ! do x

      do y=1, nr
         if ((lat >= lat_grid(y)-0.125 .and. lat <= lat_grid(y)+0.125)) then
            plat=y
            exit
         end if
      end do ! do y

   end subroutine checkgrid


#if(defined USE_NETCDF4)   
   subroutine read_forestfraction(ff_filename, ff)

      ! Imports
      use LDT_logMod, only: LDT_verify
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      character (len=255), intent(in) :: ff_filename
      real, intent(out)           :: ff(1440,720)   ! fixed 1/4 deg

      ! Local variables
      integer                                          :: ncid, ff_varid

      ! open the file for reading
      call LDT_verify(nf90_open(path=trim(ff_filename), &
           mode=NF90_NOWRITE, ncid=ncid), &
           '[ERR] Error in nf90_open for '//trim(ff_filename))

      ! read forest fraction variable
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='Forest_fraction', varid=ff_varid), &
           '[ERR] Error in nf90_inq_varid for forest fraction')
      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=ff_varid, values=ff), &
           '[ERR] Error in nf90_get_var for forest fraction')

      ! close the file
      call LDT_verify(nf90_close(ncid), &
           '[ERR] Error in nf90_close for '//trim(ff_filename))

   end subroutine read_forestfraction

#else
   ! Dummy version
   subroutine read_forestfraction(ff_filename, ff)
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      implicit none
      character (len=255), intent(in) :: ff_filename
      real, intent(out)           :: ff(1440,720)   ! fixed 1/4 deg
      write(LDT_logunit,*) &
           "[ERR] Recompile LDT with netCDF4 support!"
      write(LDT_logunit,*) &
           "[ERR] Stopping in read_forestfraction..."
      call LDT_endrun()
   end subroutine read_forestfraction
#endif

   subroutine read_forestdensity(fd_filename, fd)
      ! Imports
      use LDT_logMod, only: LDT_verify
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      character (len=255), intent(in) :: fd_filename
      real, intent(out)           :: fd(1440,720)   ! fixed 1/4 deg

      ! Local variables
      integer                                          :: ncid, fd_varid

      ! open the file for reading
      call LDT_verify(nf90_open(path=trim(fd_filename), &
           mode=NF90_NOWRITE, ncid=ncid), &
           '[ERR] Error in nf90_open for '//trim(fd_filename))

      ! read forest density variable
      call LDT_verify(nf90_inq_varid(ncid=ncid, &
           name='PTC', varid=fd_varid), &
           '[ERR] Error in nf90_inq_varid for forest density')
      call LDT_verify(nf90_get_var(ncid=ncid, &
           varid=fd_varid, values=fd), &
           '[ERR] Error in nf90_get_var for forest density')

      ! close the file
      call LDT_verify(nf90_close(ncid), &
           '[ERR] Error in nf90_close for '//trim(fd_filename))

   end subroutine read_forestdensity


   subroutine search_files(date10, gmi_in)
      ! Imports
      use USAFSI_utilMod, only: date10_julhr, julhr_date10

      ! Defaults
      implicit none

      ! Arguments
      character*10,       intent(in) :: date10
      character*255,      intent(in) :: gmi_in

      ! Local variables
      integer            :: eof, n, i, j, k
      character(len=255)               :: file_path, cmd
      character*10                   :: date10_prev
      integer                        :: hr, st_hr, julhr
      character*2                    :: st_hr_str, cnt

      ! EMK
      character*12                   :: program_name          ! NAME OF CALLING PROGRAM
      character*12                   :: routine_name          ! NAME OF THIS ROUTINE

      ! define data values
      data routine_name  / 'search_files' /

      ! FIND THE DATE/TIME GROUP OF THE PREVIOUS CYCLE.
      call date10_julhr(date10, julhr, program_name, routine_name)
      !CALL DATE10_JULHR (DATE10, JULHR)
      julhr  = julhr  - 13
      call julhr_date10 (julhr, date10_prev, program_name, routine_name)
      !CALL JULHR_DATE10 (JULHR-13, DATE10_PREV)

      read(date10(9:10), '(I2)') hr

      n=1
      !do i=1, 3
         st_hr=hr - 13           ! start time: hr-13

         if (st_hr < 0 .and. hr >= 0) then
            st_hr = st_hr + 24
            do j=st_hr, 23
               write(st_hr_str,'(I0.2)') j ! convert int to string
               write(cnt,'(I0.2)') n ! convert int to string

               !file_path = trim(gmi_in)//'1C-R.GPM.GMI.XCAL2016-C.'//date10_prev(1:8)//'-S'//st_hr_str// &
               !            '*.RT-H5 2>/dev/null > file'//cnt//'.txt'

               file_path = trim(gmi_in)//'1C-R.GPM.GMI.XCAL2016-C.'//date10_prev(1:8)//'-S'//st_hr_str// &
                           '*.HDF5 2>/dev/null > file'//cnt//'.txt'

               cmd = 'ls '//file_path
               call system(cmd)
               n=n+1
            end do

            k=0  ! check start time
         else
            k=st_hr
         end if

         do j=k, hr
            write(st_hr_str,'(I0.2)') j ! convert int to string
            write(cnt,'(I0.2)') n ! convert int to string

            !file_path = trim(gmi_in)//'1C-R.GPM.GMI.XCAL2016-C.'//date10(1:8)//'-S'//st_hr_str// &
            !               '*.RT-H5 2>/dev/null > file'//cnt//'.txt'
            
            file_path = trim(gmi_in)//'1C-R.GPM.GMI.XCAL2016-C.'//date10(1:8)//'-S'//st_hr_str// &
                           '*.HDF5 2>/dev/null > file'//cnt//'.txt'

            cmd = 'ls '//file_path
            call system(cmd)
            n=n+1
         end do
      !end do

      call system ('ls file*.txt | xargs cat > ./gmi_filelist.txt')
      call system ('find . -type f -name "file*.txt" -print0 | xargs -0 rm -rf')
   
   end subroutine search_files

end module USAFSI_xcalgmiMod
