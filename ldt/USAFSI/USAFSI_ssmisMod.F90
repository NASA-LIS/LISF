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
! MODULE: USAFSI_ssmisMod
!
! REVISION HISTORY:
!  30 Dec 2018: Yeosang Yoon; Initial Implementation
!  03 Apr 2019: Yeosang Yoon; Update codes to fit LDT format
!  16 Apr 2019: Eric Kemp; Put code into module
!  09 May 2019: Eric Kemp; Renamed LDTSI
!  13 Dec 2019: Eric Kemp; Renamed USAFSI
!  09 Mar 2019: Eric Kemp; bug fixes and tweaks addressing mixed precision.
!
! DESCRIPTION:
! Source code for the retrieval of snow depth and sea ice concentration from
! SSMIS observations
!-------------------------------------------------------------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module USAFSI_ssmisMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: USAFSI_proc_ssmis

contains

   ! Public routine for processing SSMIS data
   subroutine USAFSI_proc_ssmis(date10, ssmis_in, ssmis, option)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
      use USAFSI_utilMod
      use USAFSI_paramsMod
      use map_utils

      ! Defaults
      implicit none

      ! Arguments
      character(len=10),  intent(in)  :: date10
      character(len=255), intent(in)  :: ssmis_in
      character(len=255), intent(in)  :: ssmis
      integer,            intent(in)  :: option

      ! Local variables
      integer                         :: eof, i, j, n, nArr, nFile, x, y
      character(len=255)              :: filename, nc_filename
      integer,           dimension(:), allocatable :: surflag0, surflag, &
           satid0, satid
      character(len=10), dimension(:), allocatable :: date0, date10_arr
      real(kind=8),      dimension(:), allocatable :: lat0,lon0,tb19h0, &
           tb19v0,tb37h0,tb37v0,tb22v0,tb91v0,tb91h0
      real(kind=8),      dimension(:), allocatable :: lat,lon,tb19h,tb19v,&
           tb37h,tb37v,tb22v,tb91v,tb91h,snowdepth,ct
      integer                                      :: hemi
      real                                         :: icoord, jcoord
      INTEGER, PARAMETER                           :: POSE    = 0   ! LONGITUDE ORIENTATION FLAG; 0 = POSITIVE WEST
      real :: rlat, rlon ! EMK
      type(proj_info) :: snodep_0p25deg_proj

      call search_files (date10, ssmis_in)
      open(unit=90, file='./ssmis_filelist.txt', form='formatted', &
           status='old', action='read')

      nFile = 0
      do
         read(90,'(A)',iostat=eof) filename
         if (eof/=0) exit

         write (LDT_logunit,*) '[INFO] get ssmis obs:', trim(filename)

         call read_bufr_attributes(filename, satid0, date0, lat0, lon0, &
              surflag0, tb19h0, tb19v0, tb37h0, tb37v0, tb22v0, tb91v0, tb91h0)

         ! resize the array considering the SSMIS TDR observations
         if(nFile==0) then
            call move_alloc(satid0, satid)
            call move_alloc(date0, date10_arr)
            call move_alloc(lat0, lat)
            call move_alloc(lon0, lon)
            call move_alloc(surflag0, surflag)
            !      call move_alloc(rainflag0, rainflag)
            call move_alloc(tb19h0, tb19h)
            call move_alloc(tb19v0, tb19v)
            call move_alloc(tb37h0, tb37h)
            call move_alloc(tb37v0, tb37v)
            call move_alloc(tb22v0, tb22v)
            call move_alloc(tb91v0, tb91v)
            call move_alloc(tb91h0, tb91h)
         else
            call resize_array_int(satid, satid0)
            call resize_array_str(date10_arr, date0)
            call resize_array(lat, lat0)
            call resize_array(lon, lon0)
            call resize_array_int(surflag, surflag0)
            !       call resize_array_int(rainflag, rainflag0)
            call resize_array(tb19h, tb19h0)
            call resize_array(tb19v, tb19v0)
            call resize_array(tb37h, tb37h0)
            call resize_array(tb37v, tb37v0)
            call resize_array(tb22v, tb22v0)
            call resize_array(tb91v, tb91v0)
            call resize_array(tb91h, tb91h0)
         endif

         nFile=nFile+1

      end do
      close(90)

      ! calculate snow depth
      nArr=size(tb19h)
      if(allocated(snowdepth)) deallocate(snowdepth)
      if(allocated(ct))        deallocate(ct)
      allocate(snowdepth(nArr))
      allocate(ct(nArr))

      write (LDT_logunit,*) '[INFO] calculate snow depth from ssmis'
      call calculate_snowdepth(nArr, surflag, tb19h, tb19v, tb37h, tb37v, &
           tb22v, tb91v, tb91h, snowdepth, option, lat, lon)

      write (LDT_logunit,*) '[INFO] calculate sea ice concentration from ssmis'
      call calculate_sea_ice_concentration(nArr, surflag, lat, lon, tb19h, &
           tb19v, tb37h, tb37v, ct)

      ! write file
      filename = trim(ssmis)//'ssmis_snoice_nh.06hr.'//date10//'.txt'
      write (LDT_logunit,*) &
           '[INFO] Writing SSMIS data to ', trim(filename)
      open(unit=10, file=filename,status='unknown', action='write')
      filename = trim(ssmis)//'ssmis_snoice_sh.06hr.'//date10//'.txt'
      write (LDT_logunit,*) &
           '[INFO] Writing SSMIS data to ', trim(filename)
      open(unit=20, file=filename,status='unknown', action='write')

      call map_set(proj_code=proj_latlon, &
           lat1=-89.875, &
           lon1=-179.875, &
           dx=0.25, &
           stdlon=0.25, &
           truelat1=0.25, &
           truelat2=0., &
           idim=igrid,&
           jdim=jgrid, &
           proj=snodep_0p25deg_proj)

      do i=1,size(lat)

         if (snowdepth(i) < 0) then      ! remove unrealistic values
            snowdepth(i)=-0.1
         end if

         ! lat, lon*100, snow depth(mm)*10
         !if (surflag(i)==0 .and. snowdepth(i)>=0) then      !  land points only
         if (ct(i) >=0 .or. snowdepth(i) >=0) then
            if (lat(i) > 0) then
               hemi=1
               rlat = real(lat(i))
               rlon = real(lon(i))
               !call LLTOPS (POSE, rlat, rlon, 16, hemi, icoord, jcoord)
               call latlon_to_ij(snodep_0p25deg_proj,rlat,rlon,icoord,jcoord)
               write(10, '(A10, I3, I6, I7, 2I5, 2I6)') date10_arr(i), &
                    satid(i), nint(real(lat(i))*100), nint(real(lon(i))*100), &
                    nint(icoord), nint(jcoord), nint(ct(i)), &
                    nint(real(snowdepth(i))*10)

            else
               hemi=2
               rlat = real(lat(i))
               rlon = real(lon(i))
               !call LLTOPS (POSE, rlat, rlon, 16, hemi, icoord, jcoord)
               call latlon_to_ij(snodep_0p25deg_proj,rlat,rlon,icoord,jcoord)
               write(20, '(A10, I3, I6, I7, 2I5, 2I6)') date10_arr(i), &
                    satid(i), nint(real(lat(i))*100), nint(real(lon(i))*100), &
                    nint(icoord), nint(jcoord), nint(ct(i)), &
                    nint(real(snowdepth(i))*10)

            end if ! if (lat(i) > 0) then
         end if
         !end if ! f (surflag(i) ==  0)
      end do ! do i=1,size(lat)
      close(10)
      close(20)

      ! write netCDF file for USAFSI
      nc_filename = trim(ssmis)//'ssmis_snoice_0p25deg.'//date10//'.nc'
      write (LDT_logunit,*) &
           '[INFO] Writing SSMIS data to ', trim(nc_filename)
      call write_netcdf(nc_filename, nArr, lat,lon, snowdepth, ct)

      if(allocated(satid0)) deallocate(satid0)
      if(allocated(date0)) deallocate(date0)
      if(allocated(lat0)) deallocate(lat0)
      if(allocated(lon0)) deallocate(lon0)
      if(allocated(surflag0)) deallocate(surflag0)
      if(allocated(tb19h0)) deallocate(tb19h0)
      if(allocated(tb19v0)) deallocate(tb19v0)
      if(allocated(tb37h0)) deallocate(tb37h0)
      if(allocated(tb37v0)) deallocate(tb37v0)
      if(allocated(tb22v0)) deallocate(tb22v0)
      if(allocated(tb91v0)) deallocate(tb91v0)
      if(allocated(tb91h0)) deallocate(tb91h0)

      if(allocated(satid)) deallocate(satid)
      if(allocated(date10_arr)) deallocate(date10_arr)
      if(allocated(lat)) deallocate(lat)
      if(allocated(lon)) deallocate(lon)
      if(allocated(surflag)) deallocate(surflag)
      if(allocated(tb19h)) deallocate(tb19h)
      if(allocated(tb19v)) deallocate(tb19v)
      if(allocated(tb37h)) deallocate(tb37h)
      if(allocated(tb37v)) deallocate(tb37v)
      if(allocated(tb22v)) deallocate(tb22v)
      if(allocated(tb91v)) deallocate(tb91v)
      if(allocated(tb91h)) deallocate(tb91h)

   end subroutine USAFSI_proc_ssmis

   ! *** Remaining routines are private ***

#if (defined USE_ECCODES)
   subroutine read_bufr_attributes(filename, satid, date10_arr, lat, lon, &
        surflag, tb19h, tb19v, tb37h, tb37v, tb22v, tb91v, tb91h)

      ! Imports
      use eccodes
      use LDT_logMod, only: ldt_logunit

      ! Defaults
      implicit none

      ! Arguments
      character(len=255), intent(in) :: filename
      integer, allocatable, intent(inout) :: satid(:)
      character(len=10), allocatable, intent(inout) :: date10_arr(:)
      real(kind=8), allocatable, intent(inout) :: lat(:)
      real(kind=8), allocatable, intent(inout) :: lon(:)
      integer, allocatable, intent(inout) :: surflag(:)
      real(kind=8), allocatable, intent(inout) :: tb19h(:)
      real(kind=8), allocatable, intent(inout) :: tb19v(:)
      real(kind=8), allocatable, intent(inout) :: tb37h(:)
      real(kind=8), allocatable, intent(inout) :: tb37v(:)
      real(kind=8), allocatable, intent(inout) :: tb22v(:)
      real(kind=8), allocatable, intent(inout) :: tb91v(:)
      real(kind=8), allocatable, intent(inout) :: tb91h(:)

      ! Local variables
      integer                        :: ifile, iret, ibufr
      integer                        :: i, option, n
      integer                        :: numObs
      integer                        :: yyyy, mm, dd, hh, id
      character(len=10)              :: sbufr
      character(len=10), dimension(:), allocatable       :: sval
      integer,           dimension(:), allocatable       :: ival
      real(kind=8),              dimension(:), allocatable       :: rval

      !----- check data size
      call codes_open_file(ifile,filename,'r')
      call codes_bufr_new_from_file(ifile,ibufr,iret)

      n=0
      do while (iret/=CODES_END_OF_FILE)
         call codes_get(ibufr,'numberOfSubsets',numObs)    ! Read the total number of subsets.
         call codes_release(ibufr)
         call codes_bufr_new_from_file(ifile,ibufr,iret)
         n=n+1
      end do
      call codes_close_file(ifile) ! Close file
      !-------------------------------------------------

      if(allocated(satid)) deallocate(satid)
      if(allocated(date10_arr)) deallocate(date10_arr)
      if(allocated(lat)) deallocate(lat)
      if(allocated(lon)) deallocate(lon)
      if(allocated(surflag)) deallocate(surflag)
      if(allocated(tb19h)) deallocate(tb19h)
      if(allocated(tb19v)) deallocate(tb19v)
      if(allocated(tb37h)) deallocate(tb37h)
      if(allocated(tb37v)) deallocate(tb37v)
      if(allocated(tb22v)) deallocate(tb22v)
      if(allocated(tb91v)) deallocate(tb91v)
      if(allocated(tb91h)) deallocate(tb91h)
      allocate(satid(numObs*n))
      allocate(date10_arr(numObs*n))
      allocate(lat(numObs*n))
      allocate(lon(numObs*n))
      allocate(surflag(numObs*n))
      allocate(tb19h(numObs*n))
      allocate(tb19v(numObs*n))
      allocate(tb37h(numObs*n))
      allocate(tb37v(numObs*n))
      allocate(tb22v(numObs*n))
      allocate(tb91v(numObs*n))
      allocate(tb91h(numObs*n))

      ! real actual dataset
      call codes_open_file(ifile,trim(filename),'r')
      call codes_bufr_new_from_file(ifile,ibufr,iret)

      n=1;
      do while (iret/=CODES_END_OF_FILE)
         call codes_set(ibufr,"unpack",1);    ! i.e. unpack the data values

         if(allocated(sval)) deallocate(sval)
         if(allocated(ival)) deallocate(ival)
         if(allocated(rval)) deallocate(rval)
         allocate(sval(numObs))
         allocate(ival(numObs))
         allocate(rval(numObs))

         call codes_get(ibufr,'satelliteIdentifier',id)
         if (id==249) then
            ival(1:numObs)=16
         else if (id==285) then
            ival(1:numObs)=17
         else
            ival(1:numObs)=18
         end if
         satid(1+numObs*(n-1):numObs*n)=ival

         call codes_get(ibufr,'#1#year',yyyy)  ! Get time
         call codes_get(ibufr,'#1#month',mm)
         call codes_get(ibufr,'#1#day',dd)
         call codes_get(ibufr,'#1#hour',hh)

         write(sbufr(1:4),'(I4)') yyyy
         write(sbufr(5:6),'(I0.2)') mm
         write(sbufr(7:8),'(I0.2)') dd
         write(sbufr(9:10),'(I0.2)') hh

         sval(1:numObs)=sbufr
         date10_arr(1+numObs*(n-1):numObs*n)=sval

         call codes_get(ibufr,'#1#latitude',rval)         ! Get latitude
         lat(1+numObs*(n-1):numObs*n)=rval
         call codes_get(ibufr,'#1#longitude',rval)        ! Get longitude
         lon(1+numObs*(n-1):numObs*n)=rval

         !Get SurfaceFlag (0:Land, 1:Reserved, 2:Near coast, 3:Ice, 4:Possilbe ice, 5:Ocean, 6:Coast,
         !                  7-14: Reserved, 15: Missing value)
         call codes_get(ibufr,'surfaceFlag',ival)         ! Get SurfaceFlag
         surflag(1+numObs*(n-1):numObs*n)=ival

         call codes_get(ibufr,'#12#brightnessTemperature',rval) ! Get Tb #12: 19.35H
         tb19h(1+numObs*(n-1):numObs*n)=rval
         call codes_get(ibufr,'#13#brightnessTemperature',rval) ! Get Tb #13: 19.35V
         tb19v(1+numObs*(n-1):numObs*n)=rval

         !! EMK TEST
         !do i = 1+numObs*(n-1), numObs*n
         !   if (tb91v(i) < 0) then
         !      write(ldt_logunit,*) &
         !           '[WARN] i, sbufr, satid, lat, lon, tb91v = ',&
         !           i, sbufr, satid(i), lat(i), lon(i), tb91v(i)
         !   end if
         !end do

         call codes_get(ibufr,'#14#brightnessTemperature',rval) ! Get Tb #14: 22.235V
         tb22v(1+numObs*(n-1):numObs*n)=rval
         call codes_get(ibufr,'#15#brightnessTemperature',rval) ! Get Tb #15: 37.0H
         tb37h(1+numObs*(n-1):numObs*n)=rval
         call codes_get(ibufr,'#16#brightnessTemperature',rval) ! Get Tb #16: 37.0V
         tb37v(1+numObs*(n-1):numObs*n)=rval
         call codes_get(ibufr,'#17#brightnessTemperature',rval) ! Get Tb #17: 91.6550V
         tb91v(1+numObs*(n-1):numObs*n)=rval
         call codes_get(ibufr,'#18#brightnessTemperature',rval) ! Get Tb #18: 91.6550H
         tb91h(1+numObs*(n-1):numObs*n)=rval

         call codes_release(ibufr)  ! Release the bufr message
         call codes_bufr_new_from_file(ifile,ibufr,iret) ! Load the next bufr message

         n=n+1

      end do

      call codes_close_file(ifile) ! Close file

      ! Clean up
      if(allocated(sval)) deallocate(sval)
      if(allocated(ival)) deallocate(ival)
      if(allocated(rval)) deallocate(rval)

   end subroutine read_bufr_attributes

#else
   ! Dummy version if no ECCODES support was compiled
   subroutine read_bufr_attributes(filename, satid, date10_arr, lat, lon, &
        surflag, tb19h, tb19v, tb37h, tb37v, tb22v, tb91v, tb91h)
      use LDT_logMod, only: LDT_endrun, LDT_logunit
      implicit none
      character(len=255), intent(in) :: filename
      integer, allocatable, intent(inout) :: satid(:)
      character(len=10), allocatable, intent(inout) :: date10_arr(:)
      real(kind=8), allocatable, intent(inout) :: lat(:)
      real(kind=8), allocatable, intent(inout) :: lon(:)
      integer, allocatable, intent(inout) :: surflag(:)
      real(kind=8), allocatable, intent(inout) :: tb19h(:)
      real(kind=8), allocatable, intent(inout) :: tb19v(:)
      real(kind=8), allocatable, intent(inout) :: tb37h(:)
      real(kind=8), allocatable, intent(inout) :: tb37v(:)
      real(kind=8), allocatable, intent(inout) :: tb22v(:)
      real(kind=8), allocatable, intent(inout) :: tb91v(:)
      real(kind=8), allocatable, intent(inout) :: tb91h(:)
      write(LDT_logunit,*) &
           "[ERR] Recompile LDT with ECCODES support and run again!"
      write(LDT_logunit,*) &
           "[ERR] Stopping in read_bufr_attributes..."
      call LDT_endrun()
   end subroutine read_bufr_attributes
#endif

   subroutine calculate_snowdepth(n, surflag, tb19h, tb19v, tb37h, tb37v, &
        tb22v,tb91v,tb91h, snowdepth, option, lat, lon)

      ! option == 1: Hollinger, 1991,    SD=4445.0-17.95TB_37V
      ! option == 2: Markus (Chang et al., 1987), SD=1.58(TB_19H-TB_37H)
      ! option == 3; Foster et al., 1997

      ! Imports
      use LDT_logMod, only: ldt_logunit
      use LDT_usafsiMod, only: usafsi_settings

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in)   :: n
      integer, intent(in)   :: surflag(:)
      real(kind=8), intent(in) :: tb19h(:)
      real(kind=8), intent(in) :: tb19v(:)
      real(kind=8), intent(in) :: tb37h(:)
      real(kind=8), intent(in) :: tb37v(:)
      real(kind=8), intent(in) :: tb22v(:)
      real(kind=8), intent(in) :: tb91h(:)
      real(kind=8), intent(in) :: tb91v(:)
      real(kind=8), intent(out)  :: snowdepth(1:n)
      integer, intent(in)   :: option
      real(kind=8), intent(in) :: lat(:)
      real(kind=8), intent(in) :: lon(:)

      ! Local constants
      integer, parameter                             :: nc=1440
      integer, parameter                             :: nr=720

      ! Local variables
      integer                                        :: i
      character(len=255)                             :: ff_filename
      real(kind=8)                                   :: pd19,pd91,tt,si91, &
           scat,sc37,sc91,scx
      logical                                        :: flag
      real                                           :: lon_grid(nc), &
           lat_grid(nr), ratio
      real                                           :: ff(nc,nr)
      integer                                        :: plat, plon

      ! set 0.25 deg grid, for forest fraction
      if (option==3) then
         do i=1, nc
            lon_grid(i)=-179.875+0.25*(i-1)
         end do
         do i=1, nr
            lat_grid(i)=-89.875+0.25*(i-1)
         end do

         ff_filename=trim(usafsi_settings%ff_file)
         ! read forest fraction
         call read_forestfraction(ff_filename,ff)
      end if

      do i=1,n

         ! EMK Sanity check
         if (tb91v(i) < 0) then
            !write(ldt_logunit,*) &
            !     '[WARN] Found negative tb91v for i = ', i
            snowdepth(i) = -0.1
            cycle
         end if

         ! quality check
         pd19=tb19v(i)-tb19h(i)
         pd91=tb91v(i)-tb91h(i)
         tt=169d0+0.5d0*tb91v(i)
         sc91=tb22v(i)-tb91v(i)-3.0d0
         sc37=tb19v(i)-tb37v(i)-3.0d0
         scx=tb37v(i)-tb91v(i)-1.0d0
         scat=sc91
         if (sc37>scat) then
            scat=sc37
         end if
         if (scx>scat) then
            scat=scx
         end if

         ! NG, 2002
         if (scat > 0) then ! scattering materials
            flag=.true.
            if ((.not. real(tb22v(i)) < 260) .and. (real(pd91) > 3) .and. (.not. real(scat) > 3)) then  ! cold rain
               flag=.false.
            end if
            if ((.not. real(tb22v(i)) < 264) .or. (.not. real(tb22v(i)) < real(tt))) then  ! other rain events
               flag=.false.
            end if
            if ((.not. real(pd19) < 18) .and. (.not. real(sc37) > 10) .and. (.not. real(scx) > 10)) then ! cold desert
               flag=.false.
            end if
            if ((.not. real(pd19) < 8) .and. (.not. real(sc91) > 6) .and. (.not. real(sc37) > 2)) then ! frozen ground
               flag=.false.
            end if
            if ((.not. real(tb22v(i)) > 216) .or. ((.not. real(tb22v(i)) > 235) .and. (.not. real(pd19) < 23))) then ! gracier
               flag=.false.
            end if
         else
            flag=.false.
         end if

         !Get SurfaceFlag (0:Land, 1:Reserved, 2:Near coast, 3:Ice, 4:Possilbe ice, 5:Ocean, 6:Coast,
         !                  7-14: Reserved, 15: Missing value)
         if ((flag .eqv. .true.) .and. (surflag(i)==0 .or. surflag(i)==2)) then
            if (option==1) then
               snowdepth(i)=4445.0-17.95*tb37v(i)
            else if (option==2) then
               snowdepth(i)=1.58*(tb19h(i)-tb37h(i))*10.0  ! unit cm -> mm
               !if (snowdepth(i) <= 25) then   ! derived snow depth less than 2.5 cm will be assigned as no snow.
               !  snowdepth(i) = -0.1
               !end if
            else
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
            end if ! if (option==1) then
         else
            snowdepth(i) = -0.1
         end if
      end do

   end subroutine calculate_snowdepth

   subroutine calculate_sea_ice_concentration(n, surflag, lat, lon, tb19h, &
        tb19v, tb37h, tb37v, ct)

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in)   :: n
      integer, intent(in)   :: surflag(:)
      real(kind=8), intent(in) :: lat(:)
      real(kind=8), intent(in) :: lon(:)
      real(kind=8), intent(in) :: tb19h(:)
      real(kind=8), intent(in) :: tb19v(:)
      real(kind=8), intent(in) :: tb37h(:)
      real(kind=8), intent(in) :: tb37v(:)
      real(kind=8), intent(out) :: ct(:)

      ! Local variables
      integer                                      :: i
      real(kind=8)                                 :: pr, gr, d, cf, cm

      do i=1,n
         if (real(lat(i))>44.5 .or. real(lat(i))<-52) then

            !Get SurfaceFlag (0:Land, 1:Reserved, 2:Near coast, 3:Ice, 4:Possilbe ice, 5:Ocean, 6:Coast,
            !                  7-14: Reserved, 15: Missing value)
            !if (surflag(i)==2 .or. surflag(i)==3 .or. surflag(i)==4 .or. surflag(i)==5 .or. surflag(i)==6) then
            if (surflag(i)==3 .or. surflag(i)==4) then
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
            end if
         else
            ct(i)=-1
         end if

      end do

   end subroutine calculate_sea_ice_concentration

   subroutine resize_array(arr1, arr2)

      ! Defaults
      implicit none

      ! Arguments
      real(kind=8), allocatable, intent(inout) :: arr1(:)
      real(kind=8), allocatable, intent(in)    :: arr2(:)

      ! Local variables
      real(kind=8), dimension(:), allocatable                :: temp
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
      integer, intent(in)       :: n
      real(kind=8), intent(in) :: lat(1:n)
      real(kind=8), intent(in) :: lon(1:n)
      real(kind=8), intent(in) :: snowdepth(1:n)
      real(kind=8), intent(in) :: ct(1:n)

      ! Local constants
      integer, parameter                             :: nc=1440 ! fixed 1/4 deg
      integer, parameter                             :: nr=720
      real, parameter                                :: fillval = -1.0

      ! Local variables
      real, dimension(1:nc)                          :: lon_grid
      real, dimension(1:nr)                          :: lat_grid
      real, dimension(1:nc,1:nr)                     :: smidep, ndep, smicon, &
           ncon
      integer                                        :: plat, plon
      integer                                        :: iret, ncid, &
           smidep_varid, smicon_varid, lat_varid, lon_varid, dim_ids(2)
      integer :: i,x,y

      ! set 0.25 deg grid
      do i=1, nc
         lon_grid(i)=-179.875+0.25*(i-1)
      end do
      do i=1, nr
         lat_grid(i)=-89.875+0.25*(i-1)
      end do

      ! initialization
      smidep(1:nc,1:nr)=0.0
      ndep(1:nc,1:nr)=0.0
      smicon(1:nc,1:nr)=0.0
      ncon(1:nc,1:nr)=0.0

      do i=1,n

         ! find nearest grid location
         call checkgrid(lat_grid,lon_grid,lat(i),lon(i),plat,plon)

         ! for snow depth
         if (snowdepth(i)>=0) then
            smidep(plon,plat)=smidep(plon,plat)+snowdepth(i)
            ndep(plon,plat)=ndep(plon,plat)+1
         end if

         ! for ice concentration
         if (ct(i)>=0) then
            smicon(plon,plat)=smicon(plon,plat)+ct(i)
            ncon(plon,plat)=ncon(plon,plat)+1
         end if

      end do ! i=1,size(lat)

      ! Averaging snow depth & ice concentration
      do x=1,nc
         do y=1,nr
            if (ndep(x,y)>=1) then
               smidep(x,y)=smidep(x,y)/ndep(x,y)
            else
               smidep(x,y)=fillval   ! apply null vaule
            end if

            if (ncon(x,y)>=1) then
               smicon(x,y)=smicon(x,y)/ncon(x,y)
            else
               smicon(x,y)=fillval   ! apply null vaule
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

      ! ssmis snow depth
      call LDT_verify(nf90_def_var(ncid,"smidep",nf90_float, &
           dimids=dim_ids, varid=smidep_varid),'[ERR] nf90_def_var failed')
      call LDT_verify(nf90_put_att(ncid,smidep_varid, &
           "units","mm"),'[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,smidep_varid, &
           "long_name","ssmis snow depth of surface snow over land"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,smidep_varid, &
           "standard_name","ssmis_snow_depth_of_surface_snow"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,smidep_varid, &
           '_FillValue',fillval), &
           '[ERR] nf90_put_att failed for SMIDEP')

      ! sea ice concentration analysis
      call LDT_verify(nf90_def_var(ncid,"smicon",nf90_float, &
           dimids=dim_ids, varid=smicon_varid),'[ERR] nf90_def_var failed')
      call LDT_verify(nf90_put_att(ncid,smicon_varid, &
           "units","%"),'[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,smicon_varid, &
           "long_name","concentration of sea ice (0-100)"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,smicon_varid, &
           "standard_name","sea_ice_area_fraction"),&
           '[ERR] nf90_put_att failed')
      call LDT_verify(nf90_put_att(ncid,smicon_varid, &
           'missing_value',fillval), &
           '[ERR] nf90_put_att failed for SMICON')
      call LDT_verify(nf90_put_att(ncid,smicon_varid, &
           '_FillValue',fillval), &
           '[ERR] nf90_put_att failed for SMICON')
      call LDT_verify(nf90_put_att(ncid,smicon_varid, &
           'valid_range',(/0.,100./)), &
           '[ERR] nf90_put_att failed for SMICON')

      ! Miscellaneous header information
      call LDT_verify(nf90_put_att(ncid,nf90_global,"title", &
           "LDT SSMIS snow depth analysis"), &
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

      ! Write the SSMIS snow depth/ice concentration fields
      call LDT_verify(nf90_put_var(ncid,smidep_varid,&
           smidep(:,:), (/1,1/),(/nc,nr/)), &
           '[ERR] nf90_put_var failed for smidep')
      call LDT_verify(nf90_put_var(ncid,smicon_varid,smicon(:,:), &
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
      integer, intent(in)       :: n
      real(kind=8), intent(in) :: lat(1:n)
      real(kind=8), intent(in) :: lon(1:n)
      real(kind=8), intent(in) :: snowdepth(1:n)
      real(kind=8), intent(in) :: ct(1:n)
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
      real(kind=8), intent(in)  :: lat
      real(kind=8), intent(in)  :: lon
      integer, intent(out) :: plat
      integer, intent(out) :: plon

      ! Local variables
      integer :: x,y
      real :: rlat, rlon

      rlat = real(lat)
      rlon = real(lon)

      do x=1, nc
         if ((.not. rlon < lon_grid(x)-0.125) .and. (.not. rlon > lon_grid(x)+0.125)) then

            plon=x
            exit
         end if
      end do ! do x

      do y=1, nr
         if ((.not. rlat < lat_grid(y)-0.125) .and. (.not. rlat > lat_grid(y)+0.125)) then
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

   subroutine search_files(date10, ssmis_in)

      ! Imports
      use USAFSI_utilMod, only: date10_julhr, julhr_date10

      ! Defaults
      implicit none

      ! Arguments
      character*10,       intent(in) :: date10
      character*255,      intent(in) :: ssmis_in

      ! Local variables
      integer            :: eof, n, i, j, k
      character(len=255)               :: file_path, cmd
      character*10                   :: date10_prev
      integer                        :: hr, st_hr, julhr
      character*2                    :: satid (3)
      character*2                    :: st_hr_str, cnt

      ! EMK
      character*12                   :: program_name          ! NAME OF CALLING PROGRAM
      character*12                   :: routine_name          ! NAME OF THIS ROUTINE

      ! define data values
      data satid            / '16', '17', '18' /
      data routine_name  / 'search_files' /

      ! FIND THE DATE/TIME GROUP OF THE PREVIOUS CYCLE.
      call date10_julhr(date10, julhr, program_name, routine_name)
      !CALL DATE10_JULHR (DATE10, JULHR)
      julhr  = julhr  - 13
      call julhr_date10 (julhr, date10_prev, program_name, routine_name)
      !CALL JULHR_DATE10 (JULHR-13, DATE10_PREV)

      read(date10(9:10), '(I2)') hr

      n=1
      do i=1, 3
         st_hr=hr - 13           ! start time: hr-13

         if (st_hr < 0 .and. hr >= 0) then
            st_hr = st_hr + 24
            do j=st_hr, 23
               write(st_hr_str,'(I0.2)') j ! convert int to string
               write(cnt,'(I0.2)') n ! convert int to string

               file_path = trim(ssmis_in)//'tdrcr_f'//satid(i)//'_d'//date10_prev(1:8)//'_s'//st_hr_str// &
                    '* 2>/dev/null > file'//cnt//'.txt'

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

            file_path = trim(ssmis_in)//'tdrcr_f'//satid(i)//'_d'//date10(1:8)//'_s'//st_hr_str// &
                 '* 2>/dev/null > file'//cnt//'.txt'

            cmd = 'ls '//file_path
            call system(cmd)
            n=n+1
         end do
      end do

      call system ('ls file*.txt | xargs cat > ./ssmis_filelist.txt')
      call system ('find . -type f -name "file*.txt" -print0 | xargs -0 rm -rf')

   end subroutine search_files

end module USAFSI_ssmisMod
