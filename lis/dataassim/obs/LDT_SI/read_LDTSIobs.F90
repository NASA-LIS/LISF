!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"
#include "LIS_NetCDF_inc.h"

! Read the LDT SI observations and process for later use within the 
! DA algorithm
! TODO:  Wrap this in a module
subroutine read_LDTSIobs(n, k, OBS_State, OBS_Pert_State)

   ! Imports
   use ESMF
   use LDTSIobs_Mod, only: LDTSI_obs
   use LIS_coreMod, only: LIS_rc
   use LIS_DAobservationsMod, only: LIS_obs_domain
   use LIS_logMod, only: LIS_verify, LIS_logunit
   use LIS_timeMgrMod, only: LIS_isAlarmRinging, LIS_clock

   ! Defaults
   implicit none

   ! Arguments
   integer, intent(in) :: n 
   integer, intent(in) :: k
   type(ESMF_State), intent(inout)    :: OBS_State
   type(ESMF_State), intent(inout)    :: OBS_Pert_State

   ! Local variables
   logical             :: alarmCheck
   integer, allocatable :: assimflag(:)
   integer             :: c
   type(ESMF_Time)     :: currTime
   logical             :: data_update
   integer             :: dd
   type(ESMF_TimeInterval) :: deltaT
   logical             :: file_exists
   character(255)      :: filename
   integer, allocatable :: gid(:)
   integer             :: hh
   integer             :: ierr
   integer             :: mm
   integer             :: nc
   integer             :: nr
   character(100)      :: obsdir
   real,    pointer    :: obsl(:)
   integer             :: r
   type(ESMF_Field)    :: snowfield
   real, allocatable   :: snoanl(:,:)
   integer             :: status
   integer             :: t
   real, allocatable   :: varfield(:,:)
   integer             :: yyyy

   call ESMF_AttributeGet(OBS_State, "Data Directory",&
        obsdir, rc=status)
   call LIS_verify(status)
   call ESMF_AttributeGet(OBS_State, "Data Update Status",&
       data_update, rc=status)
   call LIS_verify(status)

   file_exists = .false.

   alarmCheck = LIS_isAlarmRinging(LIS_rc, "LDTSI read alarm")

   if (alarmCheck) then
      call ESMF_TimeIntervalSet(deltaT, s=nint(LIS_rc%ts))
      call ESMF_ClockGet(LIS_clock, currTime=currTime, rc=status)       
      currTime = currTime + deltaT

      call ESMF_TimeGet(currTime, yy=yyyy, mm=mm, dd=dd, h=hh, rc=status)
      call LDTSI_filename(filename, obsdir, yyyy, mm, dd, hh)

      inquire(file=trim(filename), exist=file_exists)

      if (.not. file_exists) then
         write(LIS_logunit,*)'[WARN] Cannot find file ',trim(filename)
      end if

   end if

   ! Jump out if we have no file to process
   if (.not. file_exists) then
      call ESMF_AttributeSet(OBS_State, "Data Update Status",&
          .false., rc=status)
      call LIS_verify(status)
      return
   end if

   ! Beyond this point we have a file to process
   call ESMF_StateGet(OBS_State, "Observation01", snowfield,&
        rc=status)
   call LIS_verify(status)

   call ESMF_FieldGet(snowfield, localDE=0, farrayPtr=obsl, rc=status)
   call LIS_verify(status)

   ! Pull the snow depth from the file
   call LDTSI_reader(filename, n, k, snoanl, ierr)   
   if (ierr .ne. 0) then
      write(LIS_logunit,*) &
           "[WARN] Could not read from LDTSI file ", trim(filename)
      call ESMF_AttributeSet(OBS_State, "Data Update Status",&
          .false., rc=status)
      call LIS_verify(status)
      return      
   end if

   ! Interpolate the data
   allocate(varfield(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k)))
   call interp_LDTSIfield(n, k, snoanl, LIS_rc%udef, varfield)

   do r = 1, LIS_rc%obs_lnr(n)
      do c = 1, LIS_rc%obs_lnc(n)
         if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1) then
            obsl(LIS_obs_domain(n,k)%gindex(c,r)) = varfield(c,r)
         end if
      end do ! c
   end do ! r
   
   call ESMF_AttributeSet(OBS_State, "Data Update Status",&
        .true. , rc=status)
   call LIS_verify(status, 'Problem with Data Update Status')

   allocate(gid(LIS_rc%obs_ngrid(k)))
   allocate(assimflag(LIS_rc%obs_ngrid(k)))
   do t = 1, LIS_rc%obs_ngrid(k)
      gid(t) = t
      if(obsl(t).ne.-9999.0) then 
         assimflag(t) = 1
      else
         assimflag(t) = 0
      endif
   end do ! t

   if(LIS_rc%obs_ngrid(k).gt.0) then 
      call ESMF_AttributeSet(snowfield, "Grid Number",&
           gid, itemCount=LIS_rc%obs_ngrid(k), rc=status)
      call LIS_verify(status)
      
      call ESMF_AttributeSet(snowfield, "Assimilation Flag",&
           assimflag, itemCount=LIS_rc%obs_ngrid(k), rc=status)
      call LIS_verify(status)
   endif

   ! Clean up
   deallocate(assimflag)
   deallocate(gid)
   deallocate(snoanl)
   deallocate(varfield)
   deallocate(LDTSI_obs(n)%rlat1)
   deallocate(LDTSI_obs(n)%rlon1)
   deallocate(LDTSI_obs(n)%n111)

contains

   ! Constructs LDTSI filename
   subroutine LDTSI_filename(filename, dir, yyyy, mm, dd, hh)
      implicit none
      character(255), intent(inout) :: filename
      character(100), intent(in) :: dir
      integer, intent(in) :: yyyy
      integer, intent(in) :: mm
      integer, intent(in) :: dd
      integer, intent(in) :: hh
      character(4) :: cyyyy
      character(2) :: cmm, cdd, chh
      write(unit=cyyyy, fmt='(i4.4)') yyyy
      write(unit=cmm, fmt='(i2.2)') mm
      write(unit=cdd, fmt='(i2.2)') dd
      write(unit=chh, fmt='(i2.2)') hh
      filename = trim(dir) // "/" &
           // "ldtsi_" &
           // trim(cyyyy) // trim(cmm) // trim(cdd) // trim(chh) &
           // ".nc"           
   end subroutine LDTSI_filename
   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   
   ! Reads LDTSI netcdf file
   subroutine LDTSI_reader(filename, n, k, snoanl, ierr)
      
      ! Imports 
      use LDTSIobs_Mod, only: LDTSI_obs
      use LIS_coreMod, only: LIS_rc
      use LIS_logMod, only: LIS_logunit, LIS_verify, LIS_endrun
      use netcdf
      
      ! Defaults
      implicit none
      
      ! Arguments
      character(255), intent(in) :: filename
      integer, intent(in) :: n
      integer, intent(in) :: k
      real, allocatable, intent(out) :: snoanl(:,:)
      integer, intent(out) :: ierr
      
      ! Local variables
      integer :: c
      integer :: dim_ids(3)
      logical :: diff
      real :: dx
      real :: dy
      integer :: i
      character(255) :: map_projection_name
      integer :: ncid
      integer :: nlat
      integer :: nlon
      real :: north_east_corner_lat
      real :: north_east_corner_lon
      integer :: ntime
      integer :: r
      integer :: snoanl_varid
      real :: south_west_corner_lat
      real :: south_west_corner_lon
      real, allocatable :: tmp(:,:,:)
      real :: tmp_griddesci(50)
      
      ierr = 1 ! Change this below
      
      write(LIS_logunit,*) &
           "[INFO] Reading LDTSI file ", trim(filename)
      
      ! Open the file for reading
      call LIS_verify(nf90_open(path=trim(filename), &
           mode=NF90_NOWRITE, &
           ncid=ncid), &
           '[ERR] Error in nf90_open for '//trim(filename))
      
      ! Get the dimension IDs
      call LIS_verify(nf90_inq_dimid(ncid=ncid,&
           name='time',&
           dimid=dim_ids(3)), &
           '[ERR] Error in nf90_inq_dimid for dimension time')
      call LIS_verify(nf90_inq_dimid(ncid=ncid,&
           name='lat',&
           dimid=dim_ids(2)), &
           '[ERR] Error in nf90_inq_dimid for dimension lat')
      call LIS_verify(nf90_inq_dimid(ncid=ncid,&
           name='lon',&
           dimid=dim_ids(1)), &
           '[ERR] Error in nf90_inq_dimid for dimension lon')
      
      ! Get the actual dimension sizes
      call LIS_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(3), &
           len=ntime), &
           '[ERR] Error in nf90_inquire_dimension for dimension time')
      call LIS_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(2), &
           len=nlat), &
           '[ERR] Error in nf90_inquire_dimension for dimension lat')
      call LIS_verify(nf90_inquire_dimension(ncid=ncid, &
           dimid=dim_ids(1), &
           len=nlon), &
           '[ERR] Error in nf90_inquire_dimension for dimension lon')

      ! Sanity checks
      if (ntime .ne. 1) then
         write(LIS_logunit, *) &
              "[ERR] Dimension mismatch"
         write(LIS_logunit, *) &
              "[ERR] Expected time = 1, found time = ", ntime
         call LIS_endrun()
      end if

      ! Fetch the snoanl varid
      call LIS_verify(nf90_inq_varid(ncid=ncid, &
           name='snoanl', &
           varid=snoanl_varid), &
           '[ERR] Error in nf90_inq_varid for snoanl')

      ! Read the snoanl variable
      allocate(tmp(nlon,nlat,1))
      call LIS_verify(nf90_get_var(ncid=ncid, &
           varid=snoanl_varid, &
           values=tmp), &
           '[ERR] Error in nf90_get_var for snoanl')
      allocate(snoanl(nlon,nlat))
      do r = 1, nlat
         do c = 1, nlon
            if (tmp(c,r,1) < 0) then
               snoanl(c,r) = LIS_rc%udef
            else
               snoanl(c,r) = tmp(c,r,1)
            end if
         end do ! c
      end do ! r
      deallocate(tmp)

      ! Fetch the map projection name 
      call LIS_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name="MAP_PROJECTION", &
           values=map_projection_name), &
           '[ERR] Error in nf90_get_att for MAP_PROJECTION')
      if (trim(map_projection_name) .ne. "EQUIDISTANT CYLINDRICAL") then
         write(LIS_logunit,*) &
              "[ERR] Unsupported map projection for LDTSI product!"
         write(LIS_logunit,*) &
              "[ERR] Expected EQUIDISTANT CYLINDRICAL, found ", &
              trim(map_projection_name)
         call LIS_endrun()
      end if

      ! Fetch the LDTSI data resolution
      call LIS_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name="DX", &
           values=dx), &
           '[ERR] Error in nf90_get_att for DX')
      call LIS_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name="DY", &
           values=dy), &
           '[ERR] Error in nf90_get_att for DY')

      ! Fetch the southwest lat/lon
      call LIS_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name="SOUTH_WEST_CORNER_LAT", &
           values=south_west_corner_lat), &
           '[ERR] Error in nf90_get_att for SOUTH_WEST_CORNER_LAT')
      call LIS_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name="SOUTH_WEST_CORNER_LON", &
           values=south_west_corner_lon), &
           '[ERR] Error in nf90_get_att for SOUTH_WEST_CORNER_LON')

      ! Fetch the northeast lat/lon
      call LIS_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name="NORTH_EAST_CORNER_LAT", &
           values=north_east_corner_lat), &
           '[ERR] Error in nf90_get_att for NORTH_EAST_CORNER_LAT')
      call LIS_verify(nf90_get_att(ncid=ncid, &
           varid=NF90_GLOBAL, &
           name="NORTH_EAST_CORNER_LON", &
           values=north_east_corner_lon), &
           '[ERR] Error in nf90_get_att for NORTH_EAST_CORNER_LON')

      ! Close the file
      call LIS_verify(nf90_close(ncid), &
           '[ERR] Error in nf90_close for '//trim(filename))

      ! Populate the griddesci array for this file
      tmp_griddesci(:) = 0
      tmp_griddesci(1) = 0 ! Lat/lon
      tmp_griddesci(2) = real(nlon)
      tmp_griddesci(3) = real(nlat)
      tmp_griddesci(4) = south_west_corner_lat
      tmp_griddesci(5) = south_west_corner_lon
      tmp_griddesci(6) = 128
      tmp_griddesci(7) = north_east_corner_lat
      tmp_griddesci(8) = north_east_corner_lon
      tmp_griddesci(9) = dx
      tmp_griddesci(10) = dy
      tmp_griddesci(11) = 64
      tmp_griddesci(20) = 64

      LDTSI_obs(n)%nc_ldt = nlon
      LDTSI_obs(n)%nr_ldt = nlat
      LDTSI_obs(n)%mi = nlon*nlat ! On LDT grid

      ! These are on LIS grid
      LDTSI_obs(n)%nc_lis = &
           nint((LIS_rc%obs_gridDesc(k,8) - LIS_rc%obs_gridDesc(k,5)) / &
                 LIS_rc%obs_gridDesc(k,9)) + 1
      LDTSI_obs(n)%nr_lis = &
           nint((LIS_rc%obs_gridDesc(k,7) - LIS_rc%obs_gridDesc(k,4)) / &
                 LIS_rc%obs_gridDesc(k,10)) + 1
      LDTSI_obs(n)%mo1 = LDTSI_obs(n)%nc_lis * LDTSI_obs(n)%nr_lis

      LDTSI_obs(n)%gridDesco(1) = 0
      LDTSI_obs(n)%gridDesco(2) = LDTSI_obs(n)%nc_lis
      LDTSI_obs(n)%gridDesco(3) = LDTSI_obs(n)%nr_lis
      LDTSI_obs(n)%gridDesco(4) = LIS_rc%obs_gridDesc(k,4)
      LDTSI_obs(n)%gridDesco(5) = LIS_rc%obs_gridDesc(k,5)
      LDTSI_obs(n)%gridDesco(6) = LIS_rc%obs_gridDesc(k,6)
      LDTSI_obs(n)%gridDesco(7) = LIS_rc%obs_gridDesc(k,7)
      LDTSI_obs(n)%gridDesco(8) = LIS_rc%obs_gridDesc(k,8)
      LDTSI_obs(n)%gridDesco(9) = LIS_rc%obs_gridDesc(k,9)
      LDTSI_obs(n)%gridDesco(10) = LIS_rc%obs_gridDesc(k,10)
      LDTSI_obs(n)%gridDesco(20) = 255

      allocate(LDTSI_obs(n)%rlat1(LDTSI_obs(n)%mo1))
      allocate(LDTSI_obs(n)%rlon1(LDTSI_obs(n)%mo1))
      allocate(LDTSI_obs(n)%n111(LDTSI_obs(n)%mo1))

      call neighbor_interp_input_withgrid(tmp_gridDesci, &
           LDTSI_obs(n)%gridDesco, &
           LDTSI_obs(n)%mo1, LDTSI_obs(n)%rlat1, &
           LDTSI_obs(n)%rlon1, LDTSI_obs(n)%n111)

      ! Normal exit
      ierr = 0
   end subroutine LDTSI_reader

#else

   ! Dummy version
   subroutine LDTSI_reader(filename, n, k, snoanl, ierr)
      use LIS_logMod, only: LIS_logunit, LIS_endrun
      implicit none
      character(255), intent(in) :: filename
      integer, intent(in) :: n
      integer, intent(in) :: k
      real, allocatable, intent(out) :: snoanl(:,:)
      integer, intent(out) :: ierr
      ierr = 1
      write(LIS_logunit,*) &
           "[ERR] LIS not compiled with netCDF support!"
      write(LIS_logunit,*) &
           "[ERR] Recompile and try again!"
      call LIS_endrun()
   end subroutine LDTSI_reader

#endif

   ! Interpolates the LDTSI data to the LIS grid
   subroutine interp_LDTSIfield(n, k, snoanl, udef, varfield)

      ! Imports
      use LDTSIobs_Mod, only: LDTSI_obs
      use LIS_coreMod, only: LIS_rc

      ! Defaults
      implicit none

      ! Arguments
      integer, intent(in) :: n
      integer, intent(in) :: k
      real, intent(in) :: snoanl(LDTSI_obs(n)%nc_ldt, LDTSI_obs(n)%nr_ldt)
      real, intent(in) :: udef
      real, intent(out) :: varfield(LIS_rc%obs_lnc(k), LIS_rc%obs_lnr(k))

      ! Local variables
      integer :: c
      real, allocatable :: gi(:)
      real, allocatable :: go(:)
      integer :: i
      logical*1, allocatable :: li(:)
      logical*1, allocatable :: lo(:)
      integer :: r
      integer :: rc

      allocate(gi(LDTSI_obs(n)%nc_ldt*LDTSI_obs(n)%nr_ldt))
      gi(:) = 0
      allocate(li(LDTSI_obs(n)%nc_ldt*LDTSI_obs(n)%nr_ldt))
      li(:) = .false.

      do r = 1, LDTSI_obs(n)%nr_ldt
         do c = 1, LDTSI_obs(n)%nc_ldt
            if (snoanl(c,r) < 0) cycle
            li(c + (r-1)*LDTSI_obs(n)%nc_ldt) = .true.
            gi(c + (r-1)*LDTSI_obs(n)%nc_ldt) = snoanl(c,r)
         end do ! c
      end do ! r

      allocate(go(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
      go(:) = 0
      allocate(lo(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
      lo(:) = .false.

      call neighbor_interp(LDTSI_obs(n)%gridDesco, &
           li, gi, lo, go, &
           LDTSI_obs(n)%mi, LDTSI_obs(n)%mo1, &
           LDTSI_obs(n)%rlat1, &
           LDTSI_obs(n)%rlon1, &
           LDTSI_obs(n)%n111, &
           udef, rc)

      do r = 1, LIS_rc%obs_lnr(k)
         do c = 1, LIS_rc%obs_lnc(k)
            varfield(c,r) = go(c + (r-1)*LIS_rc%obs_lnc(k))
         end do ! c
      end do ! r

      deallocate(li)
      deallocate(lo)
      deallocate(go)
      deallocate(gi)
   end subroutine interp_LDTSIfield

end subroutine read_LDTSIobs
