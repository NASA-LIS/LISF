!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: LDTSI_gofsMod
! 
! REVISION HISTORY:
! 01 Apr 2019  Eric Kemp  First version.
! 09 May 2019  Eric Kemp  Rename to LDTSI
!
! DESCRIPTION:
! Source code for reading US Navy GOFS data.
!-------------------------------------------------------------------------

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDTSI_gofsMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: process_gofs_sst
   public :: process_gofs_cice
   
contains

   ! Find GOFS CICE file on file system
   subroutine find_gofs_cice_file(rootdir, region, yyyy, mm, dd, hh, fh, &
        filename)

      ! Imports
      use LDT_logMod, only: LDT_logunit
      use LDT_timeMgrMod, only: LDT_get_julhr, LDT_julhr_date
      
      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: rootdir
      character*3, intent(in) :: region
      integer, intent(inout) :: yyyy
      integer, intent(inout) :: mm
      integer, intent(inout) :: dd
      integer, intent(inout) :: hh
      integer, intent(inout) :: fh
      character*255, intent(out) :: filename

      ! Local variables
      integer :: julhr, julhr_orig
      logical :: file_exists

      ! Build the file name.  Note that all GOFS CICE runs start at 12Z.
      call LDT_get_julhr(yyyy, mm, dd, 12, 0, 0, julhr)
      if (hh >= 12) then
         julhr_orig = julhr
         fh = 0
      else
         julhr_orig = julhr - 24 ! Must use previous day's run
         fh = 12
      end if
      call LDT_julhr_date(julhr_orig, yyyy, mm, dd, hh)
      call construct_gofs_cice_filename(rootdir, region, &
           yyyy, mm, dd, hh, fh, filename)
      
      write(LDT_logunit,*) &
           '------------------------------------------------------------------'
      write(LDT_logunit,*)'[INFO] *** SEARCHING FOR GOFS CICE FOR ',&
           trim(region),' REGION ***'
      inquire(file=trim(filename),exist=file_exists)
      if (file_exists) then
         write(LDT_logunit,*)'[INFO] Will use ',trim(filename)
         return        
      end if

      ! At this point, we are rolling back to earlier CICE file
      ! Start looping for earlier files
      julhr = julhr_orig
      do
         write(LDT_logunit,*)'[WARN] Cannot find ',trim(filename)
         fh = fh + 24
         julhr = julhr - 24
         if ( (julhr_orig - julhr) > 24*5) then
            write(LDT_logunit,*)&
                '[WARN] *** GIVING UP ON GOFS CICE FOR ',trim(region),' ***'
            write(LDT_logunit,*) &
                '[WARN] *** NO GOFS CICE DATA FOR ',trim(region), &
                ' AVAILABLE!!! ***'
            filename = 'NONE'
            return
         end if
         call LDT_julhr_date(julhr, yyyy, mm, dd, hh)
         
         call construct_gofs_cice_filename(rootdir, region, &
              yyyy, mm, dd, hh, fh, filename)
         inquire(file=trim(filename),exist=file_exists)
         if (file_exists) then
            write(LDT_logunit,*)'[INFO] Will use ',trim(filename)
            return        
         end if
      end do
      
   end subroutine find_gofs_cice_file

   ! Builds path to GOFS CICE netcdf file
   subroutine construct_gofs_cice_filename(rootdir, region, &
           yyyy, mm, dd, hh, fh, filename)

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: rootdir
      character*3, intent(in) :: region
      integer, intent(in) :: yyyy
      integer, intent(in) :: mm
      integer, intent(in) :: dd
      integer, intent(in) :: hh
      integer, intent(in) :: fh
      character*255, intent(out) :: filename

      ! Local variables
      character*10 :: yyyymmddhh
      character*4  :: thhh

      write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yyyy, mm, dd, hh
      write(thhh,'(a1,i3.3)') 't', fh

      filename = trim(rootdir) // '/hycom-cice_inst_' // trim(region) &
           // 'u0.08_930_' // yyyymmddhh // '_' // thhh // '.nc'

   end subroutine construct_gofs_cice_filename

   ! Find GOFS SST file on file system
   subroutine find_gofs_sst_file(rootdir, yyyy, mm, dd, hh, fh, &
        filename)

      ! Imports
      use LDT_logMod, only: LDT_logunit
      use LDT_timeMgrMod, only: LDT_get_julhr, LDT_julhr_date

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: rootdir
      integer, intent(inout) :: yyyy
      integer, intent(inout) :: mm
      integer, intent(inout) :: dd
      integer, intent(inout) :: hh
      integer, intent(inout) :: fh
      character*255, intent(inout) :: filename

      ! Local variables
      integer :: julhr, julhr_orig
      logical :: file_exists
            

      ! Build the file name.  Note that all GOFS SST runs start at 00Z.
      if (hh < 6) then
         fh = 0
      else if (hh < 12) then
         fh = 6
      else if (hh < 18) then
         fh = 12
      else
         fh = 18
      end if
      hh = 0 
      call construct_gofs_sst_filename(rootdir, &
           yyyy, mm, dd, hh, fh, filename)
      
      ! Check if file exists
      write(LDT_logunit,*) &
           '------------------------------------------------------------------'
      write(LDT_logunit,*) &
           '[INFO] *** SEARCHING FOR GOFS SST ***'
      inquire(file=trim(filename),exist=file_exists)
      if (file_exists) then
         write(LDT_logunit,*)'[INFO] Will use ',trim(filename)
         return
      end if

      ! At this point, we are rolling back to earlier SST file
      call LDT_get_julhr(yyyy, mm, dd, hh, 0, 0, julhr)
      julhr_orig = julhr

      ! Start looping for earlier files
      do
         write(LDT_logunit,*)'[WARN] Cannot find ',trim(filename)
         fh = fh - 6
         if (fh < 0) then
            fh = 24
            julhr = julhr - 24 ! Roll back to previous 00Z cycle
            ! Give up after 5 days
            if ( (julhr_orig - julhr) > 24*5) then
               write(LDT_logunit,*)"[WARN] *** GIVING UP ON GOFS SST! ***"
               write(LDT_logunit,*)"[WARN] *** NO GOFS SST AVAILABLE!!! ***"
               filename = "NONE"
               return
            end if
            call LDT_julhr_date(julhr, yyyy, mm, dd, hh)
         end if

         call construct_gofs_sst_filename(rootdir, &
              yyyy, mm, dd, hh, fh, filename)
         inquire(file=trim(filename),exist=file_exists)
         if (file_exists) then
            write(LDT_logunit,*)'[INFO] Will use ',trim(filename)
            return
         end if
      end do

   end subroutine find_gofs_sst_file

   ! Builds path to GOFS SST netcdf file
   subroutine construct_gofs_sst_filename(rootdir, &
        yyyy, mm, dd, hh, fh, filename)

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: rootdir
      integer, intent(in) :: yyyy
      integer, intent(in) :: mm
      integer, intent(in) :: dd
      integer, intent(in) :: hh
      integer, intent(in) :: fh
      character*255, intent(out) :: filename

      ! Local variables
      character*10 :: yyyymmddhh
      character*4  :: thhh

      write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') yyyy,mm,dd,hh
      write(thhh,'(a1,i3.3)') 't',fh

      filename = trim(rootdir) // "/hycom_glb_sfc_u_" // yyyymmddhh // &
           "_" // thhh // ".nc"

   end subroutine construct_gofs_sst_filename


#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   ! Read GOFS sea surface temperature and reproject to LDT grid
   subroutine process_gofs_sst(rootdir, nc, nr, landmask, sst, &
        yyyy, mm, dd, hh, fh, &
        ierr)

      ! Imports
      use LDT_coreMod, only: LDT_rc, LDT_domain
      use LDT_logMod, only: LDT_verify
      use netcdf

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: rootdir
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)
      real, intent(inout) :: sst(nc,nr)      
      integer, intent(inout) :: yyyy
      integer, intent(inout) :: mm
      integer, intent(inout) :: dd
      integer, intent(inout) :: hh
      integer, intent(inout) :: fh
      integer, intent(out) :: ierr

      ! Local parameters
      integer, parameter :: nlat = 2001
      integer, parameter :: nlon = 4500

      ! Local variables
      character*255 :: filename
      integer :: ncid, water_temp_varid
      real, allocatable :: water_temp(:,:,:,:)
      real, allocatable :: water_temp_1d(:)
      real, allocatable :: sst_1d(:)
      integer :: c, r, c1, r1
      logical*1, allocatable :: lb(:)
      logical*1, allocatable :: lo(:)
      real :: griddesci(50)
      real, allocatable :: n11(:)
      integer :: gindex
      real :: rlat

      ! Find a valid file on the file system
      call find_gofs_sst_file(rootdir, yyyy, mm, dd, hh, fh, filename)
      if (trim(filename) == "NONE") then
         ierr = 1
         return
      end if

      ! Open the file
      call LDT_verify(nf90_open(path=trim(filename), &
           mode=nf90_nowrite, &
           ncid=ncid), &
           "[ERR] Error in nf90_open for " // trim(filename))

      ! Get the varid for water_temp
      call LDT_verify(nf90_inq_varid(ncid, "water_temp", water_temp_varid), &
           "[ERR] Error in nf90_inq_varid for water_temp")

      ! Allocate the water_temp array
      allocate(water_temp(nlon, nlat, 1, 1))

      ! Pull from the GOFS file
      call LDT_verify(nf90_get_var(ncid, water_temp_varid, water_temp), &
           "[ERR] Error in nf90_get_var for water_temp")

      ! Close the file
      call LDT_verify(nf90_close(ncid), &
           "[ERR] Error in nf90_close for "// trim(filename))
      
      ! We need to interpolate to the LDT grid.  First, copy to 1D array
      allocate(water_temp_1d(nlon*nlat*1*1))
      water_temp_1d(:) = -9999.0
      allocate(lb(nlon*nlat*1*1))
      lb(:) = .false.
      do r = 1, nlat
         do c = 1, nlon
            if (water_temp(c,r,1,1) .eq. -30000) cycle
            ! Convert from Celsius to Kelvin, taking into account the scale
            ! factor and offset.  Also, wrap the data so it starts at 180W
            if (c .gt. 2250) then
               c1 = c - 2250
               r1 = r
            else
               c1 = c + 2250
               r1 = r
            end if
            water_temp_1d(c1 + (r1-1)*nlon) = &
                 (0.001*water_temp(c,r,1,1)) + 20.0 + 273.15
            lb(c1 + (r1-1)*nlon) = .true.
         end do ! c
      end do ! r
      deallocate(water_temp)
      
      ! Set up interpolation weights
      gridDesci = 0
      gridDesci(1) = 0 
      gridDesci(2) = nlon
      gridDesci(3) = nlat
      gridDesci(4) =  -80.0
      gridDesci(5) = -180.0
      gridDesci(7) =   80.0
      gridDesci(8) =  180.0
      gridDesci(6) = 128
      gridDesci(9) =    0.08
      gridDesci(10) =   0.08
      gridDesci(20) = 64
      allocate(n11(nlon*nlat))
      call upscaleByAveraging_input(gridDesci, LDT_rc%gridDesc, &
           nlon*nlat, nc*nr, n11)
      
      ! Now interpolate
      allocate(sst_1d(nc*nr))
      sst_1d = -9999.
      allocate(lo(nc*nr))
      lo(:) = .false.
      call upscaleByAveraging(nlon*nlat, nc*nr, -9999., &
           n11, lb, water_temp_1d, lo, sst_1d)

      ! Since SST is missing north of 80N, we need to set water points in
      ! this region to a reasonable value.  We follow the typical 
      ! UKMET SURF value of 271.35K.
      sst(:,:) = -1
      do r = 1, nr
         do c = 1, nc
            ! Skip land points
            if (landmask(c,r) >= 0.5) cycle

            gindex = c + (r-1)*nc
            rlat = LDT_domain(1)%lat(gindex)
            if (rlat >= 80.) then
               sst(c,r) = 271.35
            else
               if (sst_1d(gindex) > 0) then
                  sst(c,r) = sst_1d(gindex)
               end if
            end if
         end do ! c
      end do ! r

      ! Clean up
      deallocate(water_temp_1d)
      deallocate(lb)
      deallocate(lo)
      deallocate(sst_1d)
      deallocate(n11)

      ! The end
      ierr = 0
   end subroutine process_gofs_sst

#else

   ! Dummy version with no netCDF support
   subroutine process_gofs_sst(rootdir, nc, nr, landmask, sst, &
        yyyy, mm, dd, hh, fh, &
        ierr)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: rootdir
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)
      real, intent(inout) :: sst(nc,nr)      
      integer, intent(inout) :: yyyy
      integer, intent(inout) :: mm
      integer, intent(inout) :: dd
      integer, intent(inout) :: hh
      integer, intent(inout) :: fh
      integer, intent(out) :: ierr

      write(LDT_logunit,*) &
           '[ERR] LDT was compiled without netCDF support!'
      write(LDT_logunit,*) "[ERR] Recompile and try again!"
      ierr = 1
      call LDT_endrun()
   end subroutine process_gofs_sst

#endif

   ! Read GOFS sea ice and reproject to LDT grid
   subroutine process_gofs_cice(rootdir, nc, nr, landmask, icecon, &
        yyyy, mm, dd, hh, fh, &
        ierr)
      
      ! Imports
      use LDT_coreMod, only: LDT_domain

      ! Defaults
      implicit none

      ! Arguments
      character(len=*), intent(in) :: rootdir
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)
      real, intent(inout) :: icecon(nc,nr)      
      integer, intent(inout) :: yyyy
      integer, intent(inout) :: mm
      integer, intent(inout) :: dd
      integer, intent(inout) :: hh
      integer, intent(inout) :: fh
      integer, intent(out) :: ierr

      ! Local variables
      real, allocatable :: icecon_arc(:,:)
      real, allocatable :: icecon_ant(:,:)
      integer :: c, r
      integer :: gindex
      real :: rlat

      ! First handle Arctic region
      call process_gofs_cice_region('ARC', rootdir, nc, nr, landmask, &
           yyyy, mm, dd, hh, fh, icecon_arc, ierr)
      if (ierr .ne. 0) then
         if (allocated(icecon_arc)) deallocate(icecon_arc)
         return
      end if

      ! Next handle Antarctic region
      call process_gofs_cice_region('ANT', rootdir, nc, nr, landmask, &
           yyyy, mm, dd, hh, fh, icecon_ant, ierr)
      if (ierr .ne. 0) then
         if (allocated(icecon_arc)) deallocate(icecon_arc)
         if (allocated(icecon_ant)) deallocate(icecon_ant)
         return
      end if

      ! Merge the two regions together
      icecon(:,:) = -1
       do r = 1, nr
          do c = 1, nc
             ! Skip land points
             if (landmask(c,r) > 0.5) cycle

             gindex = c + (r-1)*nc
             rlat = LDT_domain(1)%lat(gindex)
             if (rlat >= 0) then
                icecon(c,r) = icecon_arc(c,r)
             else
                icecon(c,r) = icecon_ant(c,r)
             end if

          end do ! c
       end do ! r

      ! Clean up
      deallocate(icecon_arc)
      deallocate(icecon_ant)
      ierr = 0
   end subroutine process_gofs_cice

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

   ! Process a particular region of GOFS CICE data (Arctic or Antarctic
   subroutine process_gofs_cice_region(region, rootdir, nc, nr, landmask, &
        yyyy, mm, dd, hh, fh, icecon, ierr)

      ! Imports
      use LDT_coreMod, only: LDT_rc
      use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify
      use netcdf

      ! Arguments
      character*3, intent(in) :: region
      character(len=*), intent(in) :: rootdir
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)
      integer, intent(inout) :: yyyy
      integer, intent(inout) :: mm
      integer, intent(inout) :: dd
      integer, intent(inout) :: hh
      integer, intent(inout) :: fh
      real, allocatable, intent(out) :: icecon(:,:)      
      integer, intent(out) :: ierr

      ! Local parameters
      integer, parameter :: nlat_arc = 1251
      integer, parameter :: nlat_ant =  775
      integer, parameter :: nlon = 4500

      ! Local variables
      character*255 :: filename
      integer :: ncid, aice_varid
      real, allocatable :: aice(:,:,:)
      real, allocatable :: aice_1d(:)
      real, allocatable :: icecon_1d(:)
      integer :: c, r
      logical*1, allocatable :: lb(:)
      logical*1, allocatable :: lo(:)
      real :: griddesci(50)
      real, allocatable :: n11(:)
      integer :: gindex, nlat

      ! Sanity check the region
      if (region .eq. 'ARC') then
         nlat = nlat_arc
      else if (region .eq. 'ANT') then
         nlat = nlat_ant
      else
         write(LDT_logunit,*)'[ERR] Invalid GOFS region for cice: ' // region
         write(LDT_logunit,*)'[ERR] Must be either ARC or ANT'
         ierr = 1
         call LDT_endrun()
      end if

      ! Find a valid file on the file system
      call find_gofs_cice_file(rootdir, region, yyyy, mm, dd, hh, fh, filename)
      if (trim(filename) == "NONE") then
         ierr = 1
         return
      end if

      ! Open the file
      call LDT_verify(nf90_open(path=trim(filename), &
           mode=nf90_nowrite, &
           ncid=ncid), &
           "[ERR] Error in nf90_open for " // trim(filename))

      ! Get the varid for aice
      call LDT_verify(nf90_inq_varid(ncid, "aice", aice_varid), &
           "[ERR] Error in nf90_inq_varid for aice")

      ! Allocate the aice array
      allocate(aice(nlon, nlat, 1))

      ! Pull from the GOFS file
      call LDT_verify(nf90_get_var(ncid, aice_varid, aice), &
           "[ERR] Error in nf90_get_var for aice")

      ! Close the file
      call LDT_verify(nf90_close(ncid), &
           "[ERR] Error in nf90_close for "// trim(filename))

      ! We need to interpolate to the LDT grid.  First, copy to 1D array
      allocate(aice_1d(nlon*nlat*1))
      aice_1d(:) = -9999
      allocate(lb(nlon*nlat*1))
      lb(:) = .false.
      do r = 1, nlat
         do c = 1, nlon
            if (aice(c,r,1) .eq. -30000) cycle
            
            ! Take into account the scale factor and offset
            aice_1d(c + (r-1)*nlon) = &
                 aice(c,r,1)*0.0001*100
            lb(c + (r-1)*nlon) = .true.
         end do ! c
      end do ! r
      deallocate(aice)

      ! Set up interpolation weights
      if (region .eq. 'ARC') then
         gridDesci = 0 
         gridDesci(1) = 0 ! Lat/lon projection
         gridDesci(2) = nlon
         gridDesci(3) = nlat
         gridDesci(4) =   40.     ! Lower-left latitude (deg N)
         gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
         gridDesci(6) = 128     ! Not used
         gridDesci(7) =   90.0             ! Upper-right latitude (deg N)
         gridDesci(8) =  179.920043945312 ! Upper-right longitude (deg E)
         gridDesci(9) =    0.080017089844005795  ! delta-lon (deg)
         gridDesci(10) =   0.040000915527301117 ! delta-lat (deg)
         gridDesci(20) = 64  ! East-west ordering         
      else if (region .eq. 'ANT') then
         gridDesci = 0 
         gridDesci(1) = 0 ! Lat/lon projection
         gridDesci(2) = nlon
         gridDesci(3) = nlat
         gridDesci(4) =  -80.4800033569336     ! Lower-left latitude (deg N)
         gridDesci(5) = -180.0  ! Lower-left longitude (deg E)
         gridDesci(6) = 128     ! Not used
         gridDesci(7) =  -49.5200004577637 ! Upper-right latitude (deg N)
         gridDesci(8) =  179.920043945312 ! Upper-right longitude (deg E)
         gridDesci(9) =    0.080017089844005795  ! delta-lon (deg)
         gridDesci(10) =   0.040000915527400593 ! delta-lat (deg)
         gridDesci(20) = 64  ! East-west ordering
      end if
      allocate(n11(nlon*nlat))

      call upscaleByAveraging_input(gridDesci, LDT_rc%gridDesc, &
           nlon*nlat, nc*nr, n11)

      ! Now interpolate
      allocate(icecon_1d(nc*nr))
      icecon_1d(:) = -9999
      allocate(lo(nc*nr))
      lo(:) = .false.
      call upscaleByAveraging(nlon*nlat, nc*nr, -9999., &
           n11, lb, aice_1d, lo, icecon_1d)

      ! Just copy the non-missing values to the output array.  This should
      ! prevent overwriting of data outside of the GOFS polar region.
      allocate(icecon(nc,nr))
      do r = 1, nr
         do c = 1, nc
            ! Skip land points
            if (landmask(c,r) >= 0.5) cycle

            gindex = c + (r-1)*nc
            if (icecon_1d(gindex) .ne. -9999) then
               icecon(c,r) = icecon_1d(gindex)
            end if
         end do ! c
      end do ! r

      ! Clean up
      deallocate(aice_1d)
      deallocate(lb)
      deallocate(lo)
      deallocate(icecon_1d)
      deallocate(n11)

      ! The end
      ierr = 0
   end subroutine process_gofs_cice_region

#else
   ! Dummy version
   subroutine process_gofs_cice_region(region, rootdir, nc, nr, landmask, &
        icecon, yyyy, mm, dd, hh, fh, ierr)

      ! Imports                                                                
      use LDT_logMod, only: LDT_logunit, LDT_endrun

      ! Arguments
      character*3, intent(in) :: region
      character(len=*), intent(in) :: rootdir
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(in) :: landmask(nc,nr)
      real, intent(inout) :: icecon(nc,nr)      
      integer, intent(inout) :: yyyy
      integer, intent(inout) :: mm
      integer, intent(inout) :: dd
      integer, intent(inout) :: hh
      integer, intent(inout) :: fh
      integer, intent(out) :: ierr

      write(LDT_logunit,*) &
           '[ERR] LDT was compiled without netCDF support!'
      write(LDT_logunit,*) "[ERR] Recompile and try again!"
      ierr = 1
      call LDT_endrun()

   end subroutine process_gofs_cice_region

#endif

end module LDTSI_gofsMod
