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
! MODULE: USAFSI_lisMod
! 
! REVISION HISTORY:
! 01 Mar 2019  Eric Kemp  First version.
! 09 May 2019  Eric Kemp  Renamed to LDTSI.
! 13 Dec 2019  Eric Kemp  Renamed to USAFSI.
!
! DESCRIPTION:
! Source code for reading LIS 2-meter temperatures.
!-------------------------------------------------------------------------

#include "LDT_misc.h"

module USAFSI_lisMod

   ! Defaults
   implicit none
   private

   ! Public routines
   public :: read_gr2_t2

contains

#if (defined USE_GRIBAPI)
   ! Routine for reading LIS 2-m temperature
   subroutine read_gr2_t2(date10, nc, nr, sfctmp, ierr)

      ! Imports
      use grib_api
      use LDT_coreMod, only: LDT_rc, LDT_domain
      use LDT_logMod, only: LDT_logunit, LDT_verify, LDT_endrun
      use USAFSI_paramsMod

      ! Defaults
      implicit none

      ! Arguments
      character*10, intent(in) :: date10
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(out) :: sfctmp(nc,nr)
      integer, intent(out) :: ierr

      ! Local variables
      character*255 :: infilename
      logical :: file_exists
      character*12 :: routine_name
      integer :: ifile, igrib
      integer :: rc
      logical :: found_t2
      integer :: editionNumber
      integer :: discipline, parameterCategory, parameterNumber
      integer :: typeOfFirstFixedSurface, level
      integer :: ni, nj
      integer :: gridDefinitionTemplateNumber
      integer :: latitudeOfFirstGridPoint
      integer :: longitudeOfFirstGridPoint
      real :: swlat, swlon      
      integer :: latitudeOfLastGridPoint
      integer :: longitudeOfLastGridPoint
      real :: nelat, nelon      
      integer :: iDirectionIncrement
      integer :: jDirectionIncrement
      real :: dlon, dlat
      real :: gridDesci(20)
      integer, allocatable :: n11(:)
      real, allocatable :: tmp1d_in(:), tmp1d_out(:)
      logical*1, allocatable :: li(:), lo(:)
      integer :: c, r, i, j, gindex

      routine_name = 'read_gr2_t2'

      ierr = 1 ! Update below if we succeed in interpolating temperature

      ! Get the file name
      call construct_lis_grib2_filename(date10, infilename)

      ! See if the LIS GRIB2 file exists
      inquire(file=trim(infilename), exist=file_exists)
      if (.not. file_exists) then
         write(LDT_logunit,*) &
              '[WARN] Cannot find ' // trim(infilename)
         return
      end if

      ! Open the GRIB2 file
      call grib_open_file(ifile, trim(infilename), "r", rc)
      if (rc .ne. 0) then
         write(LDT_logunit,*) &
              "[WARN] Failed to open " // trim(infilename)
         return
      else
         write(LDT_logunit,*) &
              "[INFO] Reading " // trim(infilename)
      end if

      found_t2 = .false.

      do
         call grib_new_from_file(ifile, igrib, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read " // trim(infilename)
            exit
         end if
         
         ! Make sure this is a GRIB2 message
         call grib_get(igrib, "editionNumber", editionNumber, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read editionNumber from " // &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         if (editionNumber .ne. 2) then
            write(LDT_logunit,*) &
                 "[WARN] Found a non-GRIB2 message in " // trim(infilename)
            cycle ! We'll try reading the next GRIB message
         end if

         ! Check the discipline.  We need a meteorological product
         call grib_get(igrib, "discipline", discipline, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read discipline from " // trim(infilename)
            exit ! Assume corrupted file
         end if
         if (discipline .ne. 0) then
            cycle ! We'll try reading the next GRIB message
         end if

         ! Check the category.  We want temperature.
         call grib_get(igrib, "parameterCategory", parameterCategory, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read parameterCategory from "// &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         if (parameterCategory .ne. 0) then
            cycle ! We'll try reading the next GRIB message
         end if

         ! Check the number.  We want temperature
         call grib_get(igrib, "parameterNumber", parameterNumber, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read parameterNumber from "// &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         if (parameterNumber .ne. 0) then
            cycle ! We'll try reading the next GRIB message
         end if

         ! Check the first fixed surface type.  Should be specified height
         ! above ground.
         call grib_get(igrib, "typeOfFirstFixedSurface", &
              typeOfFirstFixedSurface, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read typeOfFirstFixedSurface from "// &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         if (typeOfFirstFixedSurface .ne. 103) then
            cycle ! We'll try reading the next GRIB message
         end if

         ! Check the level.  Should be 2-meters.
         call grib_get(igrib, "level", level, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read level from " // trim(infilename)
            exit ! Assume corrupted file
         end if
         if (level .ne. 2) then
            cycle ! We'll try reading the next GRIB message
         end if

         ! Check the map projection.  
         call grib_get(igrib, "gridDefinitionTemplateNumber", &
              gridDefinitionTemplateNumber, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read gridDefinitionTemplateNumber from " &
                 // trim(infilename)
            exit ! Assume corrupted file
         end if
         if (gridDefinitionTemplateNumber .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Unsupported grid type found in " // trim(infilename)
            write(LDT_logunit,*) &
                 "[WARN] Expected 0 (Latitude/longitude)"
            write(LDT_logunit,*) &
                 "[WARN] Found ", gridDefinitionTemplateNumber
            cycle ! We'll try reading the next GRIB message
         end if

         ! Get the dimensions
         call grib_get(igrib, "Ni", Ni, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read Ni from " // trim(infilename)
            exit ! Assume corrupted file
         end if
         call grib_get(igrib, "Nj", Nj, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read Nj from " // trim(infilename)
            exit ! Assume corrupted file
         end if

         ! Get the southwest latitude
         call grib_get(igrib, "latitudeOfFirstGridPoint", &
              latitudeOfFirstGridPoint, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read latitudeOfFirstGridPoint from " // &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         swlat = latitudeOfFirstGridPoint*(1.e-6)

         ! Get the southwest longitude
         call grib_get(igrib, "longitudeOfFirstGridPoint", &
              longitudeOfFirstGridPoint, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read longitudeOfFirstGridPoint from " // &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         swlon = longitudeOfFirstGridPoint*(1.e-6)
         ! LIS wants longitude range from -180 to 180, not 0 to 360
         if (swlon > 180) then
            swlon = swlon - 360.
         end if

         ! Get the northeast latitude
         call grib_get(igrib, "latitudeOfLastGridPoint", &
              latitudeOfLastGridPoint, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read latitudeOfLastGridPoint from " // &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         nelat = latitudeOfLastGridPoint*(1.e-6)

         ! Get the northeast longitude
         call grib_get(igrib, "longitudeOfLastGridPoint", &
              longitudeOfLastGridPoint, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read longitudeOfLastGridPoint from " // &
                 trim(infilename)
            exit ! Assume corrupted file
         end if
         nelon = longitudeOfLastGridPoint*(1.e-6)

         ! Get the longitude resolution
         call grib_get(igrib, "iDirectionIncrement", &
              iDirectionIncrement, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read iDirectionIncrementInDegrees from " &
                 // trim(infilename)
            exit ! Assume corrupted file
         end if
         dlon = iDirectionIncrement*(1.e-6)

         ! Get the latitude resolution
         call grib_get(igrib, "jDirectionIncrement", &
              jDirectionIncrement, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read jDirectionIncrement from " &
                 // trim(infilename)
            exit ! Assume corrupted file
         end if
         dlat = jDirectionIncrement*(1.e-6)

         ! Get the temperature field
         allocate(tmp1d_in(ni*nj))
         call grib_get(igrib, "values", tmp1d_in, rc)
         if (rc .ne. 0) then
            write(LDT_logunit,*) &
                 "[WARN] Failed to read values from " &
                 // trim(infilename)
         end if

         ! Set up grid information for interpolator
         gridDesci = 0
         gridDesci(1) = 0 
         gridDesci(2) = ni
         gridDesci(3) = nj
         gridDesci(4) =  swlat
         gridDesci(5) =  swlon
         gridDesci(7) =  nelat
         gridDesci(8) =  nelon
         gridDesci(6) = 128
         gridDesci(9) =    dlon
         gridDesci(10) =   dlat
         gridDesci(11) = 64
         gridDesci(20) = 64

         allocate(n11(ni*nj))
         call neighbor_interp_input(1, gridDesci, n11)
         
         allocate(li(ni*nj))
         li(:) = .false.
         do j = 1, nj
            do i = 1, ni
               gindex = i + (j-1)*ni
               if (tmp1d_in(gindex) < 9998) then
                  li(gindex) = .true.
               end if
            end do ! i
         end do ! j
         
         ! Now interpolate
         allocate(tmp1d_out(nc*nr))
         allocate(lo(nc*nr))
         call neighbor_interp(LDT_rc%gridDesc(1,:), li, tmp1d_in, &
              lo, tmp1d_out, ni*nj, nc*nr, &
              LDT_domain(1)%lat, LDT_domain(1)%lon, &
              n11, 9999., rc)
          
         ! Copy the interpolated temperature field
         sfctmp(:,:) = -1
         do r = 1, nr
            do c = 1, nc
               gindex = c + (r-1)*nc
               if (tmp1d_out(gindex) < 9998) then
                  sfctmp(c,r) = tmp1d_out(gindex)
               end if
            end do ! c
         end do ! r

         ! Clean up
         if (allocated(tmp1d_in)) deallocate(tmp1d_in)
         if (allocated(tmp1d_out)) deallocate(tmp1d_out)
         if (allocated(lo)) deallocate(lo)
         if (allocated(li)) deallocate(li)
         if (allocated(n11)) deallocate(n11)
         
         ! We found the data.  Break out of loop
         found_t2 = .true.
         exit
      end do

      ! Close the GRIB2 file
      call grib_close_file(ifile)

      ! Give up if we didn't find the 2-meter temperature
      if (.not. found_t2) then
         write(LDT_logunit,*) &
              "[WARN] Did not find 2-meter temperature from file" // &
              trim(infilename)
         return
      end if
      
      ! The end
      ierr = 0

   end subroutine read_gr2_t2

#else
   ! Dummy version
   subroutine read_gr2_t2(date10, nc, nr, sfctmp, ierr)
      use LDT_logMod, only: LDT_logunit, LDT_endrun
      implicit none
      character*10, intent(in) :: date10
      integer, intent(in) :: nc
      integer, intent(in) :: nr
      real, intent(out) :: sfctmp(nc,nr)
      integer, intent(out) :: ierr
      write(LDT_logunit,*) &
           '[ERR] LDT was compiled without GRIB_API or ECCODES support!'
      write(LDT_logunit,*) &
           '[ERR] Recompile and try again!'
      ierr = 1
      call LDT_endrun()
   end subroutine read_gr2_t2
#endif

   ! Builds path to LIS GRIB2 file
   subroutine construct_lis_grib2_filename(date10, filename)
      use LDT_usafsiMod, only: usafsi_settings
      implicit none
      character*10, intent(in) :: date10
      character*255, intent(out) :: filename
      filename = trim(usafsi_settings%lis_grib2_dir) &
           // "/PS." &
           // "557WW_SC." &
           // trim(usafsi_settings%security_class) // "_DI." &
           // trim(usafsi_settings%data_category) // "_GP." &
           // "LIS_GR." &
           // trim(usafsi_settings%data_res) // "_AR." &
           // trim(usafsi_settings%area_of_data) // "_PA." &
           // "LIS_DD." &
           // date10(1:8) // "_DT." &
           // date10(9:10) // "00_DF." &
           // "GR2"                     
   end subroutine construct_lis_grib2_filename

end module USAFSI_lisMod
