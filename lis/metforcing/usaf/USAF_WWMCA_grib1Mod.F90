!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"

!
! MODULE: USAF_WWMCA_grib1Mod
!
! DESCRIPTION:
! This module contains code for reading WWMCA cloud amounts, types, tops,
! and time information from Air Force GRIB1 files, with data projected
! on 16th mesh polar stereographic grids.
!
! REVISION HISTORY:
! 23 Nov 2020: Eric Kemp, Initial specification.
!

module USAF_WWMCA_grib1Mod

  ! Defaults
  implicit none
  private

  ! Create structure to share internal routines
  type :: wwmca_grib1_t
     private
     character(len=255) :: full_grib_path
     integer :: nx
     integer :: ny
     integer :: nk
     real, allocatable :: cldamt(:,:,:)
     real, allocatable :: cldtyp(:,:,:)
     real, allocatable :: cldtop(:,:,:)
     real, allocatable :: pixtim(:,:)
   contains
     procedure :: new
     procedure :: destroy
     procedure :: read_grib
     procedure :: return_binary_fields
  end type wwmca_grib1_t

  type(wwmca_grib1_t), public :: wwmca_grib1

  public :: USAF_wwmca_grib1_filename

  ! Local constants for the 16th mesh WWMCA data
  integer, parameter :: NX = 1024
  integer, parameter :: NY = 1024
  integer, parameter :: NZ = 4

contains

  ! Constructor
  subroutine new(this, full_grib_path)

    ! Defaults
    implicit none

    ! Arguments
    class(wwmca_grib1_t), intent(inout) :: this
    character(*), intent(in) :: full_grib_path

    this%full_grib_path = trim(full_grib_path)
    this%nx = NX
    this%ny = NY
    this%nk = NZ

    allocate(this%cldamt(NX, NY, NZ))
    this%cldamt = -9999

    allocate(this%cldtyp(NX, NY, NZ))
    this%cldtyp = -9999

    allocate(this%cldtop(NX, NY, NZ))
    this%cldtop = -9999

    allocate(this%pixtim(NX, NY))
    this%pixtim = -9999

  end subroutine new

  ! Destructor
  subroutine destroy(this)
    implicit none
    class(wwmca_grib1_t), intent(inout) :: this
    this%full_grib_path = "NULL"
    if (allocated(this%cldamt)) deallocate(this%cldamt)
    if (allocated(this%cldtyp)) deallocate(this%cldtyp)
    if (allocated(this%cldtop)) deallocate(this%cldtop)
    if (allocated(this%pixtim)) deallocate(this%pixtim)
  end subroutine destroy

  ! Read the data from the GRIB file
  subroutine read_grib(this, year, month, day, hour, minute, rc)

    ! Imports
#if (defined USE_GRIBAPI)
    use grib_api
#endif
    use LIS_logMod, only: LIS_logunit

    ! Defaults
    implicit none

    ! Arguments
    class(wwmca_grib1_t), intent(inout) :: this
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    integer, intent(in) :: minute
    integer, intent(inout) :: rc

    ! Locals
    logical :: file_exists
    integer :: ftn
    integer :: nmsgs
    real, allocatable :: dum1d(:)
    integer :: counter
    integer :: k
    integer :: igrib
    integer :: edition
    integer :: grid
    integer :: inx
    integer :: iny
    integer :: firstlat
    integer :: firstlon
    integer :: orient
    integer :: dx
    integer :: dy
    integer :: res_component_flags
    integer :: datadate
    integer :: datatime
    integer :: param
    integer :: leveltype
    integer :: level

    ! See if file exists
    inquire(file=trim(this%full_grib_path), exist=file_exists)
    if (.not. file_exists) then
       write(LIS_logunit,*) '[WARN] ', trim(this%full_grib_path), &
            ' does not exist!'
       rc = 1
       return
    end if

#if (defined USE_GRIBAPI)

    ! Open the file through GRIB_API
    call grib_open_file(ftn, trim(this%full_grib_path), 'r', rc)
    if (rc .ne. 0) then
       write(LIS_logunit,*)'[WARN] Cannot open ', trim(this%full_grib_path)
       rc = 1
       return
    end if

    ! Check total number of messages in this file.
    call grib_count_in_file(ftn, nmsgs, rc)
    if (rc .ne. 0) then
       write(LIS_logunit,*)'[WARN] Cannot count messages in ', &
            trim(this%full_grib_path)
       call grib_close_file(ftn)
       rc = 1
       return
    end if

    allocate(dum1d(NX*NY))
    counter = 0

    write(LIS_logunit,*)'[INFO] Reading ', trim(this%full_grib_path)
    write(LIS_logunit,'(1X, A, I8.8, A, I4.4)') &
         '[INFO] Searching for data valid ', &
         (year*10000 + month*100 + day), '/', (hour*100 + minute)

    ! Loop through the GRIB file until all fields are found, or we reach
    ! end of file
    do k = 1, nmsgs

       ! Find next message
       call grib_new_from_file(ftn, igrib, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*)'[WARN] Cannot read from ', &
               trim(this%full_grib_path)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if

       ! Check edition number of this message
       call grib_get(igrib, 'editionNumber', edition, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*)'[WARN] Cannot read editionNumber from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (edition .ne. 1) then
          write(LIS_logunit,*)'[WARN] Did not find GRIB1 message in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the grid definition -- should be 16th mesh, either northern
       ! or southern hemisphere
       call grib_get(igrib, 'gridDefinition', grid, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*)'[WARN] Cannot read gridDefinition from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (grid .ne. 212 .and. grid .ne. 213) then
          write(LIS_logunit,*)'[WARN] Did not find 16th mesh data from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the nx dimension
       call grib_get(igrib, 'Nx', inx, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*)'[WARN] Cannot read Nx from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (inx .ne. NX) then
          write(LIS_logunit,*)'[WARN] Found wrong Nx dimension in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the ny dimension
       call grib_get(igrib, 'Ny', iny, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*)'[WARN] Cannot read Ny from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (iny .ne. NY) then
          write(LIS_logunit,*)'[WARN] Found wrong Ny dimension in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the first latitude
       call grib_get(igrib, 'latitudeOfFirstGridPoint', firstlat, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read latitudeOfFirstGridPoint from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (.not. (grid .eq. 212 .and. firstlat .eq. -20826) .and. &
            .not. (grid .eq. 213 .and. firstlat .eq. 20826)) then
          write(LIS_logunit,*) &
               '[WARN] Found wrong latitudeOfFirstGridPoint in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the first longitude
       call grib_get(igrib, 'longitudeOfFirstGridPoint', firstlon, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read longitudeOfFirstGridPoint from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (.not. (grid .eq. 212 .and. firstlon .eq. 145000) .and. &
            .not. (grid .eq. 213 .and. firstlon .eq. -125000)) then
          write(LIS_logunit,*) &
               '[WARN] Found wrong longitudeOfFirstGridPoint in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the grid orientation
       call grib_get(igrib, 'orientationOfTheGrid', orient, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read orientationOfTheGrid from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (orient .ne. 100000) then
          write(LIS_logunit,*) &
               '[WARN] Found wrong orientationOfTheGrid in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the grid resolution
       call grib_get(igrib, 'DxInMetres', dx, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read DxInMetres from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (dx .ne. 23813) then
          write(LIS_logunit,*) &
               '[WARN] Found wrong DxInMetres in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the grid resolution
       call grib_get(igrib, 'DyInMetres', dy, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read DyInMetres from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (dy .ne. 23813) then
          write(LIS_logunit,*) &
               '[WARN] Found wrong DyInMetres in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the resolution and component flags
       call grib_get(igrib, 'resolutionAndComponentFlags', &
            res_component_flags, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read resolutionAndComponentFlags from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if (res_component_flags .ne. 128) then
          write(LIS_logunit,*) &
               '[WARN] Found wrong resolutionAndComponentFlags in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the valid date
       call grib_get(igrib, 'dataDate', datadate, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read dataDate from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if ((year*10000 + month*100 + day) .ne. datadate) then
          write(LIS_logunit,*) &
               '[WARN] Found wrong datadate in ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          cycle
       end if

       ! Check the valid time
       call grib_get(igrib, 'dataTime', datatime, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read dataTime from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if
       if ((hour*100 + minute) .ne. datatime) then
          ! write(LIS_logunit,*) &
          !      '[WARN] Found wrong dataTime in ', &
          !      trim(this%full_grib_path)
          ! write(LIS_logunit,*) &
          !      '[WARN] dataTime = ', datatime
          ! write(LIS_logunit,*) &
          !      '[WARN] Wanted ', hour*100 + minute
          call grib_release(igrib, rc)
          cycle
       end if

       ! At this point we are ready to check the parameter and level
       ! information.
       call grib_get(igrib, 'indicatorOfParameter', param, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read indicatorOfParameter from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if

       call grib_get(igrib, 'indicatorOfTypeOfLevel', leveltype, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read indicatorOfTypeOfLevel from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if

       call grib_get(igrib, 'level', level, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*) &
               '[WARN] Cannot read level from ', &
               trim(this%full_grib_path)
          call grib_release(igrib, rc)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if

       ! See if this field is desired
       rc = 0
       if (param .eq. 163 .and. leveltype .eq. 109 .and. &
            level .gt. 0 .and. level .lt. 5) then
          ! Cloud amount (%)
          call fetch_values(this%full_grib_path, igrib, inx, iny, dum1d, &
               this%cldamt(:, :, level), counter, rc)
       else if (param .eq. 164 .and. leveltype .eq. 109 .and. &
            level .gt. 0 .and. level .lt. 5) then
          ! Cloud type [code]
          call fetch_values(this%full_grib_path, igrib, inx, iny, dum1d, &
               this%cldtyp(:, :, level), counter, rc)
       else if (param .eq. 228 .and. leveltype .eq. 109 .and. &
            level .gt. 0 .and. level .lt. 5) then
          ! Cloud top (m)
          call fetch_values(this%full_grib_path, igrib, inx, iny, dum1d, &
               this%cldtop(:, :, level), counter, rc)
       else if (param .eq. 183 .and. leveltype .eq. 200 .and. &
            level .eq. 0) then
          ! Time of last update from base time (min)
          call fetch_values(this%full_grib_path, igrib, inx, iny, dum1d, &
               this%pixtim, counter, rc)
       end if

       if (rc .ne. 0) then
          write(LIS_logunit,*)'[WARN] Problem reading values in ', &
               trim(this%full_grib_path)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if

       ! Close the GRIB message
       call grib_release(igrib, rc)
       if (rc .ne. 0) then
          write(LIS_logunit,*)'[WARN] Problem releasing message in ', &
               trim(this%full_grib_path)
          call grib_close_file(ftn)
          deallocate(dum1d)
          rc = 1
          return
       end if

       ! We can stop reading if we found all the variables we need.
       if (counter .gt. (3*NZ + 1)) exit

    end do ! k

    call grib_close_file(ftn)
#endif

    ! Make sure we found all the fields we were looking for.
    rc = 0
    if (counter .ne. (3*NZ + 1)) then
       write(LIS_logunit,*) 'Missing WWMCA data in ', &
            trim(this%full_grib_path)
       rc = 1
    end if

    ! Clean up
    deallocate(dum1d)

  contains

    ! Internal routine. Fetches value field from GRIB message.
    subroutine fetch_values(full_grib_path, igrib, inx, iny, data1d, data2d, &
         counter, rc)

      ! Imports
#if (defined USE_GRIBAPI)
      use grib_api
#endif
      use LIS_logMod, only: LIS_logunit

      ! Defaults
      implicit none

      ! Arguments
      character(*), intent(in) :: full_grib_path
      integer, intent(in) :: igrib
      integer, intent(in) :: inx
      integer, intent(in) :: iny
      real, intent(inout) :: data1d(inx*iny)
      real, intent(inout) :: data2d(inx, iny)
      integer, intent(inout) :: counter
      integer, intent(inout) :: rc

      ! Locals
      integer :: i
      integer :: j

#if (defined USE_GRIBAPI)

      call grib_get(igrib, 'values', data1d, rc)
      if (rc .ne. 0) then
         write(LIS_logunit,*)'[WARN] Failed to read values from ', &
              trim(full_grib_path)
         rc = 1
         return
      end if

      do j = 1, iny
         do i = 1, inx
            data2d(i, j) = data1d(i + (j-1)*inx)
         end do ! i
      end do ! j

      counter = counter + 1
      rc = 0

#endif
    end subroutine fetch_values

  end subroutine read_grib

  ! Copy the WWMCA GRIB1 data into the legacy binary arrays.
  subroutine return_binary_fields(this, amounts, types, tops, times)

    ! Defaults
    implicit none

    ! Arguments
    class(wwmca_grib1_t), intent(in) :: this
    byte, intent(inout) :: amounts(4, 1024, 1024)
    byte, intent(inout) :: types(4, 1024, 1024)
    integer*2, intent(inout) :: tops(4, 1024, 1024)
    integer*4, intent(inout) :: times(1024, 1024)

    ! Locals
    integer :: i, j, k

    do k = 1, 4
       do j = 1, 1024
          do i = 1, 1024
             amounts(k, i, j) = this%cldamt(i, j, k)
          end do ! i
       end do ! j
    end do ! k

    do k = 1, 4
       do j = 1, 1024
          do i = 1, 1024
             types(k, i, j) = this%cldtyp(i, j, k)
          end do ! i
       end do ! j
    end do ! k

    do k = 1, 4
       do j = 1, 1024
          do i = 1, 1024
             tops(k, i, j) = this%cldtop(i, j, k)
          end do ! i
       end do ! j
    end do ! k

    do j = 1, 1024
       do i = 1, 1024
          times(i, j) = this%pixtim(i, j)
       end do ! j
    end do ! i

  end subroutine return_binary_fields

  ! Constructs WWMCA GRIB1 filename, which contains all WWMCA fields.
  subroutine USAF_wwmca_grib1_filename(fname, rootdir, dir, &
       use_timestamp, hemi, yr, mo, da, hr)

    ! Defaults
    implicit none

    ! Arguments
    character(*), intent(inout) :: fname
    character(*), intent(in) :: rootdir
    character(*), intent(in) :: dir
    integer, intent(in) :: use_timestamp
    integer, intent(in) :: hemi
    integer, intent(in) :: yr, mo, da, hr

    ! Locals
    character(1), parameter :: FHEMI(2) = (/'N', 'S'/)
    character(10) :: ftime1, ftime2

    write(unit=ftime2, fmt='(i4, i2.2, i2.2, i2.2)') yr, mo, da, hr

    if (use_timestamp .eq. 1) then
       write(unit=ftime1, fmt='(a1, i4, i2.2, i2.2, a1)') '/', yr, mo, da, '/'
       fname = trim(rootdir) // ftime1 // trim(dir) // '/WWMCA_' // &
            ftime2 // '_' // fhemi(hemi) // '_16_M.GR1'
    else
       fname = trim(rootdir) // trim(dir) // '/WWMCA_' // &
            ftime2 // '_' // fhemi(hemi) // '_16_M.GR1'
    end if
  end subroutine USAF_wwmca_grib1_filename

end module USAF_WWMCA_grib1Mod
