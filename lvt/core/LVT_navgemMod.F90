!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"

!------------------------------------------------------------------------------
!
! MODULE: LVT_navgemMod
!
! DESCRIPTION:
! Contains routines for reading skin temperature and sea ice from
! NAVGEM GRIB files at 0.25 deg resolution, and interpolate to
! LVT grid.  Intended to run as part of 557post mode for Air Force operations.
!
! REVISION HISTORY:
! 10 May 2021: Eric Kemp (SSAI), Initial implementaton. Basic logic for
!              pulling merged SST/skin temperature and interpolating,
!              based on sample file provided by FNMOC. Still TODO: Handling
!              sea ice (thickness and areal fraction), and finalizing
!              file name convention.
!------------------------------------------------------------------------------

module LVT_navgemMod

  ! Defaults
  implicit none
  private

  ! Public routines
  public :: LVT_get_navgem_sst_gr1_filename
  public :: LVT_fetch_navgem_sst_gr1_field
  public :: LVT_get_navgem_sst_bin_filename
  public :: LVT_fetch_navgem_sst_bin_field

contains

  subroutine construct_navgem_sst_gr1_filename(rootdir, &
       year, month, day, hour, &
       fcst_hr, filename)

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in)  :: rootdir
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    integer, intent(in) :: fcst_hr
    character(len=*), intent(out) :: filename

    ! Local variables
    character(len=10) :: yyyymmddhh
    character(len=4)  :: hhhh

    write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') year, month, day, hour
    write(hhhh,'(i4.4)') fcst_hr

    ! FIXME:  Update file name to match that provided by 557WW.  The
    ! existing code is for a sample file provided by FNMOC.
    filename = trim(rootdir) // '/NAVGEM-' &
         // yyyymmddhh &
         // '-global_1440x721-grnd_sea_temp-surface-' &
         // '00000000-00000000-fcst_ops-' &
         // hhhh // '.gr1'

  end subroutine construct_navgem_sst_gr1_filename

  subroutine LVT_get_navgem_sst_gr1_filename(filename, &
       year, month, day, hour, fcst_hr)

    ! Modules
    use LVT_coreMod, only: LVT_rc
    use LVT_logMod, only: LVT_logunit
    use LVT_timeMgrMod, only: LVT_get_julhr, LVT_julhr_date

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(inout) :: filename
    integer, intent(out) :: year
    integer, intent(out) :: month
    integer, intent(out) :: day
    integer, intent(out) :: hour
    integer, intent(out) :: fcst_hr

    ! Locals
    integer :: navgem_julhr, lvt_julhr
    logical :: file_exists

    ! EMK TEST...Hardwired attributes for sample NAVGEM data from FNMOC.
    ! year = 2021
    ! month = 04
    ! day = 13
    ! hour = 00
    ! fcst_hr = 00

    ! First guess for NAVGEM run
    call LVT_get_julhr(LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
         LVT_rc%hr, 0, 0, lvt_julhr)
    select case (LVT_rc%hr)
    case(21)
       navgem_julhr = lvt_julhr - 3
       fcst_hr = 3
    case (18)
       navgem_julhr = lvt_julhr
       fcst_hr = 0
    case (15)
       navgem_julhr = lvt_julhr - 3
       fcst_hr = 3
    case (12)
       navgem_julhr = lvt_julhr
       fcst_hr = 0
    case (9)
       navgem_julhr = lvt_julhr - 3
       fcst_hr = 3
    case (6)
       navgem_julhr = lvt_julhr
       fcst_hr = 0
    case (3)
       navgem_julhr = lvt_julhr - 3
       fcst_hr = 3
    case (0)
       navgem_julhr = lvt_julhr
       fcst_hr = 0
    case default
       write(LVT_logunit,*)'[ERR] Invalid hour for LVT postprocessing!'
       write(LVT_logunit,*)'LVT_rc%hr = ', LVT_rc%hr
       stop
    end select

    call LVT_julhr_date(navgem_julhr, year, month, day, hour)

    call construct_navgem_sst_gr1_filename('./navgem', &
         year, month, day, hour, fcst_hr, filename)

    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
       write(LVT_logunit,*)'[INFO] Will use ', trim(filename)
       return
    end if

    ! At this point, we are rolling back to earlier NAVGEM file
    ! Start looping for earlier files
    do
       write(LVT_logunit,*) '[WARN] Cannot find ', trim(filename)
       navgem_julhr = navgem_julhr - 6
       fcst_hr = fcst_hr + 6

       ! NAVGEM only produces SST out to 180 hours (as of 12 May 2021).
       ! So that is a good criteria to give up.
       if (fcst_hr .gt. 180) then
          write(LVT_logunit,*) &
               '[WARN] *** GIVING UP ON NAVGEM SST ***'
          write(LVT_logunit,*) &
               '[WARN] *** NO NAVGEM SST DATA AVAILABLE!!!'
          filename = 'NONE'
          return
       end if
       call LVT_julhr_date(navgem_julhr, year, month, day, hour)

       call construct_navgem_sst_gr1_filename('./navgem', &
            year, month, day, hour, fcst_hr, filename)
       inquire(file=trim(filename), exist=file_exists)
       if (file_exists) then
          write(LVT_logunit,*)'[INFO] Will use ', trim(filename)
          return
       end if
    end do

  end subroutine LVT_get_navgem_sst_gr1_filename

  ! Routine for fetching merged sea surface temperature/land surface
  ! field.  We refer to this as "SST" for simplicity.
  subroutine LVT_fetch_navgem_sst_gr1_field(filename, sst, gridDesc)

    ! Modules
    use grib_api
    use LVT_logMod, only: LVT_logunit, LVT_verify

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: filename
    real, allocatable, intent(out) :: sst(:)
    real, intent(out) :: gridDesc(50)

    ! Locals
    integer :: year, month, day, hour, fcst_hr
    logical :: file_exists
    integer :: ftn, iret, igrib
    integer :: nvars, ivar
    integer :: editionNumber, indicatorOfParameter, centre, table2Version, &
         generatingProcessIdentifier, indicatorOfTypeOfLevel, level, &
         dataRepresentationType, Ni, Nj, &
         latitudeOfFirstGridPoint, longitudeOfFirstGridPoint, &
         latitudeOfLastGridPoint, longitudeOfLastGridPoint, &
         iDirectionIncrement, jDirectionIncrement, scanningMode
    integer :: i_180, i, j, ilon, iNewStartLon, iNewEndLon, i_rotate
    real, allocatable :: sst1(:,:), sst2(:,:)

    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) return
    write(LVT_logunit,*)'[INFO] Processing ', trim(filename)

    call grib_open_file(ftn, trim(filename), 'r', iret)
    call LVT_verify(iret, '[ERR] Bad return from grib_open_file')

    call grib_count_in_file(ftn, nvars, iret)
    call LVT_verify(iret, '[ERR] Bad return from grib_count_in_file')

    ! NOTE: Below code assumes GRIB1.
    do ivar = 1, nvars

       ! Get next GRIB message
       call grib_new_from_file(ftn, igrib, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_new_from_file')

       ! Make sure this is a GRIB1 message
       call grib_get(igrib, 'editionNumber', editionNumber, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')
       if (editionNumber .ne. 1) then
          write(LVT_logunit,*) &
               '[WARN] Expected GRIB1 message, found ', editionNumber
          cycle
       end if

       ! Check the parameter number
       call grib_get(igrib, 'indicatorOfParameter', indicatorOfParameter, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')
       if (indicatorOfParameter .ne. 133) cycle

       ! At this point we think we have SST.  We sanity
       ! check by looking at the center, table version, etc.
       call grib_get(igrib, 'centre', centre, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')
       if (centre .ne. 58) cycle ! Not from FNMOC

       call grib_get(igrib, 'table2Version', table2Version, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')
       if (table2Version .ne. 3) cycle ! Wrong table version

       call grib_get(igrib, 'generatingProcessIdentifier', &
            generatingProcessIdentifier, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')
       if (generatingProcessIdentifier .ne. 18) cycle ! Not from NAVGEM

       call grib_get(igrib, 'indicatorOfTypeOfLevel', &
            indicatorOfTypeOfLevel, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')
       if (indicatorOfTypeOfLevel .ne. 1) cycle ! Not at ground or water sfc

       call grib_get(igrib, 'level', level, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')
       if (level .ne. 0) cycle ! Not at ground or water surface

       ! Confidence is high this is the SST field.  Now we assemble the
       ! map projection data.
       call grib_get(igrib, 'dataRepresentationType', &
            dataRepresentationType, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'Ni', Ni, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'Nj', Nj, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'latitudeOfFirstGridPoint', &
            latitudeOfFirstGridPoint, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'longitudeOfFirstGridPoint', &
            longitudeOfFirstGridPoint, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'latitudeOfLastGridPoint', &
            latitudeOfLastGridPoint, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'longitudeOfLastGridPoint', &
            longitudeOfLastGridPoint, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'iDirectionIncrement', &
            iDirectionIncrement, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'jDirectionIncrement', &
            jDirectionIncrement, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       call grib_get(igrib, 'scanningMode', &
            scanningMode, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       ! We may need to rotate the grid to fit the longitude limits of
       ! -180 to 180 E.
       i_180 = 0
       iNewStartLon = longitudeOfFirstGridPoint
       iNewEndLon = longitudeOfLastGridPoint
       ilon = longitudeOfFirstGridPoint - iDirectionIncrement
       do i = 1, Ni
          ilon = ilon + iDirectionIncrement
          if (ilon .ge. 180000 .and. &
               (ilon - iDirectionIncrement) .lt. 180000) then
             i_180 = i
             iNewStartLon = ilon - 360000
             iNewEndLon = ilon - iDirectionIncrement
             exit
          end if
       end do

       ! Construct the grid description.  Look in lis/core/LIS_PRIV_rcMod.F90
       ! for info on the grid description array.
       if (dataRepresentationType .ne. 0) then
          write(LVT_logunit,*) &
               '[WARN] NAVGEM data not on lat/lon grid, will skip...'
          cycle
       end if
       griddesc = 0
       griddesc(1) = 0 ! Lat/lon projection
       griddesc(2) = Ni
       griddesc(3) = Nj
       griddesc(4) = latitudeOfFirstGridPoint * 1.e-3
       griddesc(5) = iNewStartLon * 1.e-3 ! Possibly rotated grid
       gridDesc(6) = 128
       gridDesc(7) = latitudeOfLastGridPoint * 1.e-3
       gridDesc(8) = iNewEndLon * 1.e-3 ! Possibly rotated grid
       gridDesc(9) = iDirectionIncrement * 1.e-3
       gridDesc(10) = jDirectionIncrement * 1.e-3
       gridDesc(11) = 64
       ! If scanning mode flag bit 3 is 0, data are in E-W ordering, otherwise
       ! in N-S ordering.  See NCEP ON388 Table 8.
       ! Note that GRIB1 bit ordering is right-to-left and starts at 1, while
       ! Fortran bit ordering is left-to-right and starts at 0.  So, GRIB1
       ! bit 3 is Fortran bit 5.
       if (.not. btest(scanningMode, 5)) then
          gridDesc(20) = 64  ! E-W ordering
       else
          gridDesc(20) = 255 ! N-S ordering
       end if
       gridDesc(30) = 0 ! Lat/lon projection
       griddesc(32) = Ni
       griddesc(33) = Nj
       griddesc(34) = latitudeOfFirstGridPoint * 1.e-3
       griddesc(35) = iNewStartLon * 1.e-3 ! Possibly rotated grid
       gridDesc(36) = 128
       gridDesc(37) = latitudeOfLastGridPoint * 1.e-3
       gridDesc(38) = iNewEndLon * 1.e-3 ! Possibly rotated grid
       gridDesc(39) = iDirectionIncrement * 1.e-3
       gridDesc(40) = jDirectionIncrement * 1.e-3

       ! Get the values
       allocate(sst(Ni*Nj))
       call grib_get(igrib, 'values', sst, iret)
       call LVT_verify(iret, '[ERR] Bad return from grib_get')

       ! See if we need to rotate the data
       if (iNewStartLon .ne. longitudeOfFirstGridPoint .or. &
            iNewEndLon .ne. longitudeOfLastGridPoint) then
          allocate(sst1(Ni,Nj))
          sst1 = 0
          allocate(sst2(Ni,Nj))
          sst2 = 0

          do j = 1, Nj
             do i = 1, Ni
                sst1(i,j) = sst(i + (j-1)*Ni)
             end do
          end do
          do j = 1, Nj
             do i = 1, Ni
                i_rotate = i + i_180 - 1
                if (i_rotate .gt. Ni) then
                   i_rotate = i_rotate - Ni
                end if
                sst2(i_rotate, j) = sst1(i,j)
             end do
          end do
          do j = 1, Nj
             do i = 1, Ni
                sst(i + (j-1)*Ni) = sst2(i,j)
             end do
          end do

          deallocate(sst1)
          deallocate(sst2)
       end if

       ! If we reached this point, we have what we need.
       exit
    end do

    call grib_close_file(ftn)

  end subroutine LVT_fetch_navgem_sst_gr1_field

  subroutine construct_navgem_sst_bin_filename(rootdir, &
       year, month, day, hour, filename)

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in)  :: rootdir
    integer, intent(in) :: year
    integer, intent(in) :: month
    integer, intent(in) :: day
    integer, intent(in) :: hour
    character(len=*), intent(out) :: filename

    ! Local variables
    character(len=10) :: yyyymmddhh

    write(yyyymmddhh,'(i4.4,i2.2,i2.2,i2.2)') year, month, day, hour

    ! FIXME:  Update file name to match that provided by 557WW.  The
    ! existing code is for a sample file provided by FNMOC.
    filename = trim(rootdir) // '/navyssts_0p25deg.' &
         // yyyymmddhh &
         // '.dat'

  end subroutine construct_navgem_sst_bin_filename

  ! Select binary file with Navy 0.25deg SST data.  Based on SNODEP/USAFSI
  ! code.
  subroutine LVT_get_navgem_sst_bin_filename(filename, &
       year, month, day, hour)

    ! Modules
    use LVT_logMod, only: LVT_logunit
    use LVT_timeMgrMod, only: LVT_get_julhr, LVT_julhr_date

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(inout) :: filename
    integer, intent(out) :: year
    integer, intent(out) :: month
    integer, intent(out) :: day
    integer, intent(out) :: hour

    ! Locals
    integer :: navgem_julhr, lvt_julhr
    logical :: file_exists

    ! FIXME...Add dynamic search for nearest NAVGEM file.  The
    ! existing code is hardwired for a sample file.
    year = 2021
    month = 04
    day = 13
    hour = 12

    call construct_navgem_sst_bin_filename('./navgem', &
         year, month, day, hour, filename)

    write(LVT_logunit,*)'[INFO] *** Searching for NAVGEM file ', &
         trim(filename)
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
       write(LVT_logunit,*)'[INFO] Will use ', trim(filename)
       return
    end if

    ! FIXME...Add dynamic search for NAVGEM file
    write(LVT_logunit,*)'[ERR] Cannot find NAVGEM file!'
    stop

  end subroutine LVT_get_navgem_sst_bin_filename

  ! Binary file reader. Taken from SNODEP/USAFSI
  subroutine putget_real(buffer, iofunc, file_name, igrid, jgrid)

    ! Modules
    use LVT_logMod, only: LVT_logunit

    ! Defaults
    implicit none

    ! Arguments
    real, intent(inout) :: buffer(igrid, jgrid)
    character*1, intent(in) :: iofunc
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: igrid
    integer, intent(in) :: jgrid

    ! Locals
    integer, parameter :: msglns = 20
    character*14 :: access_type
    character*4 :: cstat
    character*90 :: message(msglns)
    integer :: istat
    integer :: istat1
    integer :: reclen
    logical :: isopen
    integer :: i

    isopen  = .false.
    message = ' '
    reclen  = igrid * jgrid * 4

    access_type = '[INFO] Opening'
    open( 2, file=trim(file_name), form='unformatted', convert='little_endian', &
         access='direct', recl=reclen, iostat=istat, err=1000 )
    isopen = .true.

    ! Read from file, abort on error
    if ( (iofunc .eq. 'r') .or. (iofunc .eq. 'R') ) then
       access_type = '[INFO] Reading'
       read( 2, rec=1, iostat=istat, err=1000 ) buffer

       ! Write to file, abort on error
    else if ( (iofunc .eq. 'w') .or. (iofunc .eq. 'W') ) then
       access_type = '[INFO] Writing'
       write( 2, rec=1, iostat=istat, err=1000 ) buffer

       ! Else abort due to invalid IOFUNC value.
    else
       go to 4000
    endif

    close(2)

    return

    ! Error handling
1000 continue
    if (isopen) close (2)
    message(1) = '[ERR] LIBROUTINE: PUTGET_REAL'
    message(2) = '[ERR] ERROR ' // access_type // ' ' // trim (file_name)
    write (cstat, '(i4)', iostat=istat1) istat
    if (istat1 .eq. 0) message(3) = 'ISTAT = ' // cstat
    do i = 1, 3
       write(LVT_logunit,*) trim(message(i))
    end do
    return

4000 continue
    close(2)
    message(1) = '[ERR] LIBROUTINE: PUTGET_REAL'
    message(2) = '[ERR] INVALID IOFUNC VALUE = '// trim(iofunc)
    do i = 1, 2
       write(LVT_logunit,*) trim(message(i))
    end do
    return

  end subroutine putget_real

  ! Routine for fetching merged sea surface temperature/land surface
  ! field.  We refer to this as "SST" for simplicity.
  subroutine LVT_fetch_navgem_sst_bin_field(filename, sst, gridDesc)

    ! Modules
    use LVT_coreMod, only: LVT_rc
    use LVT_logMod, only: LVT_logunit, LVT_verify

    ! Defaults
    implicit none

    ! Arguments
    character(len=*), intent(in) :: filename
    real, allocatable, intent(out) :: sst(:)
    real, intent(out) :: gridDesc(50)

    ! Locals
    integer, parameter :: sst_igrid = 1440
    integer, parameter :: sst_jgrid = 721
    integer, parameter :: sst_size = sst_igrid * sst_jgrid
    integer :: year, month, day, hour, fcst_hr
    logical :: file_exists
    real, allocatable :: sst2d(:,:)
    integer :: i, j, i_180, i_rotate

    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) return
    write(LVT_logunit,*)'[INFO] Processing ', trim(filename)

    allocate(sst2d(sst_igrid, sst_jgrid))
    sst2d = LVT_rc%udef
    call putget_real(sst2d, 'r', trim(filename), sst_igrid, sst_jgrid)

    ! This grid ranges from 0 to 360E, but we need to rotate to -180 to 180
    ! E.
    allocate(sst(sst_size))
    sst = -1
    i_180 = 720
    do j = 1, sst_jgrid
       do i = 1, sst_igrid
          i_rotate = i + i_180 - 1
          if (i_rotate .gt. sst_igrid) then
             i_rotate = i_rotate - sst_igrid
          end if
          sst(i_rotate + (j-1)*sst_igrid) = sst2d(i,j)
       end do
    end do
    deallocate(sst2d)

    griddesc = 0
    griddesc(1) = 0 ! Lat/lon projection
    griddesc(2) = sst_igrid
    griddesc(3) = sst_jgrid
    griddesc(4) = -90.
    griddesc(5) = -180.! Rotated grid
    gridDesc(6) = 128
    gridDesc(7) = 90.
    gridDesc(8) = 179.75 ! Rotated grid
    gridDesc(9) = 0.25
    gridDesc(10) = 0.25
    gridDesc(11) = 64
    gridDesc(20) = 64  ! E-W ordering
    gridDesc(30) = 0 ! Lat/lon projection
    griddesc(32) = sst_igrid
    griddesc(33) = sst_jgrid
    griddesc(34) = -90.
    griddesc(35) = -180. ! Rotated grid
    gridDesc(36) = 128
    gridDesc(37) = 90.
    gridDesc(38) = 179.75 ! Rotated grid
    gridDesc(39) = 0.25
    gridDesc(40) = 0.25

  end subroutine LVT_fetch_navgem_sst_bin_field

end module LVT_navgemMod
