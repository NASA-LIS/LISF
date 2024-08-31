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
! MODULE: USAFSI_galwemMod
!
! REVISION HISTORY:
! 28 Oct 2020  Eric Kemp  Initial version.
!
! DESCRIPTION:
! Source code for fetching 2-m air temperature from GALWEM.
!-------------------------------------------------------------------------

#include "LDT_misc.h"

module USAFSI_galwemMod

  ! Defaults
  implicit none
  private

  ! Public methods
  public :: USAFSI_get_galwem_t2m

contains

  subroutine USAFSI_get_galwem_t2m(n, julhr, nc, nr, t2m, rc)

    ! Imports
#if (defined USE_GRIBAPI)
    use grib_api
#endif
    use LDT_logMod, only: LDT_logunit
    use LDT_timeMgrMod, only: LDT_julhr_date
    use LDT_usafsiMod, only: usafsi_settings

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: n
    integer, intent(in) :: julhr
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(inout) :: t2m(nc,nr)
    integer, intent(inout) :: rc

    ! Local variables
    character(255) :: gribfilename
    real :: gridDesci_glb(50)
    logical :: file_exists
    logical :: found
    logical :: roll_back
    integer :: fc_hr
    integer :: file_julhr
    integer :: yr1, mo1, da1, hr1
    integer :: ftn
    integer :: ierr
    integer :: nvars
    integer :: igrib
    integer :: center
    character(100) :: gtype
    integer :: iginfo(40)
    real :: gridres_dlat
    real :: gridres_dlon
    integer :: dataDate
    integer :: dataTime
    integer :: ifguess
    integer :: jfguess
    logical :: found_inq
    integer :: k
    integer :: prod_def_tmpl_num
    integer :: param_disc_val
    integer :: param_cat_val
    integer :: param_num_val
    integer :: surface_val, surface_val_2
    integer :: level_val
    real, allocatable :: dum1d(:)
    real, allocatable :: fg_t2m(:,:)

    rc = 1 ! Update below if we succeed in interpolating temperature
    call set_galwem_griddesc(gridDesci_glb)

    ! Will search previous GALWEM cycles every six hours, up to 10 days,
    ! until we find an acceptable file.
    fc_hr = 0
    file_julhr = julhr
    call LDT_julhr_date(file_julhr, yr1, mo1, da1, hr1)
    found = .false.
    roll_back = .false.
    do while (.not. found)

       ! If necessary, roll back to a previous GALWEM cycle
       if (roll_back) then
          fc_hr = fc_hr + 6
          ! GALWEM only runs out to 240 hours.
          if (fc_hr > 240) exit
          file_julhr = file_julhr - 6
          call LDT_julhr_date(file_julhr, yr1, mo1, da1, hr1)
       end if
       roll_back = .true.

       call get_galwem_filename(gribfilename, &
            usafsi_settings%galwem_root_dir, &
            usafsi_settings%galwem_sub_dir, usafsi_settings%use_timestamp, &
            usafsi_settings%galwem_res, yr1, mo1, da1, hr1, fc_hr)

       ! Before using ECCODE/GRIB_API, see if the GRIB file exists.
       inquire(file=trim(gribfilename), exist=found_inq)
       if (.not. found_inq) then
          write(LDT_logunit,*) '[WARN] Cannot find file ' // trim(gribfilename)
          cycle
       end if

#if (defined USE_GRIBAPI)
       ! Open the GRIB file
       ierr = 0
       call grib_open_file(ftn, trim(gribfilename), 'r', ierr)
       if ( ierr .ne. 0 ) then
          call grib_close_file(ftn)
          write(LDT_logunit,*) '[WARN] Failed to open ' // trim(gribfilename)
          write(LDT_logunit,*) '[WARN] ierr = ', ierr
          flush(LDT_logunit)
          cycle
       end if

       ! Count number of messages
       write(ldt_logunit,*)'[INFO] Reading ', trim(gribfilename)
       call grib_count_in_file(ftn, nvars, ierr)
       if ( ierr .ne. 0 ) then
          write(LDT_logunit,*) '[WARN] in grib_count_in_file ' // &
               'in USAFSI_get_galwem_t2m'
          call grib_close_file(ftn)
          cycle
       endif

       ! Start searching through the messages
       do k = 1, nvars

          call grib_new_from_file(ftn, igrib, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) &
                  '[WARN] failed to read message from' // trim(gribfilename)
             call grib_close_file(ftn)
             exit
          endif

          ! Check the product definition template number
          call grib_get(igrib, 'productDefinitionTemplateNumber', &
               prod_def_tmpl_num, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'productDefinitionTemplateNumber in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          end if
          ! We want analysis or forecast at a horizontal level or in a
          ! horizontal layer at a point in time.
          if (prod_def_tmpl_num .ne. 0) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Check the discipline
          call grib_get(igrib, 'discipline', param_disc_val, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'discipline in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          end if
          ! We want a meteorological product.
          if (param_disc_val .ne. 0) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Check the parameter category
          call grib_get(igrib, 'parameterCategory', param_cat_val, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'parameterCategory in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          end if
          ! We want the temperature category.
          if (param_cat_val .ne. 0) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Check the parameter number
          call grib_get(igrib, 'parameterNumber', param_num_val, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'parameterNumber in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          end if
          ! We want temperature
          if (param_num_val .ne. 0) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Check the first surface type
          call grib_get(igrib, 'typeOfFirstFixedSurface', surface_val, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'typeOfFirstFixedSurface in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          end if
          ! We want specified height above ground (m)
          if (surface_val .ne. 103) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Check the second surface type
          call grib_get(igrib, 'typeOfSecondFixedSurface', surface_val_2, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'typeOfSecondFixedSurface in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          end if
          ! We want the second surface type to be missing, meaning the data
          ! is for a level and not a layer.
          if (surface_val_2 .ne. 255) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Check the surface level
          call grib_get(igrib, 'level', level_val, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'level in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          end if
          ! Surface level should be 2 meters
          if (level_val .ne. 2) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! At this point, we found instantaneous 2-m air temperature.
          ! Make a few more sanity checks.

          ! Check originating center
          call grib_get(igrib, 'centre', center, ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) '[WARN] in grib_get: ' // &
                  'centre in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif
          if (center .ne. 57) then
             call grib_release(igrib, ierr)
             cycle
          endif

          ! Check the grid type
          call grib_get(igrib, 'gridType', gtype, ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) '[WARN] in grid_get: ' // &
                  'gridtype in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif
          if (trim(gtype) .ne. "regular_ll") then
             call grib_release(igrib, ierr)
             cycle
          endif

          ! Check the valid date
          call grib_get(igrib, 'dataDate', dataDate, ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) '[WARN] in grid_get: ' // &
                  'dataDate in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif
          if ( yr1*10000+mo1*100+da1 .ne. dataDate) then
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Check the valid time
          call grib_get(igrib, 'dataTime', dataTime, ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) '[WARN] in grid_get: ' // &
                  'dataTime in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif
          if (  hr1*100 .ne. dataTime ) then
             call grib_release(igrib, ierr)
             cycle
          end if

          call grib_get(igrib, 'Ni', iginfo(1), ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) &
                  '[WARN] in grid_get: Ni in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif
          if (iginfo(1) .ne. gridDesci_glb(2)) then
             call grib_release(igrib, ierr)
             cycle
          end if

          call grib_get(igrib, 'Nj', iginfo(2), ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) &
                  '[WARN] in grid_get: Nj in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif
          if (iginfo(2) .ne. gridDesci_glb(3)) then
             call grib_release(igrib, ierr)
             cycle
          end if

          call grib_get(igrib, 'jDirectionIncrementInDegrees', &
               gridres_dlat, ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) '[WARN] in grid_get: ' // &
                  'jDirectionIncrementInDegrees ' // &
                  'in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif

          call grib_get(igrib, 'iDirectionIncrementInDegrees', &
               gridres_dlon, ierr)
          if ( ierr .ne. 0 ) then
             write(LDT_logunit,*) '[WARN] in grid_get: ' // &
                  'iDirectionIncrementInDegrees ' // &
                  'in USAFSI_get_galwem_t2m'
             call grib_release(igrib, ierr)
             cycle
          endif

          write(LDT_logunit,*) &
               '[INFO] FOUND 2-M AIR TEMPERATURE FROM UK UM (GALWEM) MODEL'
          write(LDT_logunit,*)'[INFO] GALWEM DELTA LAT IS ', &
               gridres_dlat,' DEGREES'
          write(LDT_logunit,*)'[INFO] GALWEM DELTA LON IS ', &
               gridres_dlon,' DEGREES'

          ! Read the field
          ifguess = iginfo(1)
          jfguess = iginfo(2)
          allocate(dum1d(ifguess*jfguess))
          dum1d(:) = 0
          call grib_get(igrib, 'values', dum1d, ierr)
          if (ierr .ne. 0) then
             write(LDT_logunit,*) &
                  '[WARN] Cannot read temperature from GRIB file ', &
                  trim(gribfilename)
             deallocate(dum1d)
             call grib_release(igrib, ierr)
             cycle
          end if

          ! Reshape to 2D array
          allocate(fg_t2m(ifguess, jfguess))
          fg_t2m(:,:) = 0
          fg_t2m = reshape(dum1d, (/ifguess, jfguess/))
          deallocate(dum1d)

          ! We're finished reading the file.
          call grib_release(igrib, ierr)
          call grib_close_file(ftn)
          found = .true.
          exit

       end do

#endif

    end do
    if (.not. found) then
       if (allocated(dum1d)) deallocate(dum1d)
       if (allocated(fg_t2m)) deallocate(fg_t2m)
       return
    end if

    ! Now we interpolate the field to the LDT grid
    call interp_galwem_t2m(n, gridDesci_glb, ifguess, jfguess, fg_t2m, &
         nc, nr, t2m)

    ! Clean up and exit
    deallocate(fg_t2m)
    rc = 0

  end subroutine USAFSI_get_galwem_t2m

  ! Private subroutine
  subroutine get_galwem_filename(filename, rootdir, dir, use_timestamp, &
       nominal_res_km, yr, mo, da, hr, fc_hr)

    ! Imports
    use LDT_logMod, only: LDT_logunit, LDT_endrun

    ! Defaults
    implicit none

    ! Arguments
    character(*), intent(inout) :: filename
    character(*), intent(in) :: rootdir
    character(*), intent(in) :: dir
    integer, intent(in) :: use_timestamp
    integer, intent(in) :: nominal_res_km
    integer, intent(in) :: yr, mo, da, hr
    integer, intent(in) :: fc_hr

    ! Local variables
    character(8) :: ftime1
    character(2) :: fhr
    character(3) :: fchr

    character(len=54) :: fname1
    character(len=20) :: fname2

    write (UNIT=fhr, FMT='(i2.2)') hr
    write (UNIT=fchr, FMT='(i3.3)') fc_hr

    ! Support GALWEM 17 or 10-km data
    if (nominal_res_km == 17) then
       fname1 = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GD_GR.C17KM_AR.GLOBAL_DD.'
    else if (nominal_res_km == 10) then
       fname1 = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GD_GR.C10KM_AR.GLOBAL_DD.'
    else
       write(LDT_logunit)'[ERR] Invalid nominal resolution for GALWEM!'
       write(LDT_logunit)'[ERR] Found ', nominal_res_km
       write(LDT_logunit)'[ERR] Only supports 17 and 10'
       call LDT_endrun()
    end if

    ! Make sure ftime1 is always initialized
    write (UNIT=ftime1, FMT='(i4, i2.2, i2.2)') yr, mo, da

    if (use_timestamp .eq. 1) then
       filename = trim(rootdir) // '/' // ftime1 // '/' // trim(dir) // &
            '/' // fname1 // ftime1 // '_CY.' // fhr // '_FH.' // &
            fchr // '_DF.GR2'
    else
       filename = trim(rootdir) // '/' // trim(dir) // '/' // &
            fname1 // ftime1 // '_CY.' // fhr // '_FH.' // fchr // '_DF.GR2'
    endif

  end subroutine get_galwem_filename


  ! Private method
  subroutine interp_galwem_t2m(n, gridDesci_glb, ifguess, jfguess, fg_field, &
       nc, nr, t2m)

    ! Imports
    use LDT_coreMod, only: LDT_rc, LDT_domain
    use LDT_logMod, only: LDT_logunit
    use LDT_usafsiMod, only: usafsi_settings

    ! Defaults
    implicit none

    ! Arguments
    integer, intent(in) :: n
    real, intent(in) :: gridDesci_glb(50)
    integer, intent(in) :: ifguess
    integer, intent(in) :: jfguess
    real, intent(in) :: fg_field(ifguess, jfguess)
    integer, intent(in) :: nc
    integer, intent(in) :: nr
    real, intent(inout) :: t2m(nc, nr)

    ! Local variables
    integer   :: mi, mo
    integer   :: k
    integer   :: i,j
    integer   :: iret
    integer   :: midway
    character(len=50) :: method
    real, allocatable :: var(:,:)
    logical*1, allocatable :: lb(:)
    logical*1, allocatable :: lo(:)
    integer, allocatable :: n11(:), n12(:), n21(:), n22(:)
    real, allocatable :: w11(:), w12(:), w21(:), w22(:)

    allocate(var(ifguess,jfguess))
    allocate(lb(ifguess*jfguess))
    allocate(lo(nc*nr))

    mi = ifguess * jfguess
    mo = nc * nr
    lb = .true.

    ! translate from (0,360) to (-180,180).
    midway = ifguess/2
    do j = 1, jfguess
       do i = 1, ifguess
          if ( (i+midway) < ifguess ) then
             var(i,j) = fg_field((i+midway),j)
          else
             var(i,j) = fg_field((i-midway+1),j)
          endif
       enddo
    enddo

    ! Set bilinear weights
    allocate(n11(nc*nr))
    allocate(n12(nc*nr))
    allocate(n21(nc*nr))
    allocate(n22(nc*nr))
    allocate(w11(nc*nr))
    allocate(w12(nc*nr))
    allocate(w21(nc*nr))
    allocate(w22(nc*nr))

    call bilinear_interp_input(n, gridDesci_glb,        &
         n11, n12, n21, n22, &
         w11, w12, w21, w22)

    ! Use bilinear interpolation
    call bilinear_interp(LDT_rc%gridDesc(n,:), lb, &
         var, lo, t2m, mi, mo, &
         LDT_domain(n)%lat, LDT_domain(n)%lon, &
         w11, w12, w21, w22, &
         n11, n12, n21, n22, -9999., iret)

    ! Clean up
    deallocate(var)
    deallocate(lb)
    deallocate(lo)
    deallocate(w11)
    deallocate(w12)
    deallocate(w21)
    deallocate(w22)
    deallocate(n11)
    deallocate(n12)
    deallocate(n21)
    deallocate(n22)

  end subroutine interp_galwem_t2m

  ! Private method
  subroutine set_galwem_griddesc(gridDesci_glb)

    ! Imports
    use LDT_logMod, only: LDT_logunit, LDT_endrun
    use LDT_usafsiMod, only: usafsi_settings

    ! Defaults
    implicit none

    ! Arguments
    real, intent(inout) :: gridDesci_glb(50)

    ! Set the weights for the interpolation.  This varies by GALWEM resolution
    if (usafsi_settings%galwem_res == 17) then
       gridDesci_glb = 0
       gridDesci_glb(1) = 0
       gridDesci_glb(2) = 1536
       gridDesci_glb(3) = 1152
       gridDesci_glb(4) = -89.9219
       gridDesci_glb(5) = -179.882813
       gridDesci_glb(6) = 128
       gridDesci_glb(7) = 89.9219
       gridDesci_glb(8) = 179.887
       gridDesci_glb(9) = 0.234378
       gridDesci_glb(10) = 0.15625
       gridDesci_glb(20) = 0
    else if (usafsi_settings%galwem_res == 10) then
       gridDesci_glb = 0
       gridDesci_glb(1) = 0
       gridDesci_glb(2) = 2560
       gridDesci_glb(3) = 1920
       gridDesci_glb(4) = -89.9531250
       gridDesci_glb(5) = -179.9296875
       gridDesci_glb(6) = 128
       gridDesci_glb(7) = 89.9531250
       gridDesci_glb(8) = 179.9296875
       gridDesci_glb(9) = 0.140625
       gridDesci_glb(10) = 0.093750
       gridDesci_glb(20) = 0
    else
       write(LDT_logunit)'[ERR] Invalid nominal resolution for GALWEM!'
       write(LDT_logunit)'[ERR] Found ', usafsi_settings%galwem_res
       write(LDT_logunit)'[ERR] Only supports 17 and 10'
       call LDT_endrun()
    end if
  end subroutine set_galwem_griddesc
end module USAFSI_galwemMod

