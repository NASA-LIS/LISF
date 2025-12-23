!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_snipMod

  ! Imports
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN

  ! Defaults
  implicit none
  private

  ! Public methods
  public :: LDT_snipInit
  public :: LDT_snipRun

  ! Public type.  These are formerly environment and namelist variables
  ! from the original SNODEP.
  type, public :: snip_t

     character*10  :: date10
     character(LDT_CONST_PATH_LEN) :: snip_dir
     character(LDT_CONST_PATH_LEN) :: usafsi_dir
     character(LDT_CONST_PATH_LEN) :: snodep_modif
     character(LDT_CONST_PATH_LEN) :: snodep_unmod
     integer :: sfcobsfmt
     character(LDT_CONST_PATH_LEN) :: sfcobs
     character(LDT_CONST_PATH_LEN) :: stmpdir
     character(LDT_CONST_PATH_LEN) :: sstdir
     character(LDT_CONST_PATH_LEN) :: static

     character(LDT_CONST_PATH_LEN) :: viirsdir
     character(LDT_CONST_PATH_LEN) :: amsr2dir

     ! Former namelist variables
     real :: clmadj
     real :: unkdep
     real :: minprt
     integer :: maxsobs
     real :: minsat
     integer :: trplat(3)
     integer :: elvlim(4)
     integer :: thresh
     integer :: arctmax
     integer :: minice
     real :: icelat(12,2)
     real :: latchk(12,2)
     integer :: maxpixage
     real :: minbare
     real :: minfrac
     logical :: useviirs

     ! option for brightness temperature data
     integer       :: TB_option          !kyh20201118

     ! option for PMW snow depth retrieval algorithms

     ! Bratseth settings
     real :: ob_err_var
     real :: back_err_var
     real :: back_err_h_corr_len
     real :: back_err_v_corr_len
     real :: elevqc_diff_threshold
     real :: skewed_backqc_threshold

     ! Other new settings
     real :: fill_climo
     character(LDT_CONST_PATH_LEN) :: source_of_ocean_data ! EMK 20240718
     character(LDT_CONST_PATH_LEN) :: gofs_sst_dir
     character(LDT_CONST_PATH_LEN) :: gofs_cice_dir
     character(LDT_CONST_PATH_LEN) :: espcd_sst_dir  ! EMK 20240718
     character(LDT_CONST_PATH_LEN) :: espcd_cice_dir ! EMK 20240718
     character(LDT_CONST_PATH_LEN) :: lis_grib2_dir
     character*20 :: security_class
     character*20 :: data_category
     character*20 :: data_res
     character*20 :: area_of_data
     character(LDT_CONST_PATH_LEN) :: galwem_root_dir
     character(LDT_CONST_PATH_LEN) :: galwem_sub_dir
     integer :: use_timestamp
     integer :: galwem_res

     ! Output file name (prefix)
     character*20 :: netcdf_prefix_snip
     character*20 :: netcdf_prefix_usafsi

     ! option for snow climatology
     integer       :: climo_option

  end type snip_t
  type(snip_t), public :: snip_settings

contains

  ! Reads SNIP-specific entries from ldt.config
  subroutine LDT_snipInit()

    ! Imports
    use ESMF
    use LDT_constantsMod, only: LDT_CONST_PATH_LEN
    use LDT_coreMod, only: LDT_config
    use LDT_logMod, only: LDT_verify, LDT_logunit, LDT_endrun

    ! Defaults
    implicit none

    ! Local variables
    character(len=LDT_CONST_PATH_LEN) :: cfg_entry
    integer :: rc
    integer :: c,r



    ! Get date10
    cfg_entry = "SNIP valid date (YYYYMMDDHH):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%date10, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get SNIP
    cfg_entry = "SNIP netcdf data directory for SNIP:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%snip_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get USAFSI
    cfg_entry = "SNIP netcdf data directory for USAFSI:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%usafsi_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get modif
    cfg_entry = "SNIP data directory for modified SNODEP:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%snodep_modif, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get unmod
    cfg_entry = "SNIP data directory for unmodified SNODEP:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%snodep_unmod, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! New sfcobs format for WIGOS
    cfg_entry = "SNIP surface obs data format:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%sfcobsfmt, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get sfcobs
    cfg_entry = "SNIP surface obs data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%sfcobs, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get stmpdir
    cfg_entry = &
         "SNIP degribbed LIS 0.25 deg sfc temperature data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%stmpdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get FNMOC SST data
    cfg_entry = "SNIP FNMOC SST GRIB1 data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%sstdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get static
    cfg_entry = "SNIP static data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%static, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get viirsdir
    cfg_entry = "SNIP VIIRS data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%viirsdir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! TODO:  Directory for reading AMSR2 AI/ML retrievals
    ! Get viirsdir
    cfg_entry = "SNIP AMSR2 snowdepth data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%amsr2dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! get option for snow climatology
    cfg_entry = "SNIP Snow Climatology:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%climo_option, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! *** Get former namelist variables ***

    ! Get clmadj
    cfg_entry = &
         "SNIP decimal fraction adjustment of snow depth towards climo:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%clmadj,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get unkdep
    cfg_entry = "SNIP default snow depth (m) when actual depth unknown:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%unkdep,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minprt
    cfg_entry = &
         "SNIP minimum snow depth (m) for which to print a diagnostic:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%minprt,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get maxsobs
    cfg_entry = "SNIP maximum number of surface observations allowed:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%maxsobs,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minsat
    cfg_entry = "SNIP AMSR2 shallow snow depth threshold (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%minsat,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get trplat
    cfg_entry = "SNIP latitudes (deg * 100) for summer climo check:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do c = 1,3
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%trplat(c),&
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end do ! c

    ! Get elvlim
    cfg_entry = "SNIP elevations (m) for summer climo check:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do c = 1,4
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%elvlim(c),&
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end do ! c

    ! Get thresh
    cfg_entry = &
         "SNIP temperature (deg K * 10) above which no snow is allowed:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%thresh,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get arctmax
    cfg_entry = &
        "SNIP max reported temperature (deg K * 10) allowed around poles:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%arctmax,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minice
    cfg_entry = &
         "SNIP minimum ice concentration (%) needed to set ice flag:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%minice,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get icelat
    cfg_entry = "SNIP high latitude thresholds (deg) for sea ice::"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do r = 1,2
       call ESMF_ConfigNextLine(LDT_config, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       do c = 1, 12
          call ESMF_ConfigGetAttribute(LDT_config, &
               snip_settings%icelat(c,r), rc=rc)
          call LDT_verify(rc, trim(cfg_entry)//" not specified")
       end do ! c
    end do ! r

    ! Get latchk
    cfg_entry = "SNIP low latitude thresholds (deg) for sea ice::"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do r = 1,2
       call ESMF_ConfigNextLine(LDT_config, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       do c = 1, 12
          call ESMF_ConfigGetAttribute(LDT_config, &
               snip_settings%latchk(c,r), rc=rc)
          call LDT_verify(rc, trim(cfg_entry)//" not specified")
       end do ! c
    end do ! r

    ! Get maxpixage
    cfg_entry = "SNIP max age of VIIRS pixels to use:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%maxpixage,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minbare
    cfg_entry = "SNIP min VIIRS fraction to mark point as bare ground:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%minbare,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minfrac
    cfg_entry = "SNIP min VIIRS fraction to mark point as snow:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%minfrac,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get useviirs
    cfg_entry = "SNIP use VIIRS snow mask:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%useviirs,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! ** New Bratseth settings **
    ! Get obs_err_var
    cfg_entry = "SNIP observation error variance (m^2):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%ob_err_var,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get back_err_var
    cfg_entry = "SNIP background error variance (m^2):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%back_err_var,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get back_err_h_corr_len
    cfg_entry = "SNIP background error horizontal correlation length (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%back_err_h_corr_len,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get back_err_v_corr_len
    cfg_entry = "SNIP background error vertical correlation length (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%back_err_v_corr_len,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get elev_diff_thresh
    cfg_entry = "SNIP elevQC difference threshold (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%elevqc_diff_threshold,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get skewed_backqc_threshold
    cfg_entry = "SNIP skewed backQC snow depth threshold (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%skewed_backqc_threshold,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Bogus value for climo if not defined at land point
    cfg_entry = "SNIP bogus climatology snow depth value (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%fill_climo,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Specify source of ocean data.
    cfg_entry = "SNIP source of ocean data:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%source_of_ocean_data, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    if (snip_settings%source_of_ocean_data .ne. "GOFS" .and. &
         snip_settings%source_of_ocean_data .ne. "ESPC-D") then
       write(LDT_logunit,*)'[ERR] Unrecognized source of ocean data'
       write(LDT_logunit,*)'[ERR] Must be GOFS or ESPC-D'
       write(LDT_logunit,*) &
            "[ERR] Update entry for 'SNIP source of ocean data:'"
       write(LDT_logunit,*)'[ERR] LDT will halt.'
       call LDT_endrun()
    end if

    if (snip_settings%source_of_ocean_data == "GOFS") then
       ! Get gofs_sst_dir
       cfg_entry = "SNIP GOFS SST data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, &
            snip_settings%gofs_sst_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")

       ! Get gofs_cice_dir
       cfg_entry = "SNIP GOFS CICE data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, &
            snip_settings%gofs_cice_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")

    else if (snip_settings%source_of_ocean_data == "ESPC-D") then

       ! Get espcd_sst_dir
       cfg_entry = "SNIP ESPC-D SST data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, &
            snip_settings%espcd_sst_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")

       ! Get espcd_cice_dir
       cfg_entry = "SNIP ESPC-D CICE data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, &
            snip_settings%espcd_cice_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")

    end if

    ! Get lis_grib2_dir
    cfg_entry = "SNIP LIS GRIB2 data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%lis_grib2_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get security class
    cfg_entry = "SNIP LIS GRIB2 security class:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%security_class, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get data category
    cfg_entry = "SNIP LIS GRIB2 data category:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%data_category, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get data resolution
    cfg_entry = "SNIP LIS GRIB2 data resolution:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%data_res, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get area of data
    cfg_entry = "SNIP LIS GRIB2 area of data:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%area_of_data, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get GALWEM root directory
    cfg_entry = "SNIP GALWEM root directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%galwem_root_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get GALWEM subdirectory
    cfg_entry = "SNIP GALWEM subdirectory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%galwem_sub_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! See if directory timestamp used for GALWEM
    cfg_entry = "SNIP GALWEM use timestamp directories:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%use_timestamp, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get GALWEM nominal resolution
    cfg_entry = "SNIP GALWEM nominal resolution (km):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%galwem_res, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get SNIP output file name (prefix)
    cfg_entry = "SNIP netcdf filename prefix for SNIP:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%netcdf_prefix_snip, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get USAFSI output file name (prefix)
    cfg_entry = "SNIP netcdf filename prefix for USAFSI:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%netcdf_prefix_usafsi, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

  end subroutine LDT_snipInit

  ! This calls the actual SNIP driver
  subroutine LDT_snipRun(n)
    implicit none
    integer, intent(in) :: n
    external :: SNIP_run
    call SNIP_run(n)
  end subroutine LDT_snipRun
end module LDT_snipMod
