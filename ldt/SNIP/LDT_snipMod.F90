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

  ! Defaults
  implicit none
  private

  ! Public methods
  public :: LDT_snipInit
  public :: LDT_snipRun

  ! Public type.  These are formerly environment and namelist variables from
  ! the original SNODEP.
  type, public :: snip_t
     ! Former environment variables
     character*10  :: date10
     character*255 :: fracdir
     character*255 :: modif
     integer :: sfcobsfmt ! EMK 20230727
     character*255 :: sfcobs
     character*255 :: ssmis
     character*255 :: gmi    !kyh20201118
     character*255 :: amsr2  !kyh20201217
     character*255 :: stmpdir
     character*255 :: sstdir ! EMK 20220113
     character*255 :: static
     character*255 :: unmod
     character*255 :: viirsdir

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
     logical :: usefrac
     logical :: useviirs

     ! option for brightness temperature data
     integer       :: TB_option          !kyh20201118

     ! option for PMW snow depth retrieval algorithms
     integer       :: ssmis_option
     character*255 :: ssmis_raw_dir
     character*255 :: gmi_raw_dir        !kyh20201118
     character*255 :: amsr2_raw_dir      !kyh20201217
     character*255 :: ff_file
     character*255 :: fd_file            !kyh20210113

     ! Bratseth settings
     real :: ob_err_var
     real :: back_err_var
     real :: back_err_h_corr_len
     real :: back_err_v_corr_len
     real :: elevqc_diff_threshold
     real :: skewed_backqc_threshold

     ! Other new settings
     real :: fill_climo
     character*255 :: source_of_ocean_data ! EMK 20240718
     character*255 :: gofs_sst_dir
     character*255 :: gofs_cice_dir
     character*255 :: espcd_sst_dir  ! EMK 20240718
     character*255 :: espcd_cice_dir ! EMK 20240718
     character*255 :: lis_grib2_dir
     character*20 :: security_class
     character*20 :: data_category
     character*20 :: data_res
     character*20 :: area_of_data
     character*255 :: galwem_root_dir
     character*255 :: galwem_sub_dir
     integer :: use_timestamp
     integer :: galwem_res

     ! Output file name (prefix)
     character*20 :: netcdf_prefix

     ! option for snow climatology
     integer       :: climo_option

  end type snip_t
  type(snip_t), public :: snip_settings

contains

  ! Reads SNIP-specific entries from ldt.config
  subroutine LDT_snipInit()

    ! Imports
    use ESMF
    use LDT_coreMod, only: LDT_config
    use LDT_logMod, only: LDT_verify, LDT_logunit, LDT_endrun

    ! Defaults
    implicit none

    ! Local variables
    character(len=255) :: cfg_entry
    integer :: rc
    integer :: c,r

    ! *** Get former environment variables ***

    ! Get date10
    cfg_entry = "SNIP valid date (YYYYMMDDHH):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%date10, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get fracdir
    cfg_entry = "SNIP fractional snow data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%fracdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get modif
    cfg_entry = "SNIP modified data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%modif, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! New sfcobs format...EMK 20230728
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

!------------------------------------------------------------------------------kyh20201118
    ! Select TB data
    cfg_entry = "SNIP brightness temperature data option:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%TB_option, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get TB data
    if (snip_settings%TB_option == 1) then ! SSMIS
       cfg_entry = "SNIP SSMIS data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%ssmis, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (snip_settings%TB_option == 2) then ! XCAL GMI
       cfg_entry = "SNIP XCAL GMI data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%gmi, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (snip_settings%TB_option == 3) then ! AMSR2
       cfg_entry = "SNIP AMSR2 data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%amsr2, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if
!------------------------------------------------------------------------------kyh20201118

    ! Get stmpdir
    cfg_entry = "SNIP surface temperature data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%stmpdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! EMK 20220113
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

    ! Get unmod
    cfg_entry = "SNIP unmodified data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%unmod, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get viirsdir
    cfg_entry = "SNIP VIIRS data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%viirsdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

!------------------------------------------------------------------------------kyh20201118
    ! *** for PMW snow depth retrieval, Yeosang Yoon
    ! get PMW raw datasets
    if (snip_settings%TB_option == 1) then ! SSMIS
       cfg_entry = "SNIP SSMIS raw data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%ssmis_raw_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (snip_settings%TB_option == 2) then ! XCAL GMI
       cfg_entry = "SNIP XCAL GMI raw data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%gmi_raw_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (snip_settings%TB_option == 3) then ! AMSR2
       cfg_entry = "SNIP AMSR2 raw data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%amsr2_raw_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if
!------------------------------------------------------------------------------kyh20201118

    ! YY: get option for PMW snow depth retrieval alogrithm
    cfg_entry = "SNIP PMW snow depth retrieval algorithm option:"  
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%ssmis_option, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! get forest fraction for algorithm 3 and 4
    if (snip_settings%ssmis_option==3 .or. snip_settings%ssmis_option==4) then   !kyh20201212
       cfg_entry = "SNIP forest fraction file:"   !YY
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%ff_file, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if

    ! get forest density for algorithm 4
    if (snip_settings%ssmis_option==4) then   !kyh20210113
       cfg_entry = "SNIP forest density file:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, snip_settings%fd_file, & 
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if

    ! get option for snow climatology
    cfg_entry = "SNIP Snow Climatology:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%climo_option,&
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
    cfg_entry = "SNIP SSMIS shallow snow depth threshold (m):"
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
    cfg_entry = "SNIP min VIIRS/CDFS-II fraction to mark point as snow:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%minfrac,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get usefrac
    cfg_entry = "SNIP use CDFS-II fractional snow data:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, snip_settings%usefrac,&
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

    ! EMK 20240718...Specify source of ocean data.
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

    ! Get output file name (prefix)
    cfg_entry = "SNIP netcdf filename prefix:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         snip_settings%netcdf_prefix, &
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
