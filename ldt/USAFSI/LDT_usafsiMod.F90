!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_usafsiMod

  ! Defaults
  implicit none
  private

  ! Public methods
  public :: LDT_usafsiInit
  public :: LDT_usafsiRun

  ! Public type.  These are formerly environment and namelist variables from
  ! the original SNODEP.
  type, public :: usafsi_t
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
     character*255 :: gofs_sst_dir
     character*255 :: gofs_cice_dir
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

  end type usafsi_t
  type(usafsi_t), public :: usafsi_settings

contains

  ! Reads USAFSI-specific entries from ldt.config
  subroutine LDT_usafsiInit()

    ! Imports
    use ESMF
    use LDT_coreMod, only: LDT_config
    use LDT_logMod, only: LDT_verify

    ! Defaults
    implicit none

    ! Local variables
    character(len=255) :: cfg_entry
    integer :: rc
    integer :: c,r

    ! *** Get former environment variables ***

    ! Get date10
    cfg_entry = "USAFSI valid date (YYYYMMDDHH):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%date10, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get fracdir
    cfg_entry = "USAFSI fractional snow data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%fracdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get modif
    cfg_entry = "USAFSI modified data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%modif, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! New sfcobs format...EMK 20230728
    cfg_entry = "USAFSI surface obs data format:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%sfcobsfmt, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get sfcobs
    cfg_entry = "USAFSI surface obs data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%sfcobs, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

!------------------------------------------------------------------------------kyh20201118
    ! Select TB data
    cfg_entry = "USAFSI brightness temperature data option:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%TB_option, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get TB data
    if (usafsi_settings%TB_option == 1) then ! SSMIS
       cfg_entry = "USAFSI SSMIS data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%ssmis, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (usafsi_settings%TB_option == 2) then ! XCAL GMI
       cfg_entry = "USAFSI XCAL GMI data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%gmi, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (usafsi_settings%TB_option == 3) then ! AMSR2
       cfg_entry = "USAFSI AMSR2 data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%amsr2, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if
!------------------------------------------------------------------------------kyh20201118

    ! Get stmpdir
    cfg_entry = "USAFSI surface temperature data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%stmpdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! EMK 20220113
    cfg_entry = "USAFSI FNMOC SST GRIB1 data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%sstdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get static
    cfg_entry = "USAFSI static data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%static, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get unmod
    cfg_entry = "USAFSI unmodified data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%unmod, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get viirsdir
    cfg_entry = "USAFSI VIIRS data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%viirsdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

!------------------------------------------------------------------------------kyh20201118
    ! *** for PMW snow depth retrieval, Yeosang Yoon
    ! get PMW raw datasets
    if (usafsi_settings%TB_option == 1) then ! SSMIS
       cfg_entry = "USAFSI SSMIS raw data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%ssmis_raw_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (usafsi_settings%TB_option == 2) then ! XCAL GMI
       cfg_entry = "USAFSI XCAL GMI raw data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%gmi_raw_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    elseif (usafsi_settings%TB_option == 3) then ! AMSR2
       cfg_entry = "USAFSI AMSR2 raw data directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%amsr2_raw_dir, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if
!------------------------------------------------------------------------------kyh20201118

    ! YY: get option for PMW snow depth retrieval alogrithm
    cfg_entry = "USAFSI PMW snow depth retrieval algorithm option:"  
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%ssmis_option, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! get forest fraction for algorithm 3 and 4
    if (usafsi_settings%ssmis_option==3 .or. usafsi_settings%ssmis_option==4) then   !kyh20201212
       cfg_entry = "USAFSI forest fraction file:"   !YY
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%ff_file, &
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if

    ! get forest density for algorithm 4
    if (usafsi_settings%ssmis_option==4) then   !kyh20210113
       cfg_entry = "USAFSI forest density file:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%fd_file, & 
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end if

    ! get option for snow climatology
    cfg_entry = "USAFSI Snow Climatology:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%climo_option,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! *** Get former namelist variables ***

    ! Get clmadj
    cfg_entry = &
         "USAFSI decimal fraction adjustment of snow depth towards climo:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%clmadj,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get unkdep
    cfg_entry = "USAFSI default snow depth (m) when actual depth unknown:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%unkdep,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minprt
    cfg_entry = &
         "USAFSI minimum snow depth (m) for which to print a diagnostic:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%minprt,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get maxsobs
    cfg_entry = "USAFSI maximum number of surface observations allowed:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%maxsobs,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minsat
    cfg_entry = "USAFSI SSMIS shallow snow depth threshold (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%minsat,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get trplat
    cfg_entry = "USAFSI latitudes (deg * 100) for summer climo check:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do c = 1,3
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%trplat(c),&
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end do ! c

    ! Get elvlim
    cfg_entry = "USAFSI elevations (m) for summer climo check:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do c = 1,4
       call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%elvlim(c),&
            rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    end do ! c

    ! Get thresh
    cfg_entry = &
         "USAFSI temperature (deg K * 10) above which no snow is allowed:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%thresh,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")


    ! Get arctmax
    cfg_entry = &
         "USAFSI max reported temperature (deg K * 10) allowed around poles:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%arctmax,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minice
    cfg_entry = &
         "USAFSI minimum ice concentration (%) needed to set ice flag:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%minice,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get icelat
    cfg_entry = "USAFSI high latitude thresholds (deg) for sea ice::"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do r = 1,2
       call ESMF_ConfigNextLine(LDT_config, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       do c = 1, 12
          call ESMF_ConfigGetAttribute(LDT_config, &
               usafsi_settings%icelat(c,r), rc=rc)
          call LDT_verify(rc, trim(cfg_entry)//" not specified")
       end do ! c
    end do ! r

    ! Get latchk
    cfg_entry = "USAFSI low latitude thresholds (deg) for sea ice::"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    do r = 1,2
       call ESMF_ConfigNextLine(LDT_config, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       do c = 1, 12
          call ESMF_ConfigGetAttribute(LDT_config, &
               usafsi_settings%latchk(c,r), rc=rc)
          call LDT_verify(rc, trim(cfg_entry)//" not specified")
       end do ! c
    end do ! r

    ! Get maxpixage
    cfg_entry = "USAFSI max age of VIIRS pixels to use:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%maxpixage,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minbare
    cfg_entry = "USAFSI min VIIRS fraction to mark point as bare ground:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%minbare,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get minfrac
    cfg_entry = "USAFSI min VIIRS/CDFS-II fraction to mark point as snow:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%minfrac,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get usefrac
    cfg_entry = "USAFSI use CDFS-II fractional snow data:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%usefrac,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get useviirs
    cfg_entry = "USAFSI use VIIRS snow mask:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%useviirs,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! ** New Bratseth settings **
    ! Get obs_err_var
    cfg_entry = "USAFSI observation error variance (m^2):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%ob_err_var,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get back_err_var
    cfg_entry = "USAFSI background error variance (m^2):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, usafsi_settings%back_err_var,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get back_err_h_corr_len
    cfg_entry = "USAFSI background error horizontal correlation length (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%back_err_h_corr_len,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get back_err_v_corr_len
    cfg_entry = "USAFSI background error vertical correlation length (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%back_err_v_corr_len,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get elev_diff_thresh
    cfg_entry = "USAFSI elevQC difference threshold (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%elevqc_diff_threshold,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get skewed_backqc_threshold
    cfg_entry = "USAFSI skewed backQC snow depth threshold (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%skewed_backqc_threshold,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Bogus value for climo if not defined at land point
    cfg_entry = "USAFSI bogus climatology snow depth value (m):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%fill_climo,&
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get gofs_sst_dir
    cfg_entry = "USAFSI GOFS SST data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%gofs_sst_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get gofs_cice_dir
    cfg_entry = "USAFSI GOFS CICE data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%gofs_cice_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get lis_grib2_dir
    cfg_entry = "USAFSI LIS GRIB2 data directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%lis_grib2_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get security class
    cfg_entry = "USAFSI LIS GRIB2 security class:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%security_class, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get data category
    cfg_entry = "USAFSI LIS GRIB2 data category:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%data_category, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get data resolution
    cfg_entry = "USAFSI LIS GRIB2 data resolution:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%data_res, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get area of data
    cfg_entry = "USAFSI LIS GRIB2 area of data:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%area_of_data, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get GALWEM root directory
    cfg_entry = "USAFSI GALWEM root directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%galwem_root_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get GALWEM subdirectory
    cfg_entry = "USAFSI GALWEM subdirectory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%galwem_sub_dir, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! See if directory timestamp used for GALWEM
    cfg_entry = "USAFSI GALWEM use timestamp directories:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%use_timestamp, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get GALWEM nominal resolution
    cfg_entry = "USAFSI GALWEM nominal resolution (km):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%galwem_res, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    ! Get output file name (prefix)
    cfg_entry = "USAFSI netcdf filename prefix:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, &
         usafsi_settings%netcdf_prefix, &
         rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

  end subroutine LDT_usafsiInit

  ! This calls the actual USAFSI driver
  subroutine LDT_usafsiRun(n)
    implicit none
    integer, intent(in) :: n
    external :: USAFSI_run
    call USAFSI_run(n)
  end subroutine LDT_usafsiRun
end module LDT_usafsiMod
