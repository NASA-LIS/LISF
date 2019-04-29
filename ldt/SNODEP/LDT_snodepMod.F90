!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_snodepMod

   ! Defaults
   implicit none
   private

   ! Public methods
   public :: LDT_snodepInit
   public :: LDT_snodepRun

   ! Public type.  These are formerly environment and namelist variables from 
   ! the original SNODEP.
   type, public :: snodep_t
      ! Former environment variables
      character*10  :: date10      
      character*100 :: fracdir
      character*100 :: modif
      character*100 :: sfcobs
      character*100 :: ssmis
      character*100 :: stmpdir
      character*100 :: static
      character*100 :: unmod
      character*100 :: viirsdir
      
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

      ! option for SSMIS snow depth retrieval algorithm
      integer       :: ssmis_option
      character*100 :: ssmis_raw_dir
      character*100 :: ff_file

      ! Bratseth settings
      real :: ob_err_var
      real :: back_err_var
      real :: back_err_h_corr_len
      real :: back_err_v_corr_len
      real :: elevqc_diff_threshold
      real :: skewed_backqc_threshold

      ! Other new settings
      real :: fill_climo
      character*100 :: gofs_sst_dir
      character*100 :: gofs_cice_dir
      character*100 :: lis_grib2_dir
      character*20 :: security_class
      character*20 :: data_category
      character*20 :: data_res
      character*20 :: area_of_data

   end type snodep_t
   type(snodep_t), public :: snodep_settings

contains

   ! Reads SNODEP-specific entries from ldt.config
   subroutine LDT_snodepInit()

      ! Imports
      use ESMF
      use LDT_coreMod, only: LDT_config
      use LDT_logMod, only: LDT_verify

      ! Defaults
      implicit none

      ! Local variables
      integer :: rc
      integer :: c,r

      ! *** Get former environment variables ***

      ! Get date10
      call ESMF_ConfigFindLabel(LDT_config,"SNODEP valid date (YYYYMMDDHH):",&
           rc=rc)
      call LDT_verify(rc, "SNODEP valid date (YYYYMMDDHH): not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%date10,rc=rc)
      call LDT_verify(rc, "SNODEP valid date (YYYYMMDDHH): not specified")

      ! Get fracdir
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP fractional snow data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP fractional snow data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%fracdir,rc=rc)
      call LDT_verify(rc, &
           "SNODEP fractional snow data directory: not specified")
      
      ! Get modif
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP modified data directory:",rc=rc)
      call LDT_verify(rc, "SNODEP modified data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%modif,rc=rc)
      call LDT_verify(rc, "SNODEP modified data directory: not specified")

      ! Get sfcobs
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP surface obs data directory:",rc=rc)
      call LDT_verify(rc, "SNODEP surface obs data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%sfcobs,rc=rc)
      call LDT_verify(rc, "SNODEP surface obs data directory: not specified")

      ! Get ssmis.  
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP SSMIS data directory:",rc=rc)
      call LDT_verify(rc, "SNODEP SSMIS data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%ssmis,rc=rc)
      call LDT_verify(rc, "SNODEP SSMIS data directory: not specified")

      ! Get stmpdir
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP surface temperature data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP surface temperature data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%stmpdir,rc=rc)
      call LDT_verify(rc, &
           "SNODEP surface temperature data directory: not specified")

      ! Get static
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP static data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP static data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%static,rc=rc)
      call LDT_verify(rc, &
           "SNODEP static data directory: not specified")

      ! Get unmod
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP unmodified data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP unmodified data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%unmod,rc=rc)
      call LDT_verify(rc, &
           "SNODEP unmodified data directory: not specified")

      ! Get viirsdir
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP VIIRS data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP VIIRS data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%viirsdir,rc=rc)
      call LDT_verify(rc, &
           "SNODEP VIIRS data directory: not specified")
      
      ! *** for SSMIS now depth retrieval, Yeosang Yoon
      ! get SSMIS raw datasets
      call ESMF_ConfigFindLabel(LDT_config, &
           "SSMIS raw data directory:", rc=rc)
      call LDT_verify(rc, &
           "SSMIS raw data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%ssmis_raw_dir, rc=rc)
      call LDT_verify(rc, &
           "SSMIS raw data directory: not specified") 

      ! get option for SSMIS snow depth retrieval alogrithm
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP SSMIS snow depth retrieval algorithm option:", rc=rc)
      call LDT_verify(rc, &
           "SNODEP SSMIS snow depth retrieval algorithm option: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%ssmis_option, rc=rc)
      call LDT_verify(rc, &
           "SNODEP SSMIS snow depth retrieval algorithm option: not specified")

      ! get forest fraction for algorithm 3
      if (snodep_settings%ssmis_option==3) then

        call ESMF_ConfigFindLabel(LDT_config, "SSMIS forest fraction file:", rc=rc)
        call LDT_verify(rc, "SSMIS forest fraction file: not specified")
        call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%ff_file, rc=rc)
        call LDT_verify(rc, "SSMIS forest fraction file: not specified")
      end if 


      ! *** Get former namelist variables ***

      ! Get clmadj
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP decimal fraction adjustment of snow depth towards climo:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP decimal fraction adjustment of snow depth towards climo: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%clmadj,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP decimal fraction adjustment of snow depth towards climo: not specified")

      ! Get unkdep
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP default snow depth (m) when actual depth unknown:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP default snow depth (m) when actual depth unknown: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%unkdep,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP default snow depth (m) when actual depth unknown: not specified")
      
      ! Get minprt
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP minimum snow depth (m) for which to print a diagnostic:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP minimum snow depth (m) for which to print a diagnostic: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%minprt,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP minimum snow depth (m) for which to print a diagnostic: not specified")

      ! Get maxsobs
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP maximum number of surface observations allowed:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP maximum number of surface observations allowed: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%maxsobs,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP maximum number of surface observations allowed: not specified")
      
      ! Get minsat
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP SSMIS shallow snow depth threshold (m):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP SSMIS shallow snow depth threshold (m): not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%minsat,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP SSMIS shallow snow depth threshold (m): not specified")

      ! Get trplat
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP latitudes (deg * 100) for summer climo check:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP latitudes (deg * 100) for summer climo check: not specified")
      do c = 1,3
         call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%trplat(c),&
              rc=rc)
         call LDT_verify(rc, &
              "SNODEP latitudes (deg * 100) for summer climo check: not specified")
      end do ! c

      ! Get elvlim
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP elevations (m) for summer climo check:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP elevations (m) for summer climo check: not specified")
      do c = 1,4
         call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%elvlim(c),&
              rc=rc)
         call LDT_verify(rc, &
              "SNODEP elevations (m) for summer climo check: not specified")
      end do ! c
      
      ! Get thresh
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP temperature (deg K * 10) above which no snow is allowed:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP temperature (deg K * 10) above which no snow is allowed: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%thresh,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP temperature (deg K * 10) above which no snow is allowed: not specified")

      ! Get arctmax
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP max reported temperature (deg K * 10) allowed around poles:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP max reported temperature (deg K * 10) allowed around poles: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%arctmax,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP max reported temperature (deg K * 10) allowed around poles: not specified")
      
      ! Get minice
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP minimum ice concentration (%) needed to set ice flag:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP minimum ice concentration (%) needed to set ice flag: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%minice,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP minimum ice concentration (%) needed to set ice flag: not specified")
      
      ! Get icelat
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP high latitude thresholds (deg) for sea ice::", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP high latitude thresholds (deg) for sea ice:: not specified")
      do r = 1,2
         call ESMF_ConfigNextLine(LDT_config,rc=rc)
         call LDT_verify(rc, &
              "SNODEP high latitude thresholds (deg) for sea ice:: not specified")
         do c = 1, 12
            call ESMF_ConfigGetAttribute(LDT_config, &
                 snodep_settings%icelat(c,r),rc=rc)
            call LDT_verify(rc, &
                 "SNODEP high latitude thresholds (deg) for sea ice:: not specified")
         end do ! c
      end do ! r

      ! Get latchk
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP low latitude thresholds (deg) for sea ice::", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP low latitude thresholds (deg) for sea ice:: not specified")
      do r = 1,2
         call ESMF_ConfigNextLine(LDT_config,rc=rc)
         call LDT_verify(rc, &
              "SNODEP low latitude thresholds (deg) for sea ice:: not specified")
         do c = 1, 12
            call ESMF_ConfigGetAttribute(LDT_config, &
                 snodep_settings%latchk(c,r), rc=rc)
            call LDT_verify(rc, &
                 "SNODEP low latitude thresholds (deg) for sea ice:: not specified")
         end do ! c
      end do ! r

      ! Get maxpixage
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP max age of VIIRS pixels to use:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP max age of VIIRS pixels to use: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%maxpixage,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP max age of VIIRS pixels to use: not specified")
      
      ! Get minbare
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP min VIIRS fraction to mark point as bare ground:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP min VIIRS fraction to mark point as bare ground: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%minbare,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP min VIIRS fraction to mark point as bare ground: not specified")
      
      ! Get minfrac
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP min VIIRS/CDFS-II fraction to mark point as snow:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP min VIIRS/CDFS-II fraction to mark point as snow: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%minfrac,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP min VIIRS/CDFS-II fraction to mark point as snow: not specified")

      ! Get usefrac
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP use CDFS-II fractional snow data:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP use CDFS-II fractional snow data: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%usefrac,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP use CDFS-II fractional snow data: not specified")

      ! Get useviirs
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP use VIIRS snow mask:", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP use VIIRS snow mask: not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%useviirs,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP use VIIRS snow mask: not specified")

      ! ** New Bratseth settings **
      ! Get obs_err_var
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP observation error variance (m^-2):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP observation error variance (m^-2): not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%ob_err_var,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP observation error variance (m^-2): not specified")

      ! Get back_err_var
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP background error variance (m^-2):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP background error variance (m^-2): not specified")
      call ESMF_ConfigGetAttribute(LDT_config,snodep_settings%back_err_var,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP background error variance (m^-2): not specified")

      ! Get back_err_h_corr_len
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP background error horizontal correlation length (m):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP background error horizontal correlation length (m):"// &
           "not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%back_err_h_corr_len,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP background error horizontal correlation length (m):"// &
           "not specified")

      ! Get back_err_v_corr_len
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP background error vertical correlation length (m):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP background error vertical correlation length (m):"// &
           "not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%back_err_v_corr_len,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP background error vertical correlation length (m):"// &
           "not specified")

      ! Get elev_diff_thresh
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP elevQC difference threshold (m):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP elevQC difference threshold  (m):"// &
           "not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%elevqc_diff_threshold,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP elevQC difference threshold (m):"// &
           "not specified")

      ! Get skewed_backqc_threshold
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP skewed backQC snow depth threshold (m):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP skewed backQC snow depth threshold (m):"// &
           "not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%skewed_backqc_threshold,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP skewed backQC snow depth threshold (m):"// &
           "not specified")

      ! Bogus value for climo if not defined at land point
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP bogus climatology snow depth value (m):", &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP bogus climatology snow depth value (m):"// &
           "not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%fill_climo,&
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP bogus climatology snow depth value (m):"// &
           "not specified")

      ! Get gofs_sst_dir
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP GOFS SST data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP GOFS SST data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%gofs_sst_dir, &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP GOFS SST data directory: not specified")

      ! Get gofs_cice_dir
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP GOFS CICE data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP GOFS CICE data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%gofs_cice_dir, &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP GOFS CICE data directory: not specified")

      ! Get lis_grib2_dir
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP LIS GRIB2 data directory:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 data directory: not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%lis_grib2_dir, &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 data directory: not specified")

      ! Get security class
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP LIS GRIB2 security class:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 security class: not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%security_class, &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 security class: not specified")

      ! Get data category
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP LIS GRIB2 data category:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 data category: not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%data_category, &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 data category: not specified")

      ! Get data resolution
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP LIS GRIB2 data resolution:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 data resolution: not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%data_res, &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 data resolution: not specified")

      ! Get area of data
      call ESMF_ConfigFindLabel(LDT_config, &
           "SNODEP LIS GRIB2 area of data:",rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 area of data: not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           snodep_settings%area_of_data, &
           rc=rc)
      call LDT_verify(rc, &
           "SNODEP LIS GRIB2 area of data: not specified")

      

   end subroutine LDT_snodepInit

   ! This calls the actual SNODEP driver
   subroutine LDT_snodepRun()
      implicit none
      call SNODEP_run()
   end subroutine LDT_snodepRun
end module LDT_snodepMod
