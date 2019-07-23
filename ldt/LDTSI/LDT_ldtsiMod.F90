!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_ldtsiMod

   ! Defaults
   implicit none
   private

   ! Public methods
   public :: LDT_ldtsiInit
   public :: LDT_ldtsiRun

   ! Public type.  These are formerly environment and namelist variables from 
   ! the original SNODEP.
   type, public :: ldtsi_t
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

   end type ldtsi_t
   type(ldtsi_t), public :: ldtsi_settings

contains

   ! Reads LDTSI-specific entries from ldt.config
   subroutine LDT_ldtsiInit()

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
      cfg_entry = "LDTSI valid date (YYYYMMDDHH):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%date10, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get fracdir      
      cfg_entry = "LDTSI fractional snow data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%fracdir, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      
      ! Get modif
      cfg_entry = "LDTSI modified data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%modif, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get sfcobs
      cfg_entry = "LDTSI surface obs data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%sfcobs, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get ssmis.  
      cfg_entry = "LDTSI SSMIS data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%ssmis, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get stmpdir
      cfg_entry = "LDTSI surface temperature data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%stmpdir, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get static
      cfg_entry = "LDTSI static data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%static, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get unmod
      cfg_entry = "LDTSI unmodified data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%unmod, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get viirsdir
      cfg_entry = "LDTSI VIIRS data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%viirsdir, rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      
      ! *** for SSMIS now depth retrieval, Yeosang Yoon
      ! get SSMIS raw datasets
      cfg_entry = "LDTSI SSMIS raw data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%ssmis_raw_dir, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! get option for SSMIS snow depth retrieval alogrithm
      cfg_entry = "LDTSI SSMIS snow depth retrieval algorithm option:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%ssmis_option, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! get forest fraction for algorithm 3
      if (ldtsi_settings%ssmis_option==3) then         
         cfg_entry = "LDTSI SSMIS forest fraction file:"         
         call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
         call LDT_verify(rc, trim(cfg_entry)//" not specified")
         call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%ff_file, &
              rc=rc)
         call LDT_verify(rc, trim(cfg_entry)//" not specified")
      end if

      ! *** Get former namelist variables ***

      ! Get clmadj
      cfg_entry = &
           "LDTSI decimal fraction adjustment of snow depth towards climo:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%clmadj,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get unkdep
      cfg_entry = "LDTSI default snow depth (m) when actual depth unknown:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%unkdep,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      
      ! Get minprt
      cfg_entry = &
           "LDTSI minimum snow depth (m) for which to print a diagnostic:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%minprt,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get maxsobs
      cfg_entry = "LDTSI maximum number of surface observations allowed:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%maxsobs,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      
      ! Get minsat
      cfg_entry = "LDTSI SSMIS shallow snow depth threshold (m):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%minsat,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get trplat
      cfg_entry = "LDTSI latitudes (deg * 100) for summer climo check:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      do c = 1,3
         call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%trplat(c),&
              rc=rc)
         call LDT_verify(rc, trim(cfg_entry)//" not specified")
      end do ! c

      ! Get elvlim
      cfg_entry = "LDTSI elevations (m) for summer climo check:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      do c = 1,4
         call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%elvlim(c),&
              rc=rc)
         call LDT_verify(rc, trim(cfg_entry)//" not specified")
      end do ! c
      
      ! Get thresh
      cfg_entry = &
           "LDTSI temperature (deg K * 10) above which no snow is allowed:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%thresh,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")


      ! Get arctmax
      cfg_entry = &
           "LDTSI max reported temperature (deg K * 10) allowed around poles:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%arctmax,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      
      ! Get minice
      cfg_entry = "LDTSI minimum ice concentration (%) needed to set ice flag:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%minice,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      

      ! Get icelat
      cfg_entry = "LDTSI high latitude thresholds (deg) for sea ice::"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      do r = 1,2
         call ESMF_ConfigNextLine(LDT_config, rc=rc)
         call LDT_verify(rc, trim(cfg_entry)//" not specified")      
         do c = 1, 12
            call ESMF_ConfigGetAttribute(LDT_config, &
                 ldtsi_settings%icelat(c,r), rc=rc)
            call LDT_verify(rc, trim(cfg_entry)//" not specified")      
         end do ! c
      end do ! r

      ! Get latchk
      cfg_entry = "LDTSI low latitude thresholds (deg) for sea ice::"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      do r = 1,2
         call ESMF_ConfigNextLine(LDT_config, rc=rc)
         call LDT_verify(rc, trim(cfg_entry)//" not specified")      
         do c = 1, 12
            call ESMF_ConfigGetAttribute(LDT_config, &
                 ldtsi_settings%latchk(c,r), rc=rc)
            call LDT_verify(rc, trim(cfg_entry)//" not specified")      
         end do ! c
      end do ! r

      ! Get maxpixage
      cfg_entry = "LDTSI max age of VIIRS pixels to use:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%maxpixage,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      

      ! Get minbare
      cfg_entry = "LDTSI min VIIRS fraction to mark point as bare ground:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%minbare,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      
      ! Get minfrac
      cfg_entry = "LDTSI min VIIRS/CDFS-II fraction to mark point as snow:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%minfrac,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      

      ! Get usefrac
      cfg_entry = "LDTSI use CDFS-II fractional snow data:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%usefrac,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      

      ! Get useviirs
      cfg_entry = "LDTSI use VIIRS snow mask:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%useviirs,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      

      ! ** New Bratseth settings **
      ! Get obs_err_var
      cfg_entry = "LDTSI observation error variance (m^2):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%ob_err_var,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")      

      ! Get back_err_var
      cfg_entry = "LDTSI background error variance (m^2):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, ldtsi_settings%back_err_var,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get back_err_h_corr_len
      cfg_entry = "LDTSI background error horizontal correlation length (m):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%back_err_h_corr_len,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get back_err_v_corr_len
      cfg_entry = "LDTSI background error vertical correlation length (m):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%back_err_v_corr_len,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get elev_diff_thresh
      cfg_entry = "LDTSI elevQC difference threshold (m):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%elevqc_diff_threshold,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get skewed_backqc_threshold
      cfg_entry = "LDTSI skewed backQC snow depth threshold (m):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%skewed_backqc_threshold,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Bogus value for climo if not defined at land point
      cfg_entry = "LDTSI bogus climatology snow depth value (m):"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%fill_climo,&
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get gofs_sst_dir
      cfg_entry = "LDTSI GOFS SST data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%gofs_sst_dir, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get gofs_cice_dir
      cfg_entry = "LDTSI GOFS CICE data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%gofs_cice_dir, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get lis_grib2_dir
      cfg_entry = "LDTSI LIS GRIB2 data directory:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%lis_grib2_dir, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get security class
      cfg_entry = "LDTSI LIS GRIB2 security class:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%security_class, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get data category
      cfg_entry = "LDTSI LIS GRIB2 data category:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%data_category, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get data resolution
      cfg_entry = "LDTSI LIS GRIB2 data resolution:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%data_res, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")

      ! Get area of data
      cfg_entry = "LDTSI LIS GRIB2 area of data:"
      call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
      call ESMF_ConfigGetAttribute(LDT_config, &
           ldtsi_settings%area_of_data, &
           rc=rc)
      call LDT_verify(rc, trim(cfg_entry)//" not specified")
     
   end subroutine LDT_ldtsiInit

   ! This calls the actual LDTSI driver
   subroutine LDT_ldtsiRun()
      implicit none
      call LDTSI_run()
   end subroutine LDT_ldtsiRun
end module LDT_ldtsiMod
