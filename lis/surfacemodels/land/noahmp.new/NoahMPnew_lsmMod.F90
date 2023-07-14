!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!#include "LIS_misc.h"
module NoahMPnew_lsmMod
!BOP
!
! !MODULE: NoahMPnew_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of NoahMPnew model. It also provides the entry method
! for the initialization of NoahMPnew-specific variables. The derived
! data type {\tt NoahmpNew\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the NoahMPnew restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for NoahMPnew
! \item[LDT\_ncvar\_vegetype]
!   LDT NetCDF variable name for vegetation type
! \item[LDT\_ncvar\_soiltype]
!   LDT NetCDF variable name for soil type
! \item[LDT\_ncvar\_tbot]
!   LDT NetCDF variable name for deep soil temperature
! \item[LDT\_ncvar\_planting]
!   LDT NetCDF variable name for planting date
! \item[LDT\_ncvar\_harvest]
!   LDT NetCDF variable name for harvest date
! \item[LDT\_ncvar\_season\_gdd]
!   LDT NetCDF variable name for growing season GDD
! \item[LDT\_ncvar\_soilcL1]
!   LDT NetCDF variable name for soil texture in layer 1
! \item[LDT\_ncvar\_soilcL2]
!   LDT NetCDF variable name for soil texture in layer 2
! \item[LDT\_ncvar\_soilcL3]
!   LDT NetCDF variable name for soil texture in layer 3
! \item[LDT\_ncvar\_soilcL4]
!   LDT NetCDF variable name for soil texture in layer 4
! \item[LDT\_ncvar\_irfract]
!   LDT NetCDF variable name for total irrigation fraction
! \item[LDT\_ncvar\_sifract]
!   LDT NetCDF variable name for sprinkler irrigation fraction
! \item[LDT\_ncvar\_mifract]
!   LDT NetCDF variable name for micro/drip irrigation fraction
! \item[LDT\_ncvar\_fifract]
!   LDT NetCDF variable name for flood irrigation fraction
! \item[LDT\_ncvar\_tdfract]
!   LDT NetCDF variable name for tile drainage fraction
! \item[ts]
!   NoahMPnew model time step in second
! \item[ts\_soil]
!   NoahMPnew model soil time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[noahmpnew]
!  NoahMPnew model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[dz8w]
!   reference height of temperature and humidity
! \item[dt]
!   timestep
! \item[dt_soil]
!   timestep for soil process
! \item[sldpth]
!   thickness of soil layers
! \item[nsoil]
!   number of soil layers
! \item[nsnow]
!   maximum number of snow layers (e.g. 3)
! \item[noahmp\_tbl\_name]
!   NoahMP parameter table
! \item[landuse\_scheme\_name]
!   Landuse classification scheme
! \item[dveg\_opt]
!   dynamic vegetation, 1->off; 2->on); with opt\_crs=1
! \item[crs\_opt]
!   canopt stomatal resistance (1->Ball-Berry; 2->Jarvis)
! \item[btr\_opt]
!   soil moisture factor for stomatal resistance(1->Noah;2->CLM;3->SSiB)
! \item[runsfc\_opt]
!   surface runoff (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XinAnJiang; 8->DynamicVIC)
! \item[runsub\_opt]
!   subsurface runoff (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS; 5->MMF; 6->VIC; 7->XinAnJiang; 8->DynamicVIC)
! \item[sfc\_opt]
!   surface layer drag coeff (CH \& CM) (1->M-O; 2->Chen97)
! \item[frz\_opt]
!   supercooled liquid water (1->NY06; 2->Koren99)
! \item[tksno\_opt]
!   snow thermal conductivity (1->Yen1965; 2->Anderson1976; 3->constant; 4->Verseghy1991; 5->Yen1981)
! \item[inf\_opt]
!   frozen soil permeability (1->NY06; 2->Koren99)
! \item[rad\_opt]
!   radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
! \item[alb\_opt]
!   snow surface albedo (1->BATS; 2->CLASS)
! \item[snf\_opt]
!   rainfall \& snowfall (1->Jordan91; 2->BATS; 3->Noah)
! \item[tbot\_opt]
!   lower boundary of soil temperature
! \item[stc\_opt]
!   snow/soil temperature time scheme
! \item[gla\_opt]
!   glacier option (1->phase change; 2->simple)
! \item[rsf\_opt]
!   surface resistance(1->Sakaguchi/Zeng;2->Seller;3->mod Sellers;4->1+snow)
! \item[soil\_opt]
!   soil configuration option
! \item[pedo\_opt]
!   soil pedotransfer function option
! \item[crop\_opt]
!   crop model option (0->none; 1->Liu2016)
! \item[irr\_opt]
!    irrigation scheme option (0->none; 1->always on; 2->trigger by planting/harvest dates; 3->trigger by LAI)
! \item[irrm\_opt]
!    irrigation method option (0->fraction from input; 1->sprinkler; 2->micro/drip; 3->flood)
! \item[tdrn\_opt]
!    tile drainage option (0->none; 1->simple drainage; 2->Hooghoudt's scheme)
! \item[urban\_opt]
!   urban physics option (0->none; 1->SLUCM; 2->BEP; 3->BEP_BEM)
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  10/25/18: Shugong Wang, Zhuo Wang Initial implementation for LIS 7 and NoahMP401
!  05/01/23: Cenlin He, update to work with refactored Noah-MP (v5.0 and later)

! !USES:
    use NoahMPnew_module
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: NoahMPnew_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: NoahmpNew_struc
!EOP
    type, public :: NoahMPnew_type_dec
        character(len=LIS_CONST_PATH_LEN) :: rfile
        character*256      :: rformat
        !-------------------------------------------------------------------------
        ! Parameter file names
        !-------------------------------------------------------------------------
        character*128      :: LDT_ncvar_vegetype
        character*128      :: LDT_ncvar_soiltype
        character*128      :: LDT_ncvar_shdfac_monthly
        character*128      :: LDT_ncvar_tbot
        character*128      :: LDT_ncvar_planting
        character*128      :: LDT_ncvar_harvest
        character*128      :: LDT_ncvar_season_gdd
        character*128      :: LDT_ncvar_soilcomp
        character*128      :: LDT_ncvar_soilcL1
        character*128      :: LDT_ncvar_soilcL2
        character*128      :: LDT_ncvar_soilcL3
        character*128      :: LDT_ncvar_soilcL4
        character*128      :: LDT_ncvar_irfract
        character*128      :: LDT_ncvar_sifract
        character*128      :: LDT_ncvar_mifract
        character*128      :: LDT_ncvar_fifract
        character*128      :: LDT_ncvar_tdfract
        character*128      :: LDT_ncvar_fdepth
        character*128      :: LDT_ncvar_eqzwt
        character*128      :: LDT_ncvar_rechclim
        character*128      :: LDT_ncvar_riverbed

        !-------------------------------------------------------------------------
        ! ts, Count, rstInterval, outInterval
        !-------------------------------------------------------------------------
        real               :: ts
        real               :: ts_soil
        integer            :: count
        real               :: rstInterval
        integer            :: outInterval
        integer            :: forc_count
        !-------------------------------------------------------------------------
        ! Initial Model State for cold start
        !-------------------------------------------------------------------------
        real               :: init_tskin
        real               :: init_sneqv
        real               :: init_snowh
        real               :: init_canwat
        real, pointer      :: init_tslb(:)
        real, pointer      :: init_smc(:)
        real               :: init_zwt
        real               :: init_wa
        real               :: init_wt
        real               :: init_lai
        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        real               :: dz8w
        real               :: dt
        real               :: dt_soil
        real, pointer      :: sldpth(:)
        integer            :: nsoil
        integer            :: nsnow
        character(len=256) :: noahmp_tbl_name
        character(len=256) :: landuse_scheme_name
        integer            :: dveg_opt
        integer            :: crs_opt
        integer            :: btr_opt
        integer            :: runsfc_opt
        integer            :: runsub_opt
        integer            :: sfc_opt
        integer            :: frz_opt
        integer            :: inf_opt
        integer            :: rad_opt
        integer            :: alb_opt
        integer            :: snf_opt
        integer            :: tksno_opt
        integer            :: tbot_opt
        integer            :: stc_opt
        integer            :: gla_opt
        integer            :: sndpth_gla_opt
        integer            :: rsf_opt
        integer            :: soil_opt
        integer            :: pedo_opt
        integer            :: crop_opt
        integer            :: irr_opt
        integer            :: irrm_opt
        integer            :: tdrn_opt
        integer            :: infdv_opt
        integer            :: urban_opt
        type(NoahMPnewdec), pointer :: noahmpnew(:)
    end type NoahMPnew_type_dec

    type(NoahMPnew_type_dec), pointer :: NoahmpNew_struc(:)
 
contains 

!BOP
!
! !ROUTINE: NoahMPnew_ini
! \label{NoahMPnew_ini}
!
! !INTERFACE:
    subroutine NoahMPnew_ini()
! !USES:
        use LIS_coreMod, only : LIS_rc
        use LIS_logMod, only : LIS_verify, LIS_logunit
        use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
            LIS_update_timestep, LIS_registerAlarm
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for NoahMPnew-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for NoahMPnew from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[NoahMPnew\_readcrd](\ref{NoahMPnew_readcrd})\\
!    reads the runtime options for NoahMPnew model
!  \end{description}
!EOP
        implicit none        
        integer  :: n, t     
        integer  :: status   
        character*3 :: fnest ! EMK for RHMin

        ! allocate memory for nest 
        allocate(NoahmpNew_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call NoahMPnew_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(NoahmpNew_struc(n)%noahmpnew(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            ! allocate memory for multilevel spatial parameter
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%shdfac_monthly(12))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%soilcomp(8))
             enddo

            ! allocate memory for state variables
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%smc(NoahmpNew_struc(n)%nsoil))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%sh2o(NoahmpNew_struc(n)%nsoil))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%tslb(NoahmpNew_struc(n)%nsoil))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%tsno(NoahmpNew_struc(n)%nsnow))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%zss(NoahmpNew_struc(n)%nsnow + NoahmpNew_struc(n)%nsoil))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%snowice(NoahmpNew_struc(n)%nsnow))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%snowliq(NoahmpNew_struc(n)%nsnow))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%smoiseq(NoahmpNew_struc(n)%nsoil))
                allocate(NoahmpNew_struc(n)%noahmpnew(t)%accetrani(NoahmpNew_struc(n)%nsoil))
            enddo

            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                NoahmpNew_struc(n)%noahmpnew(t)%lwdown    = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%swdown    = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%psurf     = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%prcp      = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%tair      = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%qair      = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%wind_e    = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%wind_n    = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%sfcheadrt = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irnumsi   = 0
                NoahmpNew_struc(n)%noahmpnew(t)%irnummi   = 0
                NoahmpNew_struc(n)%noahmpnew(t)%irnumfi   = 0
                NoahmpNew_struc(n)%noahmpnew(t)%irwatsi   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irwatmi   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irwatfi   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irsivol   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irmivol   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irfivol   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%ireloss   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irrsplh   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%qtdrain   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accssoil  = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accqinsur = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accqseva  = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accetrani = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accdwater = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accprcp   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accecan   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accetran  = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%accedir   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%sfcrunoff = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%udrrunoff = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%deeprech  = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%rech      = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%acsnom    = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%acsnow    = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%irfract   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%sifract   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%mifract   = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%fifract   = 0.0
                !ag(05Jan2021)
                NoahmpNew_struc(n)%noahmpnew(t)%rivsto = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%fldsto = 0.0
                NoahmpNew_struc(n)%noahmpnew(t)%fldfrc = 0.0
            enddo ! end of tile (t) loop

            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            NoahmpNew_struc(n)%forc_count = 0

            call LIS_update_timestep(LIS_rc, n, NoahmpNew_struc(n)%ts)

            call LIS_registerAlarm("NoahMPnew model alarm",&
                                   NoahmpNew_struc(n)%ts, &
                                   NoahmpNew_struc(n)%ts)

            ! CH2023: add soil timestep that is allowed to be different from main timestep
            call LIS_registerAlarm("NoahMPnew model alarm",&
                                   NoahmpNew_struc(n)%ts, &
                                   NoahmpNew_struc(n)%ts_soil)

            call LIS_registerAlarm("NoahMPnew restart alarm", &
                                   NoahmpNew_struc(n)%ts,&
                                   NoahmpNew_struc(n)%rstInterval)

            ! EMK Add alarm to reset tair_agl_min for RHMin.  This should 
            ! match the output interval, since that is used for calculating 
            ! Tair_F_min.            
            write(fnest,'(i3.3)') n
            call LIS_registerAlarm("NoahMPnew RHMin alarm "//trim(fnest),&
                 NoahmpNew_struc(n)%ts,&
                 LIS_sfmodel_struc(n)%outInterval)
            if (LIS_sfmodel_struc(n)%outInterval .gt. 86400 .or. &
                 trim(LIS_sfmodel_struc(n)%outIntervalType) .eq. "dekad") then
               write(LIS_logunit,*) &
                    '[WARN] If RHMin is selected for output, please reset ', &
                    'surface model output interval to no more than 24 hours.'
            end if

            ! Initialize min/max values to implausible values.
            NoahmpNew_struc(n)%noahmpnew(:)%tair_agl_min = 999.0
            NoahmpNew_struc(n)%noahmpnew(:)%rhmin = 999.0

            !------------------------------------------------------------------------
            ! TODO: setup number of soil moisture/temperature layers and depth here  
            !------------------------------------------------------------------------
            ! TODO: set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = NoahmpNew_struc(n)%nsoil
            ! TODO: set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = NoahmpNew_struc(n)%nsoil
            allocate(LIS_sfmodel_struc(n)%lyrthk(NoahmpNew_struc(n)%nsoil))
            !LIS_sfmodel_struc(n)%lyrthk(:) = NoahmpNew_struc(n)%sldpth(:)
            !EMK...Output soil layer thicknesses in centimeters for 
            !consistency with other LSMs.  
            LIS_sfmodel_struc(n)%lyrthk(:) = &
                 100*NoahmpNew_struc(n)%sldpth(:)
            LIS_sfmodel_struc(n)%ts = NoahmpNew_struc(n)%ts
        enddo
    end subroutine NoahMPnew_ini
end module NoahMPnew_lsmMod
