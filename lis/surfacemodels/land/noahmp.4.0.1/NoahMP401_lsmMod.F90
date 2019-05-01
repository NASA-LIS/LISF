!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0     
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
!#include "LIS_misc.h"
module NoahMP401_lsmMod
!BOP
!
! !MODULE: NoahMP401_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of NoahMP401 model. It also provides the entry method
! for the initialization of NoahMP401-specific variables. The derived
! data type {\tt NoahMP401\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the NoahMP401 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for NoahMP401
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
! \item[ts]
!   NoahMP401 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[noahmp401]
!  NoahMP401 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[dz8w]
!   reference height of temperature and humidity
! \item[dt]
!   timestep
! \item[sldpth]
!   thickness of soil layers
! \item[nsoil]
!   number of soil layers
! \item[nsnow]
!   maximum number of snow layers (e.g. 3)
! \item[landuse\_tbl\_name]
!   Noah model landuse parameter table
! \item[soil\_tbl\_name]
!   Noah model soil parameter table
! \item[gen\_tbl\_name]
!   Noah model general parameter table
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
! \item[run\_opt]
!   runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
! \item[sfc\_opt]
!   surface layer drag coeff (CH \& CM) (1->M-O; 2->Chen97)
! \item[frz\_opt]
!   supercooled liquid water (1->NY06; 2->Koren99)
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
!   crop model option (0->none; 1->Liu et al.; 2->Gecros)
! \item[urban\_opt]
!   urban physics option
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  10/25/18: Shugong Wang, Zhuo Wang Initial implementation for LIS 7 and NoahMP401
!
! !USES:
    use NoahMP401_module

    implicit none

    PRIVATE
    !-------------------------------------------------------------------------
    ! PUBLIC MEMBER FUNCTIONS
    !-------------------------------------------------------------------------
    public :: NoahMP401_ini
    !-------------------------------------------------------------------------
    ! PUBLIC TYPES
    !-------------------------------------------------------------------------
    public :: NoahMP401_struc
!EOP
    type, public :: NoahMP401_type_dec
        character*256      :: rfile
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
        !-------------------------------------------------------------------------
        ! ts, Count, rstInterval, outInterval
        !-------------------------------------------------------------------------
        real               :: ts
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
        real, pointer      :: init_gecros_state(:)
        !-------------------------------------------------------------------------
        ! Constant Parameter
        !-------------------------------------------------------------------------
        real               :: dz8w
        real               :: dt
        real, pointer      :: sldpth(:)
        integer            :: nsoil
        integer            :: nsnow
        character(len=256) :: soil_tbl_name
        character(len=256) :: gen_tbl_name
        character(len=256) :: noahmp_tbl_name
        character(len=256) :: landuse_scheme_name
        integer            :: dveg_opt
        integer            :: crs_opt
        integer            :: btr_opt
        integer            :: run_opt
        integer            :: sfc_opt
        integer            :: frz_opt
        integer            :: inf_opt
        integer            :: rad_opt
        integer            :: alb_opt
        integer            :: snf_opt
        integer            :: tbot_opt
        integer            :: stc_opt
        integer            :: gla_opt
        integer            :: rsf_opt
        integer            :: soil_opt
        integer            :: pedo_opt
        integer            :: crop_opt
        integer            :: urban_opt
        type(NoahMP401dec), pointer :: noahmp401(:)
    end type NoahMP401_type_dec

    type(NoahMP401_type_dec), pointer :: NOAHMP401_struc(:)
 
contains 

!BOP
!
! !ROUTINE: NoahMP401_ini
! \label{NoahMP401_ini}
!
! !INTERFACE:
    subroutine NoahMP401_ini()
! !USES:
        use LIS_coreMod, only : LIS_rc
        use LIS_logMod, only : LIS_verify
        use LIS_timeMgrMod, only : LIS_clock,  LIS_calendar, &
            LIS_update_timestep, LIS_registerAlarm
        use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for NoahMP401-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for NoahMP401 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[NoahMP401\_readcrd](\ref{NoahMP401_readcrd})\\
!    reads the runtime options for NoahMP401 model
!  \end{description}
!EOP
        implicit none        
        integer  :: n, t     
        integer  :: status   

        ! allocate memory for nest 
        allocate(NOAHMP401_struc(LIS_rc%nnest))
 
        ! read configuation information from lis.config file
        call NoahMP401_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(NOAHMP401_struc(n)%noahmp401(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            ! allocate memory for multilevel spatial parameter
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(NOAHMP401_struc(n)%noahmp401(t)%shdfac_monthly(12))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%soilcomp(8))
             enddo

            ! allocate memory for state variables
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(NOAHMP401_struc(n)%noahmp401(t)%smc(NOAHMP401_struc(n)%nsoil))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%sh2o(NOAHMP401_struc(n)%nsoil))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%tslb(NOAHMP401_struc(n)%nsoil))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%tsno(NOAHMP401_struc(n)%nsnow))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%zss( NOAHMP401_struc(n)%nsnow + NOAHMP401_struc(n)%nsoil))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%snowice(NOAHMP401_struc(n)%nsnow))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%snowliq(NOAHMP401_struc(n)%nsnow))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%smoiseq(NOAHMP401_struc(n)%nsoil))
                allocate(NOAHMP401_struc(n)%noahmp401(t)%gecros_state(60))
            enddo

            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                NOAHMP401_struc(n)%noahmp401(t)%lwdown = 0.0
                NOAHMP401_struc(n)%noahmp401(t)%swdown = 0.0
                NOAHMP401_struc(n)%noahmp401(t)%psurf = 0.0
                NOAHMP401_struc(n)%noahmp401(t)%prcp = 0.0
                NOAHMP401_struc(n)%noahmp401(t)%tair = 0.0
                NOAHMP401_struc(n)%noahmp401(t)%qair = 0.0
                NOAHMP401_struc(n)%noahmp401(t)%wind_e = 0.0
                NOAHMP401_struc(n)%noahmp401(t)%wind_n = 0.0
            enddo ! end of tile (t) loop

            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            NOAHMP401_struc(n)%forc_count = 0

            call LIS_update_timestep(LIS_rc, n, NOAHMP401_struc(n)%ts)

            call LIS_registerAlarm("NoahMP401 model alarm",&
                                   NOAHMP401_struc(n)%ts, &
                                   NOAHMP401_struc(n)%ts)

            call LIS_registerAlarm("NoahMP401 restart alarm", &
                                   NOAHMP401_struc(n)%ts,&
                                   NOAHMP401_struc(n)%rstInterval)
            !------------------------------------------------------------------------
            ! TODO: setup number of soil moisture/temperature layers and depth here  
            !------------------------------------------------------------------------
            ! TODO: set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = NOAHMP401_struc(n)%nsoil
            ! TODO: set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = NOAHMP401_struc(n)%nsoil
            allocate(LIS_sfmodel_struc(n)%lyrthk(NOAHMP401_struc(n)%nsoil))
            LIS_sfmodel_struc(n)%lyrthk(:) = NOAHMP401_struc(n)%sldpth(:)
            LIS_sfmodel_struc(n)%ts = NOAHMP401_struc(n)%ts
        enddo
    end subroutine NoahMP401_ini
end module NoahMP401_lsmMod
