!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: Crocus81_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to
! control the operation of Crocus81 model. It also provides the entry method
! for the initialization of Crocus81-specific variables. The derived
! data type {\tt Crocus81\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
! \item[rfile]
!  name of the Crocus81 restart file
! \item[rformat]
!  format of restart file (binary or netcdf) for Crocus81
! \item[LDT\_ncvar\_GLACIER\_BOOL]
!   LDT NetCDF variable name for !True = Over permanent snow and ice, initialise WGI=WSAT, Hsnow$>$=10m and allow 0.8$<$SNOALB$<$0.85
! False = No specific treatment
! \item[LDT\_ncvar\_SLOPE]
!   LDT NetCDF variable name for angle between the normal to the surface and the vertical  (MN:  replaced PDIRCOSZW with slope and computed the cosine in the driver)
! \item[LDT\_ncvar\_SOILCOND]
!   LDT NetCDF variable name for soil thermal conductivity (W m-1 K-1)
! \item[LDT\_ncvar\_PERMSNOWFRAC]
!   LDT NetCDF variable name for Fraction of permanet snow/ice
! \item[LDT\_ncvar\_SLOPE\_DIR]
!   LDT NetCDF variable name for !Typical slope aspect in the grid  (deg from N clockwise)
! \item[LDT\_ncvar\_SAND]
!   LDT NetCDF variable name for Soil sand fraction (-)
! \item[LDT\_ncvar\_SILT]
!   LDT NetCDF variable name for Soil silt fraction (-)
! \item[LDT\_ncvar\_CLAY]
!   LDT NetCDF variable name for Soil clay fraction (-)
! \item[LDT\_ncvar\_POROSITY]
!   LDT NetCDF variable name for Soil porosity (m3 m-3)
! \item[use\_monthly\_albedo\_map]
! if usemonalb == .true., then the alb value passed to lsmcrocus will be used as the background snow-free albedo term.  
! if usemonalb == .false., then alb will be set to 0.2 
! \item[ts]
!   Crocus81 model time step in second
! \item[count]
!   variable to keep track of the number of timesteps before an output
! \item[rstInterval]
!   restart writing interval
! \item[outInterval]
!   output writing interval
! \item[crocus81]
!  Crocus81 model specific variables
! \item[forc\_count]
!   counter of forcing data
! \item[nsnow]
!   number of snow layer
! \item[nimpur]
!   number of impurtites
! \item[SNOWRES\_opt]
!   SNOWRES\_opt  = ISBA-SNOW3L turbulant exchange option
!   'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!   'RIL' = Limit Richarson number under very stable
! conditions (currently testing)
!   'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus
! \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  10/18/19: Mahdi Navari, Shugong Wang Initial implementation for LIS 7 and Crocus81
!  19 Jan 2021: Mahdi Navari, Edited to initialize forcing variables to zeros
!   21 Jan 2021: Mahdi Navari, Ground temperature removed form the lis.config, 
!                        and this subroutine edited to reflect the edit. For the stand-alone
!                        version, the value of TG was set to 273.15 in the Crocus81_main.F90
!
module Crocus81_lsmMod

! !USES:
  use Crocus81_module
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
  
  PRIVATE
  !-------------------------------------------------------------------------
  ! PUBLIC MEMBER FUNCTIONS
  !-------------------------------------------------------------------------
  public :: Crocus81_ini
  !-------------------------------------------------------------------------
  ! PUBLIC TYPES
  !-------------------------------------------------------------------------
  public :: Crocus81_struc
  !EOP
  type, public :: Crocus81_type_dec
     character(len=LIS_CONST_PATH_LEN) :: rfile
     character*256      :: rformat
     !-------------------------------------------------------------------------
     ! Parameter file names
     !-------------------------------------------------------------------------
     !character*128      :: LDT_ncvar_GLACIER_BOOL
     !character*128      :: LDT_ncvar_TG
     character*128      :: LDT_ncvar_SLOPE
     character*128      :: LDT_ncvar_ALB
     !character*128      :: LDT_ncvar_SOILCOND
     character*128      :: LDT_ncvar_PERMSNOWFRAC
     character*128      :: LDT_ncvar_SLOPE_DIR
     character*128      :: LDT_ncvar_SAND
     character*128      :: LDT_ncvar_SILT
     character*128      :: LDT_ncvar_CLAY
     character*128      :: LDT_ncvar_POROSITY   
     !-------------------------------------------------------------------------
     ! ts, Count, rstInterval, outInterval
     !-------------------------------------------------------------------------
     real               :: ts
     integer            :: count
     real               :: rstInterval
     integer            :: outInterval
     integer            :: forc_count
     integer            :: isba_param_count ! MN added to read parameter file isba
     !-------------------------------------------------------------------------
     ! Initial Model State for cold start
     !-------------------------------------------------------------------------
     REAL, pointer      :: init_SNOWSWE(:)
     REAL, pointer      :: init_SNOWRHO(:)
     REAL, pointer      :: init_SNOWHEAT(:)
     REAL               :: init_SNOWALB
     REAL, pointer      :: init_SNOWGRAN1(:)
     REAL, pointer      :: init_SNOWGRAN2(:)
     REAL, pointer      :: init_SNOWHIST(:)
     REAL, pointer      :: init_SNOWAGE(:)
     REAL, pointer      :: init_SNOWLIQ(:)
     REAL, pointer      :: init_SNOWTEMP(:)
     REAL, pointer      :: init_SNOWDZ(:)
     REAL               :: init_GRNDFLUX
     REAL               :: init_SNDRIFT
     REAL               :: init_RI_n
     REAL               :: init_CDSNOW
     REAL               :: init_USTARSNOW
     REAL               :: init_CHSNOW
     REAL               :: init_SNOWMAK_dz
     !-------------------------------------------------------------------------
     ! Constant Parameter
     !-------------------------------------------------------------------------
     integer            :: nsnow
     integer            :: nimpur
     CHARACTER(len=3)   :: SNOWRES_opt
     LOGICAL            :: OMEB_BOOL
     LOGICAL            :: GLACIER_BOOL ! it is a spatial param. It will be computed 
     ! for each grid cell in the driver using PERMSNOWFRAC    
     CHARACTER(len=3)   :: HIMPLICIT_WIND_opt
     REAL               :: PTSTEP
     REAL               :: TG
     REAL               :: UREF
     !REAL               :: SLOPE  
     REAL               :: ZREF
     REAL               :: Z0NAT
     REAL               :: Z0EFF
     REAL               :: Z0HNAT
     !REAL               :: ALB
     logical            :: use_monthly_albedo_map
     REAL               :: D_G
     !REAL               :: PERMSNOWFRAC          
     CHARACTER(len=4)   :: SNOWDRIFT_opt
     LOGICAL            :: SNOWDRIFT_SUBLIM_BOOL
     LOGICAL            :: SNOW_ABS_ZENITH_BOOL
     CHARACTER(len=3)   :: SNOWMETAMO_opt
     CHARACTER(len=3)   :: SNOWRAD_opt
     LOGICAL            :: ATMORAD_BOOL
     REAL, pointer      :: IMPWET(:)
     REAL, pointer      :: IMPDRY(:)
     CHARACTER(len=3)   :: SNOWFALL_opt
     CHARACTER(len=3)   :: SNOWCOND_opt
     CHARACTER(len=3)   :: SNOWHOLD_opt
     CHARACTER(len=3)   :: SNOWCOMP_opt
     CHARACTER(len=3)   :: SNOWZREF_opt
     LOGICAL            :: SNOWCOMPACT_BOOL
     LOGICAL            :: SNOWMAK_BOOL
     LOGICAL            :: SNOWTILLER_BOOL
     LOGICAL            :: SELF_PROD_BOOL
     LOGICAL            :: SNOWMAK_PROP_BOOL
     LOGICAL            :: PRODSNOWMAK_BOOL
     !REAL               :: SLOPE_DIR 
     LOGICAL            :: Partition_total_precip_BOOL     
     type(Crocus81dec), pointer :: crocus81(:)
  end type Crocus81_type_dec
  
  type(Crocus81_type_dec), pointer :: CROCUS81_struc(:)
  
contains 
  
!BOP
!
! !ROUTINE: Crocus81_ini
! \label{Crocus81_ini}
!
! !INTERFACE:
  subroutine Crocus81_ini(kk)
! !USES:
    use ESMF
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod
    use LIS_surfaceModelDataMod
    use LIS_lsmMod

! !DESCRIPTION:
!
!  This routine creates the datatypes and allocates memory for Crocus81-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for Crocus81 from the configuration file.
!
!  The routines invoked are:
!  \begin{description}
!   \item[Crocus81\_readcrd](\ref{Crocus81_readcrd})\\
!    reads the runtime options for Crocus81 model
!  \end{description}
!EOP
    implicit none        
    integer, intent(in) :: kk
    
    integer  :: n, t     
    character*3             :: fnest  !MN  Bug in the toolkit 
    integer  :: status   
    
    type(ESMF_ArraySpec) :: arrspec1
    type(ESMF_Field)     :: gtField, XWGIField, XWGField, sweField, snwdField
    
    
    ! allocate memory for nest 
    allocate(CROCUS81_struc(LIS_rc%nnest))
    
    ! read configuation information from lis.config file
    call Crocus81_readcrd()

        do n=1, LIS_rc%nnest
            ! allocate memory for all tiles in current nest 
            allocate(CROCUS81_struc(n)%crocus81(LIS_rc%npatch(n, LIS_rc%lsm_index)))
            !------------------------------------------------------------------------
            ! allocate memory for vector variables passed to model interfaces        
            ! TODO: check the following allocation statements carefully!
            !------------------------------------------------------------------------
            ! allocate memory for state variables
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWSWE(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWRHO(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWHEAT(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWGRAN1(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWGRAN2(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWHIST(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWAGE(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWLIQ(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWTEMP(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%SNOWDZ(CROCUS81_struc(n)%nsnow))
                allocate(CROCUS81_struc(n)%crocus81(t)%ALB(12))
                
               ! Crocus81_struc(n)%crocus81(t)%tg = Crocus81_struc(n)%tg
            enddo
            ! initialize forcing variables to zeros
            do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
                CROCUS81_struc(n)%crocus81(t)%ta = 0.0
                CROCUS81_struc(n)%crocus81(t)%qa = 0.0
                CROCUS81_struc(n)%crocus81(t)%wind_e = 0.0
                CROCUS81_struc(n)%crocus81(t)%wind_n = 0.0
                CROCUS81_struc(n)%crocus81(t)%rrsnow = 0.0
                CROCUS81_struc(n)%crocus81(t)%srsnow = 0.0
                CROCUS81_struc(n)%crocus81(t)%lw_rad = 0.0
                CROCUS81_struc(n)%crocus81(t)%sw_rad = 0.0
                CROCUS81_struc(n)%crocus81(t)%pps = 0.0
            enddo ! end of tile (t) loop

            !------------------------------------------------------------------------
            ! Model timestep Alarm
            !------------------------------------------------------------------------
            CROCUS81_struc(n)%forc_count = 0
            CROCUS81_struc(n)%isba_param_count = 1 ! MN  isba

            call LIS_update_timestep(LIS_rc, n, CROCUS81_struc(n)%ts)

            write(fnest,'(i3.3)') n !MN  Bug in the toolkit 
            call LIS_registerAlarm("CROCUS81 model alarm "//trim(fnest),&
                                   CROCUS81_struc(n)%ts, &
                                   CROCUS81_struc(n)%ts)

            call LIS_registerAlarm("CROCUS81 restart alarm "//trim(fnest),&
                                   CROCUS81_struc(n)%ts,&
                                   CROCUS81_struc(n)%rstInterval)
            !------------------------------------------------------------------------
            ! TODO: setup number of soil moisture/temperature layers and depth here  
            !------------------------------------------------------------------------
            ! TODO: set number of soil moisture layers in surface model
            LIS_sfmodel_struc(n)%nsm_layers = 1     
            ! TODO: set number of soil temperature layers in surface model
            LIS_sfmodel_struc(n)%nst_layers = 1 
            if(.not.allocated(LIS_sfmodel_struc(n)%lyrthk)) then 
               allocate(LIS_sfmodel_struc(n)%lyrthk(1))
            endif

            LIS_sfmodel_struc(n)%lyrthk(1) = CROCUS81_struc(n)%D_G ! MN added 

            LIS_sfmodel_struc(n)%ts = CROCUS81_struc(n)%ts

            !Create fields for LSM2SUBLSM exchanges
            call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
                 rc=status)
            call LIS_verify(status, &
                 "ESMF_ArraySpecSet failed in Crocus81_in")

            gtField = ESMF_FieldCreate(&
                 grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
                 arrayspec=arrspec1, &
                 name="Ground temperature",&
                 rc=status)
            call LIS_verify(status,&
                 'ESMF_FieldCreate failed in Crocus81_ini')
            
            call ESMF_StateAdd(LIS_LSM2SUBLSM_State(n,kk),&
                 (/gtField/),rc=status)
            call LIS_verify(status,&
                 'ESMF_StateAdd failed in Crocus81_ini')

            XWGIField = ESMF_FieldCreate(&
                 grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
                 arrayspec=arrspec1, &
                 name="soil volumetric frozen water content",&
                 rc=status)
            call LIS_verify(status,&
                 'ESMF_FieldCreate failed in Crocus81_ini')

            call ESMF_StateAdd(LIS_LSM2SUBLSM_State(n,kk),&
                 (/XWGIField/),rc=status)
            call LIS_verify(status,&
                 'ESMF_StateAdd failed in Crocus81_ini')

            XWGField = ESMF_FieldCreate(&
                 grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
                 arrayspec=arrspec1, &
                 name="soil volumetric liquid water content",&
                 rc=status)
            call LIS_verify(status,&
                 'ESMF_FieldCreate failed in Crocus81_ini')

            call ESMF_StateAdd(LIS_LSM2SUBLSM_State(n,kk),&
                 (/XWGField/),rc=status)
            call LIS_verify(status,&
                 'ESMF_StateAdd failed in Crocus81_ini')

            sweField = ESMF_FieldCreate(&
                 grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
                 arrayspec=arrspec1, &
                 name="Total SWE",&
                 rc=status)
            call LIS_verify(status,&
                 'ESMF_FieldCreate failed in Crocus81_ini')

            snwdField = ESMF_FieldCreate(&
                 grid=LIS_vecPatch(n,LIS_rc%lsm_index),&
                 arrayspec=arrspec1, &
                 name="Total snowdepth",&
                 rc=status)
            call LIS_verify(status,&
                 'ESMF_FieldCreate failed in Crocus81_ini')
            
            call ESMF_StateAdd(LIS_SUBLSM2LSM_State(n,kk),&
                 (/sweField/),rc=status)
            call ESMF_StateAdd(LIS_SUBLSM2LSM_State(n,kk),&
                 (/snwdField/),rc=status)
            call LIS_verify(status,&
                 'ESMF_StateAdd failed in Crocus81_ini')

        enddo
    end subroutine Crocus81_ini
end module Crocus81_lsmMod
