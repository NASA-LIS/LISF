!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module clm2_lsmMod
!BOP
!
! !MODULE: clm2_lsmMod
!
! !DESCRIPTION:
!
! This module provides the definition of derived data type used to 
! control the operation of CLM. It also provides the entry method
! for the initialization of CLM-specific variables. The derived
! data type {\tt clm\_struc} includes the variables that specify
! the runtime options and other control variables as described below:
!
! \begin{description}
!  \item[clm\_rfile]
!    name of the CLM restart file
!  \item[clm\_vfile]
!    name of the CLM vegetation parameter lookup table
!  \item[clm\_chtfile]
!    name of the CLM canopy heights lookup table
!  \item[count]
!    variable to keep track of the number of timesteps before an output
!  \item[clmopen]
!    variable to keep track of opened files
!  \item[numout]
!    number of output times 
!  \item[clm\_ism]
!    initial soil moisture for a cold start run
!  \item[clm\_it]
!    initial soil temperature for a cold start run
!  \item[clm\_iscv]
!    initial snow mass for a cold start run
!  \item[outInterval]
!    output writing interval
!  \item[rstInterval]
!   restart writing interval
!  \item[clm]
!   CLM specific variables
! \end{description} 
!
! !USES:
  use clm2type
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: clm2_lsm_ini
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: clm2_struc
!EOP

  type, public :: clm_type_dec
     character(len=LIS_CONST_PATH_LEN) :: clm_rfile 
     character(len=LIS_CONST_PATH_LEN) :: clm_vfile 
     character(len=LIS_CONST_PATH_LEN) :: clm_chtfile 
     integer      :: count
     integer      :: clmopen
     integer      :: numout
     real         :: clm_ism          !clm intial soil moisture
     real         :: clm_it           !clm initial soil temperature
     real         :: clm_iscv         !clm initial snow mass
     real         :: rstinterval
     integer      :: forc_count
     real         :: ts
     type (clm1d), allocatable :: clm(:)     
  end type clm_type_dec

  type(clm_type_dec), allocatable :: clm2_struc(:)
  save
  
contains

!BOP
! 
! !ROUTINE: clm2_lsm_ini
! \label{clm2_lsm_ini}
! 
! !INTERFACE:
  subroutine clm2_lsm_ini()
! !USES:
    use clm2_varpar
    use pft_varcon     
    use clm2_shr_orb_mod,  only : clm2_shr_orb_params
    use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
         LIS_update_timestep, LIS_registerAlarm
    use LIS_logMod, only : LIS_logunit, LIS_verify
    use LIS_coreMod, only : LIS_rc
    use LIS_precisionMod

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for CLM-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for CLM from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readclm2crd](\ref{readclm2crd})
!    reads the runtime options for CLM
!   \item[clm\_varder\_init](\ref{clm2_lsm_init})
!    set initial values to CLM variables
!   \item[iniTimeConst](\ref{iniTimeConst})
!    initialize time invariant parameters
!  \end{description}
!EOP
    implicit none

    integer :: n,i
    integer                 :: numpft_adj
    integer                 :: status
    character*3   :: fnest
  
    
    allocate(clm2_struc(LIS_rc%nnest))
    call readclm2crd()

    do n=1,LIS_rc%nnest
       allocate (clm2_struc(n)%clm(LIS_rc%npatch(n,LIS_rc%lsm_index)))
       write(LIS_logunit,*)'MSG: clm2_lsm_ini -- allocating clm'

       maxpatch_pft = LIS_rc%nvegtypes
       npatch_urban = maxpatch_pft + 1
       npatch_lake = npatch_urban + 1
       npatch_wet = npatch_lake + 1
       npatch_gla = npatch_wet + 1
       maxpatch = npatch_gla

       call clm2_lsm_init(n)

!------------------------------------------------------------------------
! Model timestep Alarm
!------------------------------------------------------------------------
       call LIS_update_timestep(LIS_rc, n, clm2_struc(n)%ts)

       write(fnest,'(i3.3)') n    

       call LIS_registerAlarm("CLM2 model alarm "//trim(fnest),&
            clm2_struc(n)%ts,&
            clm2_struc(n)%ts)
       
       call LIS_registerAlarm("CLM2 restart alarm "//trim(fnest),&
            clm2_struc(n)%ts,&
            clm2_struc(n)%rstInterval)

       clm2_struc(n)%numout = 0 
       clm2_struc(n)%count = 0 
       clm2_struc(n)%clmopen = 0 

       call clm2_shr_orb_params

       write (LIS_logunit,*) 'Attempting to initialize the land model .....'
       write (LIS_logunit,*)

! ----------------------------------------------------------------------
! Read list of PFTs and their corresponding parameter values
! ----------------------------------------------------------------------
       write (LIS_logunit,*) 'Attempting to read PFT physiological data ..'

       open(111, file=clm2_struc(n)%clm_vfile,form='formatted',&
            status='old')

!       if(LIS_rc%inc_water_pts) then
!          numpft_adj = numpft-2
!       else
          numpft_adj = numpft-3
!       endif
       do i = 1, numpft_adj

#if(defined PFT_MODE) 
          read (111,*)  pftname(i),              &
               z0mr(i)   , displar(i), dleaf(i)  , c3psn(i)  , &
               vcmx25(i) , mp(i)     , qe25(i)   , rhol(i,1) , &
               rhol(i,2) , rhos(i,1) , rhos(i,2) , taul(i,1) , &
               taul(i,2) , taus(i,1) , taus(i,2) , xl(i)     , &
               roota_par(i), rootb_par(i)

#else                    
          read (111,60)  pftname(i),              &
               z0mr(i)   , displar(i), dleaf(i)  , c3psn(i)  , &
               vcmx25(i) , mp(i)     , qe25(i)   , rhol(i,1) , &
               rhol(i,2) , rhos(i,1) , rhos(i,2) , taul(i,1) , &
               taul(i,2) , taus(i,1) , taus(i,2) , xl(i)     , &
               roota_par(i), rootb_par(i)
#endif
       end do
       close(111)
60     format(A27,13x,f5.3,1x,f4.2,1x,f4.2,1x,f2.0,1x,f3.0,1x,f2.0,1x,&
            7(f4.2,1x),2(f5.3,1x),f5.2,1x,f4.1,1x,f3.1)

! ----------------------------------------------------------------------
! Initialize time invariant variables as subgrid vectors of length [numpatch] 
! ----------------------------------------------------------------------
       write (LIS_logunit,*) ('Attempting to initialize time invariant variables')
       
       call iniTimeConst (n)
       write (LIS_logunit,*) ('Successfully initialized time invariant variables')
       write (LIS_logunit,*)
    enddo
  end subroutine clm2_lsm_ini

!BOP
! 
! !ROUTINE: clm2_lsm_init
! \label{clm2_lsm_init}
! 
! !INTERFACE:
  subroutine clm2_lsm_init(n)
! !USES:
    use LIS_coreMod, only : LIS_rc
    use LIS_precisionMod
    use clm2_varpar,   only : nlevsoi, nlevsno, numrad
! !DESCRIPTION: 
! 
! Initializes clm variables
! 
!  The arguments are: 
! \begin{description}
!  \item[n] 
!   index of the nest
! \end{description}
!EOP
    integer :: k, n    
    do k = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
       clm2_struc(n)%clm(k)%kpatch  = bigint        
       clm2_struc(n)%clm(k)%itypveg = bigint        
       clm2_struc(n)%clm(k)%itypwat = bigint        
       clm2_struc(n)%clm(k)%isoicol = bigint        

! level values

       clm2_struc(n)%clm(k)%dz(-nlevsno+1:0) = inf       !snow layer thickness (m)
       clm2_struc(n)%clm(k)%z (-nlevsno+1:0) = inf       !snow layer depth (m)
       clm2_struc(n)%clm(k)%zi(-nlevsno+0:0) = inf       !snow layer interfaces (m)  
       clm2_struc(n)%clm(k)%dz(1:nlevsoi)    = inf       !soil layer thickness (m)
       clm2_struc(n)%clm(k)%z (1:nlevsoi)    = inf       !soil layer depth (m)
       clm2_struc(n)%clm(k)%zi(1:nlevsoi)    = inf       !soil layer interfaces (m)  

! soil physical properties

       clm2_struc(n)%clm(k)%bsw   (1:nlevsoi) = inf      !Clapp and Hornberger "b"
       clm2_struc(n)%clm(k)%watsat(1:nlevsoi) = inf      !volumetric soil water at saturation (porosity)
       clm2_struc(n)%clm(k)%hksat (1:nlevsoi) = inf      !hydraulic conductivity at saturation (mm H2O /s)
       clm2_struc(n)%clm(k)%sucsat(1:nlevsoi) = inf      !minimum soil suction (mm)
       clm2_struc(n)%clm(k)%csol  (1:nlevsoi) = inf      !heat capacity, soil solids (J/m**3/Kelvin)
       clm2_struc(n)%clm(k)%tkmg  (1:nlevsoi) = inf      !thermal conductivity, soil minerals  [W/m-K]  (new)
       clm2_struc(n)%clm(k)%tkdry (1:nlevsoi) = inf      !thermal conductivity, dry soil       (W/m/Kelvin)
       clm2_struc(n)%clm(k)%tksatu(1:nlevsoi) = inf      !thermal conductivity, saturated soil [W/m-K]  (new)
       clm2_struc(n)%clm(k)%rootfr(1:nlevsoi) = inf      !fraction of roots in each soil layer
       clm2_struc(n)%clm(k)%rootr (1:nlevsoi) = inf      !effective fraction of roots in each soil layer

! leaf constants

!       clm2_struc(n)%clm(k)%dewmx = inf            !Maximum allowed dew [mm]

! hydraulic constants of soil      

!       clm2_struc(n)%clm(k)%wtfact = inf           !Fraction of model area with high water table
!       clm2_struc(n)%clm(k)%trsmx0 = inf           !Max transpiration for moist soil+100% veg. [mm s-1]

! roughness lengths     

!       clm2_struc(n)%clm(k)%zlnd   = inf           !Roughness length for soil [m] (new)             
!       clm2_struc(n)%clm(k)%zsno   = inf           !Roughness length for snow [m] (new)             
!       clm2_struc(n)%clm(k)%csoilc = inf           !Drag coefficient for soil under canopy [-] (new)

! numerical finite-difference

!       clm2_struc(n)%clm(k)%cnfac   = inf          !Crank Nicholson factor (between 0 and 1) (new)
!       clm2_struc(n)%clm(k)%capr    = inf          !Tuning factor to turn first layer T into surface T (new)  
!       clm2_struc(n)%clm(k)%ssi     = inf          !Irreducible water saturation of snow (new)
!       clm2_struc(n)%clm(k)%wimp    = inf          !Water impremeable if porosity less than wimp (new)
!       clm2_struc(n)%clm(k)%pondmx  = inf          !Ponding depth (mm) (new)
!       clm2_struc(n)%clm(k)%smpmax  = inf          !wilting point potential in mm (new)
!       clm2_struc(n)%clm(k)%smpmin  = inf          !restriction for min of soil potential (mm) (new)

! water and energy balance check

       clm2_struc(n)%clm(k)%begwb      = inf       !water mass begining of the time step
       clm2_struc(n)%clm(k)%endwb      = inf       !water mass end of the time step
!       clm2_struc(n)%clm(k)%errh2o     = inf       !water conservation error (mm H2O)
!       clm2_struc(n)%clm(k)%errsoi     = inf       !soil/lake energy conservation error (W/m**2)
!       clm2_struc(n)%clm(k)%errseb     = inf       !surface energy conservation error (W/m**2)
!       clm2_struc(n)%clm(k)%errsol     = inf       !solar radiation conservation error (W/m**2)
!       clm2_struc(n)%clm(k)%errlon     = inf       !longwave radiation conservation error (W/m**2)
!       clm2_struc(n)%clm(k)%acc_errseb = 0.        !accumulation of surface energy balance error
!       clm2_struc(n)%clm(k)%acc_errh2o = 0.        !accumulation of water balance error

!*************************************************************************
! subgrid patch version of atm model input
!*************************************************************************
       clm2_struc(n)%clm(k)%forc_cosz  = inf         !cosine of solar zenith angle
       clm2_struc(n)%clm(k)%forc_t     = inf         !atmospheric temperature (Kelvin)
       clm2_struc(n)%clm(k)%forc_u     = inf         !atmospheric wind speed in east direction (m s-1)
       clm2_struc(n)%clm(k)%forc_v     = inf         !atmospheric wind speed in north direction (m s-1)
       clm2_struc(n)%clm(k)%forc_q     = inf         !atmospheric specific humidity (kg kg-1)
       clm2_struc(n)%clm(k)%forc_hgt   = inf         !atmospheric reference height (m) 
       clm2_struc(n)%clm(k)%forc_hgt_u = inf         !observational height of wind [m] (new)
       clm2_struc(n)%clm(k)%forc_hgt_t = inf         !observational height of temperature [m] (new)
       clm2_struc(n)%clm(k)%forc_hgt_q = inf         !observational height of humidity [m] (new)
       clm2_struc(n)%clm(k)%forc_pbot  = inf         !atmospheric pressure (Pa)
       clm2_struc(n)%clm(k)%forc_th    = inf         !atmospheric potential temperature (Kelvin)
!       clm2_struc(n)%clm(k)%forc_vp    = inf         !atmospheric vapor pressure (Pa)
       clm2_struc(n)%clm(k)%forc_rho   = inf         !density (kg/m**3)
!       clm2_struc(n)%clm(k)%forc_co2   = inf         !atmospheric CO2 concentration (Pa)
!       clm2_struc(n)%clm(k)%forc_o2    = inf         !atmospheric O2 concentration (Pa)
       clm2_struc(n)%clm(k)%forc_lwrad = inf         !downward infrared (longwave) radiation (W/m**2)
!       clm2_struc(n)%clm(k)%forc_psrf  = inf         !surface pressure (Pa)
       clm2_struc(n)%clm(k)%forc_solad(1:numrad) = inf !direct beam radiation (vis=forc_sols , nir=forc_soll )
       clm2_struc(n)%clm(k)%forc_solai(1:numrad) = inf !diffuse radiation     (vis=forc_solsd, nir=forc_solld)
!clm+
       clm2_struc(n)%clm(k)%forc_rain  = inf         !rain rate [mm s-1]
       clm2_struc(n)%clm(k)%forc_snow  = inf         !snow rate [mm s-1]
!clm-
!       clm2_struc(n)%clm(k)%forc_swc1    = inf       !Layer 1 (0-10 cm) soil water content (% capacity for GEOS 
!       clm2_struc(n)%clm(k)%forc_sdepth  = inf       !Model liquid equivalent snow depth (kg m-2) (same for GEOS)

!*************************************************************************
! biogeophys
!*************************************************************************

! Surface solar radiation 

       clm2_struc(n)%clm(k)%rssun  = inf        !sunlit stomatal resistance (s/m)
       clm2_struc(n)%clm(k)%rssha  = inf        !shaded stomatal resistance (s/m)
       clm2_struc(n)%clm(k)%psnsun = inf        !sunlit leaf photosynthesis (umol CO2 /m**2/ s) 
       clm2_struc(n)%clm(k)%psnsha = inf        !shaded leaf photosynthesis (umol CO2 /m**2/ s)
       clm2_struc(n)%clm(k)%laisun = inf        !sunlit leaf area
       clm2_struc(n)%clm(k)%laisha = inf        !shaded leaf area
!       clm2_struc(n)%clm(k)%ndvi   = inf        !Normalized Difference Vegetation Index
       clm2_struc(n)%clm(k)%sabg   = inf        !solar radiation absorbed by ground (W/m**2)
       clm2_struc(n)%clm(k)%sabv   = inf        !solar radiation absorbed by vegetation (W/m**2)
       clm2_struc(n)%clm(k)%fsa    = inf        !solar radiation absorbed (total) (W/m**2)
!       clm2_struc(n)%clm(k)%fsr    = inf        !solar radiation reflected (W/m**2)

! Surface energy fluxes

       clm2_struc(n)%clm(k)%taux           = inf !wind stress: e-w (kg/m s-1**2)
       clm2_struc(n)%clm(k)%tauy           = inf !wind stress: n-s (kg/m s-1**2)
       clm2_struc(n)%clm(k)%eflx_lwrad_out = inf !emitted infrared (longwave) radiation (W/m**2) 
       clm2_struc(n)%clm(k)%eflx_lwrad_net = inf !net infrared (longwave) rad (W/m**2) [+  = to atm]
       clm2_struc(n)%clm(k)%eflx_sh_tot    = inf !total sensible heat flux (W/m**2) [+ to atm]
       clm2_struc(n)%clm(k)%eflx_sh_veg    = inf !sensible heat flux from leaves (W/m**2) [+ to atm]
       clm2_struc(n)%clm(k)%eflx_sh_grnd   = inf !sensible heat flux from ground (W/m**2) [+ to atm]
       clm2_struc(n)%clm(k)%eflx_lh_tot    = inf !total latent heat flux (W/m8*2)  [+ to atm] 
!       clm2_struc(n)%clm(k)%eflx_lh_vege   = inf !veg evaporation heat flux (W/m**2) [+ to atm]
!       clm2_struc(n)%clm(k)%eflx_lh_vegt   = inf !veg transpiration heat flux (W/m**2) [+ to atm]
!       clm2_struc(n)%clm(k)%eflx_lh_grnd   = inf !ground evaporation heat flux (W/m**2) [+ to atm]   
       clm2_struc(n)%clm(k)%eflx_soil_grnd = inf !soil heat flux (W/m**2) [+  = into soil]
!       clm2_struc(n)%clm(k)%eflx_snomelt   = inf !snow melt heat flux (W/m**2)

! velocities

!       clm2_struc(n)%clm(k)%u10 = inf            !10-m wind (m s-1)
!       clm2_struc(n)%clm(k)%fv  = inf            !friction velocity (m s-1)

! Temperatures

       clm2_struc(n)%clm(k)%t_veg   = inf                     !vegetation temperature (Kelvin)
       clm2_struc(n)%clm(k)%t_grnd  = inf                     !ground temperature (Kelvin)
       clm2_struc(n)%clm(k)%t_rad   = inf                     !radiative temperature (Kelvin)
       clm2_struc(n)%clm(k)%t_ref2m = inf                     !2 m height surface air temperature (Kelvin)
       clm2_struc(n)%clm(k)%t_soisno(-nlevsno+1:0) = 0      !snow temperature (Kelvin)
       clm2_struc(n)%clm(k)%t_soisno(1:nlevsoi)    = 0      !soil temperature (Kelvin)
!       clm2_struc(n)%clm(k)%t_lake(1:nlevlak)      = inf      !lak temperature (Kelvin)
!       clm2_struc(n)%clm(k)%dt_veg  = spval                   !change in t_veg, last iteration (Kelvin)
!       clm2_struc(n)%clm(k)%dt_grnd = spval                   !change in t_grnd, last iteration (Kelvin)

! Soil properties   

!       clm2_struc(n)%clm(k)%btran = inf                       !transpiration wetness factor (0 to 1) 

!*************************************************************************
! biogeochem
!*************************************************************************

!       clm2_struc(n)%clm(k)%fpsn  = inf         !photosynthesis (umol CO2 /m**2 /s)
!       clm2_struc(n)%clm(k)%frm   = inf         !total maintenance respiration (umol CO2 /m**2/s)
!       clm2_struc(n)%clm(k)%frmf  = inf         !leaf maintenance respiration  (umol CO2 /m**2 /s)
!       clm2_struc(n)%clm(k)%frms  = inf         !stem maintenance respiration  (umol CO2 /m**2 /s)
!       clm2_struc(n)%clm(k)%frmr  = inf         !root maintenance respiration  (umol CO2 /m**2 /s)
!       clm2_struc(n)%clm(k)%frg   = inf         !growth respiration (umol CO2 /m**2 /s)
!       clm2_struc(n)%clm(k)%fmicr = inf         !microbial respiration (umol CO2 /m**2 /s)
!       clm2_struc(n)%clm(k)%fco2  = inf         !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
!       clm2_struc(n)%clm(k)%dmi   = inf         !total dry matter production (ug /m**2 /s)

! dust model

!       clm2_struc(n)%clm(k)%flx_mss_vrt_dst(1:ndst) = inf !surface dust emission (kg/m**2/s) [ + = to atm]
!       clm2_struc(n)%clm(k)%vwc_thr                 = inf !threshold soil moisture based on clay content
!       clm2_struc(n)%clm(k)%mss_frc_cly_vld         = inf ![frc] Mass fraction clay limited to 0.20
!       clm2_struc(n)%clm(k)%mbl_bsn_fct             = inf

! voc model

!       clm2_struc(n)%clm(k)%vocflx(1:nvoc) = inf          !VOC flux [ug C m-2 h-1]

!*************************************************************************
! hydrology
!*************************************************************************

       clm2_struc(n)%clm(k)%qflx_infl       = inf                !Infiltration (mm H2O /s) 
       clm2_struc(n)%clm(k)%qflx_surf       = inf                !surface runoff (mm H2O /s) 
       clm2_struc(n)%clm(k)%qflx_drain      = inf                !sub-surface runoff (mm H2O /s) 
       clm2_struc(n)%clm(k)%qflx_top_soil   = inf                !net water input into soil from top (mm s-1)
       clm2_struc(n)%clm(k)%qflx_evap_soi   = inf                !soil evaporation (mm H2O/s) (+ = to atm)
       clm2_struc(n)%clm(k)%qflx_evap_veg   = inf                !vegetation evaporation (mm H2O/s) (+ = to atm)
       clm2_struc(n)%clm(k)%qflx_tran_veg   = inf                !vegetation transpiration (mm H2O/s) (+ = to atm)
       clm2_struc(n)%clm(k)%qflx_snomelt    = inf                !snow melt (mm H2O /s)
       clm2_struc(n)%clm(k)%qflx_evap_tot   = inf                !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
!       clm2_struc(n)%clm(k)%qflx_prec_intr  = inf                !interception of precipitation [mm s-1]
!       clm2_struc(n)%clm(k)%qflx_prec_grnd  = inf                !water onto ground including canopy runoff [kg/(m2 s)]
       clm2_struc(n)%clm(k)%qflx_rain_grnd  = inf                !rain on ground after interception (mm H2O/s) [+]
!       clm2_struc(n)%clm(k)%qflx_snow_grnd  = inf                !snow on ground after interception (mm H2O/s) [+]
       clm2_struc(n)%clm(k)%qflx_evap_grnd  = inf                !ground surface evaporation rate (mm H2O/s) [+]
       clm2_struc(n)%clm(k)%qflx_dew_grnd   = inf                !ground surface dew formation (mm H2O /s) [+]
       clm2_struc(n)%clm(k)%qflx_sub_snow   = inf                !sublimation rate from snow pack (mm H2O /s) [+]
       clm2_struc(n)%clm(k)%qflx_dew_snow   = inf                !surface dew added to snow pack (mm H2O /s) [+]
       clm2_struc(n)%clm(k)%qflx_snowcap    = inf                !excess precipitation due to snow capping (mm H2O /s) [+]
       clm2_struc(n)%clm(k)%qflx_qrgwl      = 0                  !qflx_surf at glaciers, wetlands, lakes
       clm2_struc(n)%clm(k)%h2osno          = inf                !snow water (mm H2O / m**2)
       clm2_struc(n)%clm(k)%h2ocan          = inf                !canopy water (mm H2O / m**2)
       clm2_struc(n)%clm(k)%h2osoi_liq(-nlevsno+1:0) = inf       !snow liquid water (kg m-2) (new)
       clm2_struc(n)%clm(k)%h2osoi_ice(-nlevsno+1:0) = inf       !snow ice lens (kg m-2) (new)
       clm2_struc(n)%clm(k)%h2osoi_liq(1:nlevsoi) = inf          !soil liquid water (kg m-2) (new)
       clm2_struc(n)%clm(k)%h2osoi_ice(1:nlevsoi) = inf          !soil ice lens (kg m-2) (new)
       clm2_struc(n)%clm(k)%h2osoi_vol(1:nlevsoi) = inf          !volumetric soil water (0<=h2osoi_vol<=watsat) [m^3 m-3]
       clm2_struc(n)%clm(k)%snowdp          = inf                !snow height (m) 
       clm2_struc(n)%clm(k)%snowage         = inf                !non dimensional snow age [-] (new)
!       clm2_struc(n)%clm(k)%t_snow          = inf                !average snow temperature
!       clm2_struc(n)%clm(k)%snowice         = inf                !average snow ice lens
!       clm2_struc(n)%clm(k)%snowliq         = inf                !average snow liquid water
       clm2_struc(n)%clm(k)%h2osno_old      = inf                !snow mass for previous time step (kg m-2) (new)
       clm2_struc(n)%clm(k)%frac_veg_nosno  = bigint             !fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
       clm2_struc(n)%clm(k)%frac_veg_nosno_alb = bigint          !fraction of vegetation not covered by snow (0 OR 1 now) [-] (new)
       clm2_struc(n)%clm(k)%frac_sno        = inf                !fraction of ground covered by snow (0 to 1) 
       clm2_struc(n)%clm(k)%frac_iceold(-nlevsno+1:nlevsoi)=inf  !fraction of ice relative to the total water (new)
!       clm2_struc(n)%clm(k)%rsw             = inf                !soil water content for root zone
       clm2_struc(n)%clm(k)%eff_porosity    = inf                !effective porosity
!       clm2_struc(n)%clm(k)%sfact           = inf                !term for implicit correction to evaporation
!       clm2_struc(n)%clm(k)%sfactmax        = inf                !maximim of "sfact"

       clm2_struc(n)%clm(k)%imelt(-nlevsno+1:nlevsoi) = bigint   !flag for melting (=1), freezing (=2), Not=0 (new)        

!*************************************************************************
! surfacealbedo (for next time step)
!*************************************************************************

       clm2_struc(n)%clm(k)%parsun          = inf !average absorbed PAR for sunlit leaves (W/m**2)
!       clm2_struc(n)%clm(k)%parsha          = inf !average absorbed PAR for shaded leaves (W/m**2)
!       clm2_struc(n)%clm(k)%albd(1:numrad)  = inf !surface albedo (direct)                     
!       clm2_struc(n)%clm(k)%albi(1:numrad)  = inf !surface albedo (diffuse)                    
       clm2_struc(n)%clm(k)%albgrd(1:numrad)= inf !ground albedo (direct)                      
       clm2_struc(n)%clm(k)%albgri(1:numrad)= inf !ground albedo (diffuse)                     
       clm2_struc(n)%clm(k)%fabd(1:numrad)  = inf !flux absorbed by veg per unit direct flux   
       clm2_struc(n)%clm(k)%fabi(1:numrad)  = inf !flux absorbed by veg per unit diffuse flux  
       clm2_struc(n)%clm(k)%ftdd(1:numrad)  = inf !down direct flux below veg per unit dir flx 
       clm2_struc(n)%clm(k)%ftid(1:numrad)  = inf !down diffuse flux below veg per unit dir flx
       clm2_struc(n)%clm(k)%ftii(1:numrad)  = inf !down diffuse flux below veg per unit dif flx
       clm2_struc(n)%clm(k)%fsun            = inf !sunlit fraction of canopy
       clm2_struc(n)%clm(k)%surfalb         = inf !instantaneous all-wave surface albedo
!       clm2_struc(n)%clm(k)%surfalb_old     = inf !previous instantaneous all-wave surface albedo
!       clm2_struc(n)%clm(k)%snoalb_old      = inf !previous instantaneous all-wave snow albedo

!*************************************************************************
! ecosysdynamics
!*************************************************************************

       clm2_struc(n)%clm(k)%displa        = inf !displacement height [m]
       clm2_struc(n)%clm(k)%z0m           = inf !roughness length, momentum [m]
       clm2_struc(n)%clm(k)%tlai          = inf !one-sided leaf area index, no burying by snow
       clm2_struc(n)%clm(k)%tsai          = inf !one-sided stem area index, no burying by snow
       clm2_struc(n)%clm(k)%elai          = inf !one-sided leaf area index with burying by snow
       clm2_struc(n)%clm(k)%esai          = inf !one-sided stem area index with burying by snow
       
!       clm2_struc(n)%clm(k)%lai_t1_f      = inf
!       clm2_struc(n)%clm(k)%lai_t2_f      = inf
!       clm2_struc(n)%clm(k)%sai_t1_f      = inf
!       clm2_struc(n)%clm(k)%sai_t2_f      = inf
!       clm2_struc(n)%clm(k)%top_t1_f      = inf
!       clm2_struc(n)%clm(k)%top_t2_f      = inf
!       clm2_struc(n)%clm(k)%bot_t1_f      = inf
!       clm2_struc(n)%clm(k)%bot_t2_f      = inf       
       
       clm2_struc(n)%clm(k)%fwet          = inf !fraction of canopy that is wet (0 to 1)
       clm2_struc(n)%clm(k)%fdry          = inf !fraction of foliage that is green and dry [-] (new)
       clm2_struc(n)%clm(k)%hbot          = inf !canopy bottom height [m]
       clm2_struc(n)%clm(k)%htop          = inf !canopy top height [m]
! next set of variables for use with the DGVM
!       clm2_struc(n)%clm(k)%agdd0 = inf
!       clm2_struc(n)%clm(k)%agdd5 = inf
!       clm2_struc(n)%clm(k)%agddtw = inf
!       clm2_struc(n)%clm(k)%agdd = inf
!       clm2_struc(n)%clm(k)%t10 = inf
!       clm2_struc(n)%clm(k)%t_mo = inf
!       clm2_struc(n)%clm(k)%t_mo_min = inf
!       clm2_struc(n)%clm(k)%fnpsn10 = inf
!       clm2_struc(n)%clm(k)%prec365 = inf
!       clm2_struc(n)%clm(k)%agdd20 = inf
!       clm2_struc(n)%clm(k)%tmomin20 = inf
!       clm2_struc(n)%clm(k)%t10min = inf
!       clm2_struc(n)%clm(k)%tsoi25 = inf
       clm2_struc(n)%clm(k)%annpsn = inf
       clm2_struc(n)%clm(k)%annpsnpot = inf
!       clm2_struc(n)%clm(k)%dphen = inf
!       clm2_struc(n)%clm(k)%leafon = inf
!       clm2_struc(n)%clm(k)%leafof = inf
!       clm2_struc(n)%clm(k)%nind = inf
!       clm2_struc(n)%clm(k)%lm_ind = inf
!       clm2_struc(n)%clm(k)%sm_ind = inf
!       clm2_struc(n)%clm(k)%hm_ind = inf
!       clm2_struc(n)%clm(k)%rm_ind = inf
!       clm2_struc(n)%clm(k)%lai_ind = inf
!       clm2_struc(n)%clm(k)%fpcinc = inf
!       clm2_struc(n)%clm(k)%fpcgrid = inf
!       clm2_struc(n)%clm(k)%crownarea = inf
!       clm2_struc(n)%clm(k)%bm_inc = inf
!       clm2_struc(n)%clm(k)%afmicr = inf
!       clm2_struc(n)%clm(k)%firelength = inf
       clm2_struc(n)%clm(k)%wf = inf
!       clm2_struc(n)%clm(k)%litterag = inf
!       clm2_struc(n)%clm(k)%litterbg = inf
!       clm2_struc(n)%clm(k)%cpool_fast = inf
!       clm2_struc(n)%clm(k)%cpool_slow = inf
!       clm2_struc(n)%clm(k)%k_fast_ave = inf
!       clm2_struc(n)%clm(k)%k_slow_ave = inf
!       clm2_struc(n)%clm(k)%litter_decom_ave = inf

!*************************************************************************
! terms due to splitting the code into Biogeophys1 and Biogeophys2
!*************************************************************************

!       clm2_struc(n)%clm(k)%cgrnd  = inf ! deriv of soil energy flux wrt to soil temp [w/m2/k]
!       clm2_struc(n)%clm(k)%cgrndl = inf ! deriv of soil sensible heat flux wrt soil temp [w/m2/k]
!       clm2_struc(n)%clm(k)%cgrnds = inf ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
!       clm2_struc(n)%clm(k)%tg     = inf ! ground surface temperature [K]
!       clm2_struc(n)%clm(k)%tssbef(-nlevsno:nlevsoi) = inf ! soil/snow temperature before update
!       clm2_struc(n)%clm(k)%qg     = inf ! ground specific humidity [kg kg-1]
!       clm2_struc(n)%clm(k)%dqgdT  = inf ! d(qg)/dT
!       clm2_struc(n)%clm(k)%emg    = inf ! ground emissivity
!       clm2_struc(n)%clm(k)%emv    = inf ! vegetation emissivity
!       clm2_struc(n)%clm(k)%htvp   = inf ! latent heat of vapor of water (or sublimation) [j/kg]
!       clm2_struc(n)%clm(k)%z0mg   = inf ! roughness length over ground, momentum [m]
!       clm2_struc(n)%clm(k)%z0hg   = inf ! roughness length over ground, sensible heat [m]
!       clm2_struc(n)%clm(k)%z0qg   = inf ! roughness length over ground, latent heat [m]
!       clm2_struc(n)%clm(k)%z0mv   = inf ! roughness length over vegetation, momentum [m]
!       clm2_struc(n)%clm(k)%z0hv   = inf ! roughness length over vegetation, sensible heat [m]
!       clm2_struc(n)%clm(k)%z0qv   = inf ! roughness length over vegetation, latent heat [m]
!       clm2_struc(n)%clm(k)%beta   = inf ! coefficient of convective velocity [-]
!       clm2_struc(n)%clm(k)%zii    = inf ! convective boundary height [m]
!       clm2_struc(n)%clm(k)%thm    = inf ! intermediate variable (forc_t+0.0098*forc_hgt_t)
!       clm2_struc(n)%clm(k)%thv    = inf ! virtual potential temperature (kelvin)
!       clm2_struc(n)%clm(k)%ur     = inf ! wind speed at reference height [m s-1]
!       clm2_struc(n)%clm(k)%dlrad  = inf ! downward longwave radiation below the canopy [W m-2]
!       clm2_struc(n)%clm(k)%ulrad  = inf ! upward longwave radiation above the canopy [W m-2]
!       clm2_struc(n)%clm(k)%qmelt  = inf ! snow melt [mm s-1]
       


    end do  ! end of patch loop
  end subroutine clm2_lsm_init

end module clm2_lsmMod
