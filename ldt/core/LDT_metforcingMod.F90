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
module LDT_metforcingMod

!BOP
!
!  !MODULE: LDT_metforcingMod
! 
!  !DESCRIPTION:
!   The code in this file provides interfaces to manage different
!   meteorological forcing analyses
! 
!  \subsubsection{Overview} 
!   This module defines interface plugins for the incorporation
!   of various meteorological forcing analyses. The analyses implemented
!   by extending these interfaces should contain all (and possibly more)
!   of the basic meteorological forcing variables. 
!   Analyses containing one or more of the basic
!   variables, but not the entire set, should be implmented as a
!   supplemental forcing. Further, any analysis that does not cover the 
!   entire globe should be implemented as a supplemental forcing. 
!   The following is a list of the basic meterological forcing variables. 
!   
!   \begin{itemize}
!   \item{\tt T 2m: Temperature interpolated to 2m ($K$)}  
!   \item{\tt Q 2m: Instantaneous specific humidity interpolated to 2m ($kg/kg$)}  
!   \item{\tt SWdown: Downward shortwave flux at the ground ($W/m^2$)}  
!   \item{\tt LWdown: Downward shortwave flux at the ground ($W/m^2$)}  
!   \item{\tt U 10m: Instantaneous zonal wind interpolated to 10m ($m/s$)}
!   \item{\tt V 10m: Instantaneous meridional wind interpolated to 10m ($m/s$)}
!   \item{\tt Psurf: Instantaneous surface pressure ($Pa$)}
!   \item{\tt rainf: Total precipitation ($mm/s$)}
!   \item{\tt rainf\_c: Convective precipitation ($mm/s$)}
!   \item{\tt snowf: Total snowfall ($mm/s$)}
!   \item{\tt PET: Potential ET for FEWSNET ($mm/s$)}
!   \item{\tt RefET: Reference ET for FEWSNET ($mm/s$)}
!   \item{\tt CAPE: Convective Available Potential Energy from NLDAS-2 ($J/kg$)}
!   \end{itemize}
!
!  !REVISION HISTORY: 
!  14 Nov 2002    Sujay Kumar  Initial Specification
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  26 May 2011:   Soni Yatheendradas: Added Potential ET forcings for FEWSNET 
!  14 Mar 2014:   David Mocko: Added CAPE as a forcing variable from NLDAS-2.
!                              Moved CRainf (convective rainfall forcing)
!                                 from LDT_MOC_RAINFCONV to LDT_MOC_CRAINFFORC.
!                              Added units of [kg/m^2] for PET and CRainf.
! 
  use ESMF
  use LDT_FORC_AttributesMod
  use LDT_spatialDownscalingMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_coreMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none 

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_metforcingInit      ! initialize met forcing setup
  public :: LDT_get_met_forcing     ! retrieve data and interpolate 
                                    !  spatially and temporally
  public :: LDT_perturb_forcing     ! perturbs the met forcing variables
  public :: LDT_metforcing_reset    ! resets required forcing variables
  public :: LDT_metforcing_finalize ! cleanup allocated structures
  public :: LDT_output_met_forcing
  public :: LDT_resetForcingVars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public   :: LDT_forc      
  public   :: LDT_FORC_State
  public   :: LDT_FORC_Base_State
  public   :: LDT_FORC_Pert_State
!EOP

  type, public :: forc_dec_type
     real, allocatable :: metdata1(:,:) ! previous forcing values
     real, allocatable :: metdata2(:,:) ! next forcing values
     real, allocatable :: metdata3(:,:) ! uber-next forcing values
     real, allocatable :: modelelev(:)
  end type forc_dec_type

  type(forc_dec_type), allocatable :: LDT_forc(:,:)
  type(ESMF_State),    allocatable :: LDT_FORC_State(:)
  type(ESMF_State),    allocatable :: LDT_FORC_Base_State(:,:)
  type(ESMF_State),    allocatable :: LDT_FORC_Pert_State(:)
contains

!BOP
! !ROUTINE: LDT_metforcingInit
! \label{LDT_metforcingInit}
!
! !INTERFACE: 
! !Private name: call using LDT_metforcingInit
  subroutine LDT_metforcingInit

! !USES:
    use LDT_metforcing_pluginMod, only : LDT_metforcing_plugin
!
! !DESCRIPTION:
!
! This routine sets up the structures to include meteorological 
! forcing analyses for an offline run. The registry that 
! defines the implemented forcing schemes are invoked followed
! by the routines to initialize the specific instance of the forcing scheme. 
! 
! The methods invoked are:  
! \begin{description}
!  \item[LDT\_metforcing\_plugin](\ref{LDT_metforcing_plugin}) \newline
!    sets up function table registries for implemented 
!    met forcing analyses
!  \item[create\_forcing\_structures](\ref{create_forcing_structures}) \newline
!    create and allocate memory for required structures
!  \item[initmetforc](\ref{initmetforc}) \newline
!    invokes the generic method in the registry to define the
!    native domain of the met forcing scheme
!  \item[forcingPerturbSetup](\ref{forcingPerturbSetup}) \newline
!    setup and allocate structures required for forcing 
!    perturbations. 
! \end{description}
!EOP
    integer :: n, m, k
    integer :: rc
    character(10) :: time
    integer :: nensem
! _______________________________________________________

    if( LDT_rc%nmetforc > 0 ) then

      write(LDT_logunit,*)" - - - - - - - - Meteorological Forcing Datasets - - - - - - - -"

       do n=1,LDT_rc%nnest
          if(LDT_rc%metforc_blend_alg.eq."ensemble") then
             if(LDT_rc%nensem(n).lt.LDT_rc%nmetforc) then
                write(LDT_logunit,*) 'The number of ensembles can not be smaller than '
                write(LDT_logunit,*) 'the number of forcings in the ensemble of forcings mode'
                write(LDT_logunit,*) 'Program stopping.... '
                call LDT_endrun()
             endif
             if(mod(LDT_rc%nensem(n),LDT_rc%nmetforc).ne.0) then
                write(LDT_logunit,*) 'The number of ensembles must be a multipleof  '
                write(LDT_logunit,*) 'the number of forcings in the ensemble of forcings mode'
                write(LDT_logunit,*) 'Program stopping.... '
                call LDT_endrun()
             else
              ! Equally divide the ensembles ...
                LDT_rc%nperforc = LDT_rc%nensem(n)/LDT_rc%nmetforc
             endif
          endif
       enddo

       allocate( LDT_forc(LDT_rc%nnest,LDT_rc%nmetforc) )

    !- Register Met Forcing routines:
       call LDT_metforcing_plugin

    !- Continue, running for metforcing field processing:
       call ESMF_ConfigGetAttribute(LDT_config,time,&
            label="Processed metforcing output interval:",rc=rc)
       call LDT_verify(rc,'Processed metforcing output interval: not defined')
       
       call LDT_parseTimeString(time,LDT_rc%metForcOutInterval)
       
       call LDT_registerAlarm("LDT metforcing output alarm",&
            LDT_rc%ts,&
            LDT_rc%metForcOutInterval)
       
       call LDT_update_timestep(LDT_rc, 1, LDT_rc%metForcOutInterval)

       write(LDT_logunit,*) " LDT Processed metforcing output interval: ",&
             LDT_rc%metForcOutInterval
       
    !- Create and allocate memory for the forcing variables:
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%forcvarlistfile,&
            label="Forcing variables list file:",rc=rc)
       call LDT_verify(rc,"Forcing variables list file: not defined")


    !- Initialize precipitation climatology arrays/routines:
!       call LDT_init_pcpclimo()

    !- Read in each meteorological forcing grid info and arrays:

     ! Initialize all config-file specified metforcing datasets:
       do m = 1, LDT_rc%nmetforc
          write(LDT_logunit,fmt='(a22,i2,a1,a20)')" MetForcing Selected: ",&
                m,",  ",trim(LDT_rc%metforc(m))
          call initmetforc(trim(LDT_rc%metforc(m))//char(0),m)
       enddo

       nensem = 0
       do m=1,LDT_rc%nmetforc
          nensem = nensem + LDT_rc%met_nensem(m)
       enddo

    !- Create metforcing names and arrays:
       call create_forcing_structures()

    !- Perturb selected forcing variables:
!       call forcingPerturbSetup()

    endif

  end subroutine LDT_metforcingInit


!BOP
! !ROUTINE: create_forcing_structures
! \label{create_forcing_structures}
!
! !INTERFACE:
  subroutine create_forcing_structures()
! !USES:

!
! !DESCRIPTION:
! This subroutine creates and allocates the memory for the forcing 
! variables for an offline run. The ESMF object for 
! forcing state is created with placeholders for all
! meteorological forcing variables
! 
!EOP
    integer              :: n,m,i
    character*100        :: temp
    character*1          :: nestid(2)
    integer              :: tnvars
    logical              :: file_exists
    integer              :: status
    type(ESMF_Config)    :: forcConfig

    LDT_rc%nf = 0
    do i=1, LDT_rc%nmetforc
       if(LDT_rc%nf < LDT_rc%met_nf(i)) then 
          LDT_rc%nf = LDT_rc%met_nf(i)
       endif
    enddo

    allocate(LDT_FORC_State(LDT_rc%nnest))
    allocate(LDT_FORC_Base_State(LDT_rc%nnest,LDT_rc%nmetforc))

    inquire(file=LDT_rc%forcvarlistfile,exist=file_exists) 
    
    if(.not.file_exists) then 
       write(LDT_logunit,*) "Forcing Variables list file, "
       write(LDT_logunit,*) trim(LDT_rc%forcvarlistfile),", does not exist ..."
       write(LDT_logunit,*) "Program stopping ..."
       call LDT_endrun()
    endif

    write(LDT_logunit,*) "Opening Forcing Variables list file: ",&
         trim(LDT_rc%forcvarlistfile)
    
    forcConfig = ESMF_ConfigCreate(rc=status)
    
    call ESMF_ConfigLoadFile(forcConfig,trim(LDT_rc%forcvarlistfile),rc=status)
  
    tnvars = 0 

    call ESMF_ConfigFindLabel(forcConfig,"Tair:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Tair,&
         "Near Surface Air Temperature", tnvars, status)

    call ESMF_ConfigFindLabel(forcConfig,"Qair:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Qair,&
         "Near Surface Specific Humidity", tnvars,status)
   
    call ESMF_ConfigFindLabel(forcConfig,"SWdown:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_SWdown,&
         "Incident Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWdirect:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_SWdirect,&
         "Incident Direct Surface Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWdiffuse:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_SWdiffuse,&
         "Incident Diffuse Surface Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"LWdown:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_LWdown,&
         "Incident Longwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind_E:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Wind_E,&
         "Eastward Wind", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind_N:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Wind_N,&
         "Northward Wind", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Psurf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Psurf,&
         "Surface Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Rainf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Rainf,&
         "Rainfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Snowf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Snowf,&
         "Snowfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CRainf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_CRainf,&
         "Convective Rainfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"TotalPrecip:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_TotalPrecip,&
         "Total Precipitation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Forc_Hgt:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Forc_Hgt,&
         "Height of Forcing Variables", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Ch:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Ch,&
         "Surface Exchange Coefficient for Heat", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Cm:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Cm,&
         "Surface Exchange Coefficient for Momentum", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Q2sat:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Q2sat,&
         "Saturated Mixing Ratio", tnvars,status)


    call ESMF_ConfigFindLabel(forcConfig,"Emiss:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Emiss,&
         "Surface Emissivity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Cosz:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Cosz,&
         "Cosine of Zenith Angle", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Albedo:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Alb,&
         "Surface Albedo", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"PARDR:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Pardr,&
         "Surface downward PAR direct flux", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"PARDF:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Pardf,&
         "Surface downward PAR diffuse flux", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWnet:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_SWnet,&
         "Net downward shortwave flux", tnvars,status)

! Begin for FEWSNET
    call ESMF_ConfigFindLabel(forcConfig,"PET:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_PET,&
         "Potential ET", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"RefET:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_RefET,&
         "Reference ET", tnvars,status)
! End for FEWSNET

! CAPE available from NLDAS-2
    call ESMF_ConfigFindLabel(forcConfig,"CAPE:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_CAPE,&
         "Convective Available Potential Energy", tnvars,status)

! for CRTM
    call ESMF_ConfigFindLabel(forcConfig,"LPressure:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_lpressure,&
         "Level Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"O3:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_o3,&
         "Ozone Concentration", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Xice:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_XICE,&
         "Sea Ice Mask", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"QSFC:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_QSFC,&
         "Surface Specific Humidity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CHS2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_CHS2,&
         "2m Surface Exchange Coefficient for Heat", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CQS2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_CQS2,&
         "2m Surface Exchange Coefficient for Moisture", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"T2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_T2,&
         "2m Air Temperature", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Q2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_Q2,&
         "2m Specific Humidity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"TH2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_TH2,&
         "2m Potential Temperature", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"TMN:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_TMN,&
         "Soil Temperature at Lower Boundary", tnvars,status)

!<for vic>
    call ESMF_ConfigFindLabel(forcConfig,"Snowflag:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_SNOWFLAG,&
         "Snowflag", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Density:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_DENSITY,&
         "Atmospheric Density", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"VaporPress:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_VAPORPRESS,&
         "Vapor Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"VaporPressDeficit:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_VAPORPRESSDEFICIT,&
         "Vapor Pressure Deficit", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind:",rc=status)
    call get_forcingvar_attributes(forcConfig,LDT_FORC_WIND,&
         "Wind Speed", tnvars,status)
!</for vic>

    if(tnvars.gt.LDT_rc%nf) then 
       LDT_rc%nf = tnvars
    endif

    do n=1,LDT_rc%nnest       
       write(unit=temp,fmt='(i2.2)') n
       read(unit=temp,fmt='(2a1)') nestid
       
       LDT_FORC_State(n) = ESMF_StateCreate(name=&
            "Forcing State"//nestid(1)//nestid(2),&
            rc=status)
       call LDT_verify(status, &
            'error in ESMF_StateCreate:LDT_FORC_State in create_forcing_structures')       

       do m=1,LDT_rc%nmetforc
          LDT_FORC_Base_State(n,m) = ESMF_StateCreate(name=&
               "Forcing State"//nestid(1)//nestid(2),&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateCreate:LDT_FORC_BaseState in create_forcing_structures')       
       enddo

       call add_forcing_fields(n,"Tair",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Tair)
       call add_forcing_fields(n,"Qair",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Qair)
       call add_forcing_fields(n,"SWdown",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_SWdown)
       call add_forcing_fields(n,"SWdirect",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_SWdirect)
       call add_forcing_fields(n,"SWdiffuse",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_SWdiffuse)
       call add_forcing_fields(n,"LWdown",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_LWdown)
       call add_forcing_fields(n,"Uwind",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Wind_E)
       call add_forcing_fields(n,"Vwind",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Wind_N)
       call add_forcing_fields(n,"Psurf",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Psurf)
       call add_forcing_fields(n,"Rainf",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Rainf)
       call add_forcing_fields(n,"Snowf",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Snowf)
       call add_forcing_fields(n,"CRainf",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_CRainf)
       call add_forcing_fields(n,"TotalPrecip",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_TotalPrecip)

       call add_forcing_fields(n,"Forc_Hgt",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Forc_Hgt)
       call add_forcing_fields(n,"CH",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Ch)
       call add_forcing_fields(n,"CM",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Cm)
       call add_forcing_fields(n,"Q2SAT",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Q2sat)
       call add_forcing_fields(n,"EMISS",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Emiss)
       call add_forcing_fields(n,"COSZ",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Cosz)
       call add_forcing_fields(n,"Albedo",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Alb)
       call add_forcing_fields(n,"PARDR",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Pardr)
       call add_forcing_fields(n,"PARDF",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Pardf)
       call add_forcing_fields(n,"SWnet",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_SWnet)

       call add_forcing_fields(n,"XICE",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Xice)

       call add_forcing_fields(n,"QSFC",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_QSFC)
       call add_forcing_fields(n,"CHS2",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_CHS2)
       call add_forcing_fields(n,"CQS2",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_CQS2)
       call add_forcing_fields(n,"T2",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_T2)
       call add_forcing_fields(n,"Q2",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_Q2)
       call add_forcing_fields(n,"TH2",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_TH2)
       call add_forcing_fields(n,"TMN",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_TMN)

       call add_forcing_fields(n,"LPRESSURE",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_lpressure)
       call add_forcing_fields(n,"O3",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_forc_o3)

!<for vic>
       call add_forcing_fields(n,"SNOWFLAG",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_SNOWFLAG)
       call add_forcing_fields(n,"DENSITY",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_DENSITY)
       call add_forcing_fields(n,"VAPORPRESS",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_VAPORPRESS)
       call add_forcing_fields(n,"VAPORPRESSDEFICIT",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_VAPORPRESSDEFICIT)
       call add_forcing_fields(n,"WIND",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_WIND)
!</for vic>
       call add_forcing_fields(n,"PET",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_PET)
       call add_forcing_fields(n,"REFET",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_RefET)
       call add_forcing_fields(n,"CAPE",LDT_FORC_State(n),&
            LDT_FORC_Base_State(n,:),LDT_FORC_CAPE)
    enddo
    
    call ESMF_ConfigDestroy(forcConfig)

    do n=1,LDT_rc%nnest      
       do m=1,LDT_rc%nmetforc

          if(LDT_rc%met_ecor(m).ne."none") then 
             allocate(LDT_forc(n,m)%modelelev(LDT_rc%ngrid(n)))
          endif
          if(LDT_rc%met_nf(m).gt.0) then 
             allocate(LDT_forc(n,m)%metdata1(LDT_rc%met_nf(m),LDT_rc%ngrid(n)))
             allocate(LDT_forc(n,m)%metdata2(LDT_rc%met_nf(m),LDT_rc%ngrid(n)))
             if(LDT_rc%met_tinterp(m).eq."trilinear") then 
                allocate(LDT_forc(n,m)%metdata3(LDT_rc%met_nf(m),&
                     LDT_rc%ngrid(n)))
                LDT_forc(n,m)%metdata3 = 0.0
             endif

             LDT_forc(n,m)%metdata1 = 0.0
             LDT_forc(n,m)%metdata2 = 0.0

          endif
       enddo
    enddo

  end subroutine create_forcing_structures

#if 0
!BOP
! 
! !ROUTINE: forcingPerturbSetup
! \label{forcingPerturbSetup}
! 
! !INTERFACE: 
  subroutine forcingPerturbSetup
! !USES: 
    use LDT_perturbMod 

! !DESCRIPTION: 
!  
!  This routine sets up and allocates structures required 
!  for forcing pertubations. The ESMF state objects for perturbed
!  forcing fields and the perturbations themselves are created. 
!  The perturbation related attributed are also assigned in this
!  routine. 
!EOP
    integer                 :: n 
    character*1             :: nestid(2)
    character*100           :: temp
    integer                 :: i,j
    integer                 :: status
    character*100, allocatable  :: pertobjs(:)
    integer,       allocatable  :: order(:)
    real   ,       allocatable  :: ccorr(:,:)
    character*40, allocatable   :: varname(:)    
    real   , allocatable        :: varmax(:)
    real   , allocatable        :: varmin(:)
    real   , allocatable        :: ssdev(:)
    type(pert_dec_type)     :: forc_pert
    type(ESMF_Field)        :: pertField
    type(ESMF_Field)        :: varField(LDT_rc%nf)
    type(ESMF_Field)        :: varpertField(LDT_rc%nf)
    type(ESMF_ArraySpec)    :: arrspec1
    integer                 :: ftn 

    if(LDT_rc%perturb_forcing.ne."none") then 
       allocate(LDT_FORC_Pert_State(LDT_rc%nnest))

       do n=1,LDT_rc%nnest
          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid
          
          LDT_FORC_Pert_State(n) = ESMF_StateCreate(name=&
               "Forcing Perturbations"//&
               nestid(1)//nestid(2),  rc=status)
          call LDT_verify(status, 'force pert state create')
          write(LDT_logunit,*) 'Successfully created forcing perturbations states ..'
       enddo

       do n=1,LDT_rc%nnest
          write(LDT_logunit,*) 'Opening Forcing Attributes file ',&
               LDT_rc%forcattribfile
          
          ftn = LDT_getNextUnitNumber()
          open(ftn,file=trim(LDT_rc%forcattribfile), status='old')
          read(ftn,*)
          read(ftn,*) LDT_rc%nforcepert
          read(ftn,*)
          allocate(varname(LDT_rc%nforcepert))
          allocate(varmax(LDT_rc%nforcepert))
          allocate(varmin(LDT_rc%nforcepert))

          allocate(forc_pert%vname(LDT_rc%nforcepert))
          allocate(forc_pert%perttype(LDT_rc%nforcepert))
          allocate(forc_pert%ssdev(LDT_rc%nforcepert))
          allocate(forc_pert%stdmax(LDT_rc%nforcepert))
          allocate(forc_pert%zeromean(LDT_rc%nforcepert))
          allocate(forc_pert%tcorr(LDT_rc%nforcepert))
          allocate(forc_pert%xcorr(LDT_rc%nforcepert))
          allocate(forc_pert%ycorr(LDT_rc%nforcepert))
          allocate(forc_pert%ccorr(LDT_rc%nforcepert,LDT_rc%nforcepert))
          
          do i=1,LDT_rc%nforcepert
             read(ftn,fmt='(a40)') varname(i)
             read(ftn,*) varmin(i),varmax(i)
             write(LDT_logunit,*) varname(i),varmin(i),varmax(i)
          enddo
          call LDT_releaseUnitNumber(ftn)

          call LDT_readPertAttributes(LDT_rc%nforcepert, &
               LDT_rc%forcpertattribfile,forc_pert)
      
          call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
               rc=status)
          call LDT_verify(status,'arrspecset in metforcingMod')
             
          do i=1,LDT_rc%nforcepert
             
             call ESMF_StateGet(LDT_FORC_state(n),&
                  trim(varname(i)), varField(i),rc=status)
             call LDT_verify(status,'stateget in metforcingMod')
             call ESMF_AttributeSet(varField(i),"Max Value",&
                  varmax(i),rc=status)
             call LDT_verify(status,'attributeset(max) in metforcingMod')
             call ESMF_AttributeSet(varField(i),"Min Value",&
                  varmin(i),rc=status)
             call LDT_verify(status,'attributeset(min) in metforcingMod')    
                
             varpertField(i) = ESMF_FieldCreate(grid=LDT_vecTile(n),&
                  arrayspec=arrspec1,&
                  name=trim(varname(i)),rc=status)
             call LDT_verify(status,'fieldcreate in metforcingMod')
                          
             call ESMF_StateAdd(LDT_FORC_Pert_State(n),(/varpertField(i)/),&
                  rc=status)
             call LDT_verify(status, &
                  'error in ESMF_StateAdd in perturbSetup')
          end do
          deallocate(varname)
          deallocate(varmax)
          deallocate(varmin)

          allocate(pertobjs(LDT_rc%nforcepert))
          allocate(order(LDT_rc%nforcepert))
          allocate(ccorr(LDT_rc%nforcepert,LDT_rc%nforcepert))

          call ESMF_StateGet(LDT_FORC_Pert_State(n), &
               itemNameList=pertobjs, rc=status)
          call LDT_verify(status, &
               'error in ESMF_StateGet in perturbSetup')

          order = -1
          do i=1,LDT_rc%nforcepert
             do j=1,LDT_rc%nforcepert
                if(forc_pert%vname(j).eq.pertobjs(i)) then 
                   order(i) = j
                   exit;
                endif
             enddo             
          enddo

          do i=1,LDT_rc%nforcepert             
             do j=1,LDT_rc%nforcepert
                ccorr(i,j) = forc_pert%ccorr(order(i),order(j))
             enddo
          enddo

          do i=1,LDT_rc%nforcepert             

             call ESMF_StateGet(LDT_FORC_Pert_State(n), &
                  pertobjs(i),pertField,&
                  rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,"Perturbation Type",&
                  forc_pert%perttype(i), rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')

             if(LDT_rc%ngrid(n).gt.0) then
                ssdev = forc_pert%ssdev(order(i))

                call ESMF_AttributeSet(pertField,"Standard Deviation",&
                     ssdev,itemCount=LDT_rc%ngrid(n),rc=status)
                call LDT_verify(status,&
                     'error in ESMF_AttributeSet in forcingPerturbSetup')
             endif

             call ESMF_AttributeSet(pertField,"Std Normal Max",&
                  forc_pert%stdmax(i),rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')

             call ESMF_AttributeSet(pertField,"Ensure Zero Mean",&
                  forc_pert%zeromean(i),rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,&
                  "Temporal Correlation Scale",&
                  forc_pert%tcorr(i),rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,"X Correlation Scale",&
                  forc_pert%xcorr(i),rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,"Y Correlation Scale",&
                  forc_pert%ycorr(i),rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,&
                  "Cross Correlation Strength",&
                  ccorr(i,:),&
                  itemCount=LDT_rc%nforcepert,&
                  rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
          enddo
          deallocate(pertobjs)
          deallocate(order)
          deallocate(ccorr)
          deallocate(ssdev)
       enddo

       call perturbinit(trim(LDT_rc%perturb_forcing)//char(0), 1)
       call perturbsetup(trim(LDT_rc%perturb_forcing)//char(0),&
            1, 1, LDT_FORC_State,&
            LDT_FORC_Pert_State)

    endif
  end subroutine forcingPerturbSetup
#endif


!BOP
! !ROUTINE: LDT_get_met_forcing
! \label{LDT_get_met_forcing}
!
! !INTERFACE:
  subroutine LDT_get_met_forcing(n)
! !USES:

! !ARGUMENTS:
    integer,intent(in) :: n

! !DESCRIPTION:
! 
! This routine issues the calls to retrieve the forcing variables from
! the specific forcing scheme. The retrieval call is expected to read
! the data, perform any spatial interpolation and other transformations
! needed to grid the data to the running domain and resolution. 
! This invocation is followed by the call to temporally interpolate
! the data to the current model timestep. The temporal interpolation 
! is performed based on the two consecutive forcing analyses. 
! A zenith angle based interpolation is performed for the temporal 
! disaggregation of downward shortwave radiation. Finally any 
! topographic corrections and downscaling to the variables are 
! applied. 
! 
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!  \end{description}
! 
! The methods invoked are:  
! \begin{description}
!  \item[retrieveforcing](\ref{retrieveforcing}) \newline
!    invokes the generic method in the registry to retrieve the 
!    met forcing data
!  \item[timeinterp](\ref{timeinterp}) \newline
!    invokes the generic method in the registry to perform
!    temporal interpolation
!  \item[LDT\_lapseRateCorrection](\ref{LDT_lapseRateCorrection}) \newline
!    method to apply lapse-rate based topographical corrections
!  \item[LDT\_microMetCorrection](\ref{LDT_microMetCorrection}) \newline
!    method to apply topographical corrections (not currently supported)
! \end{description}
!EOP

    integer :: m

    if(LDT_rc%nmetforc.gt.0) then 
       do m=1,LDT_rc%nmetforc

          call retrievemetforc(trim(LDT_rc%metforc(m))//char(0),n,m)

          call timeinterpmetforc(trim(LDT_rc%metforc(m))//char(0),n,m)

        ! Perform topographic corrections and/or downscaling:
          if(LDT_rc%met_ecor(m) .eq. "lapse-rate") then
             call LDT_lapseRateCorrection(n,LDT_forc(n,m)%modelelev,&
                  LDT_FORC_Base_State(n,m))
#if 0
          elseif(LDT_rc%met_ecor(m).eq."slope-aspect") then
             call LDT_slopeAspectCorrection(n,LDT_FORC_Base_State(n,m))
          elseif(LDT_rc%met_ecor(m).eq."lapse-rate and slope-aspect") then
             call LDT_lapseRateCorrection(n, LDT_forc(n,m)%modelelev,&
                  LDT_FORC_Base_State(n,m))
             call LDT_slopeAspectCorrection(n, LDT_FORC_Base_State(n,m))
#endif
          end if
       enddo

     ! Blending algorithms (overlay, forcing ensembles, bias correction..)
       if( LDT_rc%metforc_blend_alg.eq."overlay" ) then ! simple overlays
          call overlayForcings(n)
       elseif(LDT_rc%metforc_blend_alg.eq."ensemble") then !forcing ensembles
          call ensembleForcings(n)
       endif

    endif

  end subroutine LDT_get_met_forcing

!BOP
! !ROUTINE: overlayForcings
! \label{overlayForcings}
!
! !INTERFACE: 
  subroutine overlayForcings(n)

    implicit none
! !ARGUMENTS:     
    integer, intent(in)    :: n 
!
! !DESCRIPTION:
!  This subroutine blends different forcing inputs by sequentially 
!  overlaying them on top of each other in the order in which they 
!  are specified. 
!
!EOP
    integer                :: fobjcount
    integer                :: i,t,m
    type(ESMF_Field)       :: mrgField, baseField
    integer                :: status, status1, status2
    real, pointer          :: forcdata_base(:), forcdata_mrg(:)
    character*100, allocatable :: forcobjs(:)
    
    call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed for objcount in overlayForcings')
    
    allocate(forcobjs(fobjcount))
    
    call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed for forcobjs in overlayForcings')
    
    do i=1,fobjcount
       call ESMF_StateGet(LDT_FORC_State(n),forcobjs(i),mrgField,&
            rc=status)
       call LDT_verify(status, 'ESMF_StateGet failed for '//trim(forcobjs(i)))
       call ESMF_FieldGet(mrgField,localDE=0,farrayPtr=forcdata_mrg, &
            rc=status)
       call LDT_verify(status,'ESMF_FieldGet failed for forcdata_mrg in overlayForcings')
       
    !- Loop over stored metforcing dataset types:
       do m=1,LDT_rc%nmetforc
          call ESMF_StateGet(LDT_FORC_Base_State(n,m),forcobjs(i),baseField,&
               rc=status1)

          if( status1.eq.0 ) then 
             call ESMF_FieldGet(baseField,localDE=0,farrayPtr=forcdata_base, &
                  rc=status2)
             call LDT_verify(status2,'ESMF_FieldGet (basefield) failed in overlayForcings')
             
             do t=1,LDT_rc%ntiles(n)
                ! Ignore and not merge values that are undefined:
                if(forcdata_base(t).ne.-9999.0) then 
                   forcdata_mrg(t) = forcdata_base(t)
                endif
             enddo
          endif
       enddo
       
    enddo

    deallocate(forcobjs)

  end subroutine overlayForcings

!BOP
! !ROUTINE: ensembleForcings
! \label{ensembleForcings}
!
! !INTERFACE: 
  subroutine ensembleForcings(n)

    implicit none
! !ARGUMENTS:     
    integer, intent(in)    :: n 
!
! !DESCRIPTION:
!  This subroutine blends different forcing inputs by assigning them 
!  to different ensemble members. 
!
!EOP
    integer                :: fobjcount
    integer                :: i,t,m,tid,tid1,tid2,index1
    type(ESMF_Field)       :: mrgField, baseField
    integer                :: status, status1, status2
    real, pointer          :: forcdata_base(:), forcdata_mrg(:)
    character*100, allocatable :: forcobjs(:)
    
    call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed for in ensembleForcings')
    
    allocate(forcobjs(fobjcount))

    call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
    call LDT_verify(status,'ESMF_StateGet failed in ensembleForcings')
    
    do i=1,fobjcount
       call ESMF_StateGet(LDT_FORC_State(n),forcobjs(i),mrgField,&
            rc=status)
       call LDT_verify(status, 'ESMF_StateGet failed for '//trim(forcobjs(i)))
       call ESMF_FieldGet(mrgField,localDE=0,farrayPtr= forcdata_mrg, &
            rc=status)
       call LDT_verify(status,'ESMF_FieldGet failed for in ensembleForcings')
       
       do m=1,LDT_rc%nmetforc
          call ESMF_StateGet(LDT_FORC_Base_State(n,m),forcobjs(i),baseField,&
               rc=status1)
          if(status1.eq.0) then 
             call ESMF_FieldGet(baseField,localDE=0,farrayPtr= forcdata_base, &
                  rc=status2)
             call LDT_verify(status2,&
                  'ESMF_FieldGet failed in LDT_metForcingMod')
             if(m.ne.1) then      
                do t=1, LDT_rc%ntiles(n)/LDT_rc%nensem(n)
                   tid1 =(t-1)*LDT_rc%nensem(n)+(m-1)*LDT_rc%nperforc+ 1
                   tid2 = tid1 + ( LDT_rc%nperforc-1 )
                   do tid=tid1,tid2
!                      index1 =LDT_domain(n)%tile(tid)%index     
!                      if(forcdata_base(index1).ne.LDT_rc%udef) then 
                      if(forcdata_base(tid).ne.LDT_rc%udef) then 
                         forcdata_mrg(tid)=forcdata_base(tid)
                      endif
                   enddo
                enddo
             else
                forcdata_mrg = forcdata_base
             endif
          endif
       enddo
       
    enddo

    deallocate(forcobjs)

  end subroutine ensembleForcings


!BOP
! !ROUTINE: LDT_perturb_forcing
! \label{LDT_perturb_forcing}
! 
! !INTERFACE: 
  subroutine LDT_perturb_forcing(n)
! !USES: 

! !ARGUMENTS: 
    integer,  intent(IN)   :: n 
! 
! !DESCRIPTION: 
!  This routine computes the forcing perturbations and then applies them to 
!   the forcing fields.  This routine also invokes the diagnose routine to map 
!   the forcing fields to the history writer. 
! 
!  \begin{description}
!  \item[diagnoseForcingOutput](\ref{diagnoseForcingOutput}) \newline
!    routine to map a forcing fields to the LDT history writer
!  \end{description}
!EOP
    integer                   :: k 
    real                      :: curr_time
    integer                   :: status
    integer                   :: i, j,t,m,t_unpert
    integer                   :: objcount
    integer                   :: fobjcount
    integer                   :: offset1
    character(len=LDT_CONST_PATH_LEN)             :: fname
    real,             pointer     :: forcvar(:)    
    real,             pointer     :: pertdata(:)
    character*100,    allocatable :: forcobjs(:)
    character*100,    allocatable :: stateobjs(:)
    type(ESMF_Field), allocatable :: forcField(:)
    type(ESMF_Field), allocatable :: stateField(:)
    integer                   :: perttype
    real                      :: maxval
    real                      :: minval
    real                      :: delta
    
#if 0 
    k = 1
    if(LDT_rc%perturb_forcing.ne."none") then 
       curr_time = float(LDT_rc%hr)*3600+60*float(LDT_rc%mn)+float(LDT_rc%ss)
       if(mod(curr_time,real(LDT_rc%pertforcinterval)).eq.0)then
          
          call perturbmethod(trim(LDT_rc%perturb_forcing)//char(0),1,n,k, &
               LDT_FORC_State(n), LDT_FORC_Pert_State(n))
          
          call ESMF_StateGet(LDT_FORC_Pert_State(n),itemCount=objcount,rc=status)
          call LDT_verify(status, &
               'error in ESMF_StateGet in LDT_perturb_forcing')
          
          call ESMF_StateGet(LDT_FORC_State(n),itemCount=fobjcount,rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateGet in LDT_perturb_forcing')
          
          allocate(stateobjs(objcount))
          allocate(stateField(objcount))
          allocate(forcobjs(fobjcount))
          allocate(forcField(objcount))
          
          call ESMF_StateGet(LDT_FORC_Pert_State(n),itemNameList=stateobjs,rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateGet in LDT_perturb_forcing')
          
          call ESMF_StateGet(LDT_FORC_State(n),itemNameList=forcobjs,rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateGet in LDT_perturb_forcing')
          
          do i=1,objcount
             call ESMF_StateGet(LDT_FORC_Pert_State(n),stateobjs(i),stateField(i),&
                  rc=status)
             call ESMF_FieldGet(stateField(i),localDE=0,farrayPtr= pertdata, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet in LDT_perturb_forcing')
             
             call ESMF_AttributeGet(stateField(i),"Perturbation Type",&
                  perttype,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeGet in LDT_perturb_forcing')
             
             fname = stateobjs(i)

             call ESMF_StateGet(LDT_FORC_State(n),trim(fname),forcField(i),rc=status)
             if(status.eq.0) then 
                
                call ESMF_FieldGet(forcField(i),localDE=0,farrayPtr= forcvar,rc=status)
                call LDT_verify(status,&
                     'error in ESMF_FieldGet in LDT_perturb_forcing')
                call ESMF_AttributeGet(forcField(i), "Max Value",maxval,rc=status)
                call LDT_verify(status,&
                     'error in ESMF_AttributeGet in LDT_perturb_forcing')
                
                call ESMF_AttributeGet(forcField(i), "Min Value",minval,rc=status)
                call LDT_verify(status,&
                     'error in ESMF_AttributeGet in LDT_perturb_forcing')
                
#if 0
                do t=1,LDT_rc%ntiles(n)
                   if(perttype.eq.0) then 
                      forcvar(t) = forcvar(t)+ pertdata(t)
                   elseif(perttype.eq.1) then 
                      forcvar(t) = forcvar(t)*pertdata(t)
                   endif
                   if(forcvar(t).lt.minval) forcvar(t) = minval
                   if(forcvar(t).gt.maxval) forcvar(t) = maxval                 
                enddo
#endif
                do j=1,LDT_rc%ntiles(n)/LDT_rc%nensem(n)
                   do m=1,LDT_rc%nensem(n)
                      t = (j-1)*LDT_rc%nensem(n)+m
                      if(perttype.eq.0) then
                         if(LDT_rc%pert_bias_corr.eq.1) then
                            if(m.ne.LDT_rc%nensem(n)) then
                               forcvar(t) = forcvar(t)+ pertdata(t)
                            endif
                         else
                            forcvar(t) = forcvar(t)+ pertdata(t)
                         endif
                      elseif(perttype.eq.1) then
                         if(LDT_rc%pert_bias_corr.eq.1) then
                            if(m.ne.LDT_rc%nensem(n)) then
                               forcvar(t) = forcvar(t)*pertdata(t)
                            endif
                         else
                            forcvar(t) = forcvar(t)*pertdata(t)
                         endif
                      endif
                   enddo
                enddo
                if(LDT_rc%pert_bias_corr.eq.1) then
                   do j=1,LDT_rc%ntiles(n)/LDT_rc%nensem(n)
                      t_unpert = j*LDT_rc%nensem(n)
                      delta = 0.0
                      do m=1,LDT_rc%nensem(n)-1
                         t = (j-1)*LDT_rc%nensem(n)+m
                         if(m.ne.LDT_rc%nensem(n)) then
                            delta = delta + (forcvar(t)-forcvar(t_unpert))
                         endif
                      enddo

                      delta = delta/(LDT_rc%nensem(n)-1)
                      do m=1,LDT_rc%nensem(n)-1
                         t = (j-1)*LDT_rc%nensem(n)+m
                         forcvar(t) = forcvar(t) - delta
                      enddo
                   enddo
                endif
                do t=1,LDT_rc%ntiles(n)
                   if(forcvar(t).lt.minval) forcvar(t) = minval
                   if(forcvar(t).gt.maxval) forcvar(t) = maxval
                enddo

             endif
          enddo
          
          deallocate(forcField)
          deallocate(forcobjs)
          deallocate(stateobjs)
          deallocate(stateField)
       endif
    endif

#endif

 !- Diagnose calls to map the forcing variables to the 
 !  history output writer:
    call diagnoseForcingOutput(n)

 end subroutine LDT_perturb_forcing

!BOP
! !ROUTINE: LDT_metforcing_reset
!  \label{LDT_metforcing_reset}
! 
! !INTERFACE: 
  subroutine LDT_metforcing_reset
! !USES: 
! 
! !DESCRIPTION: 
!  This subroutine issues the invocation to reset variables for 
!  specific forcing scheme. 
!  Involves the C function-table, resetmetforc, which is found in
!   core$/LDT_metforcing_FTable.c$.
!
!EOP
    integer  :: m

    ! Call the C function-table to reset each forcing dataset variables/parameters,
    !  where needed.
    if(LDT_rc%nmetforc.gt.0) then 
       do m=1,LDT_rc%nmetforc
          call resetmetforc(trim(LDT_rc%metforc(m))//char(0),m)
       enddo
    endif

  end subroutine LDT_metforcing_reset

!BOP
! !ROUTINE: LDT_metforcing_finalize
! \label{LDT_metforcing_finalize}
! 
! !INTERFACE: 
  subroutine LDT_metforcing_finalize
! !USES:

!
! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures in the specific instance of the 
!  forcing scheme. 
! 
! The methods invoked are:  
! \begin{description}
!  \item[forcingfinalize](\ref{forcingfinalize}) \newline
!    invokes the generic method in the registry to cleanup 
!    the allocated structures for the met forcing scheme. 
! \end{description}
!EOP
    integer  :: m

    print *, " I AM HERE (inside core/LDT_metforcingMod)!"      ! KRA - 08/17/2015

    if(LDT_rc%nmetforc.gt.0) then 
       do m=1,LDT_rc%nmetforc
          call finalmetforc(trim(LDT_rc%metforc(m))//char(0),m)
          print *, m, LDT_rc%nmetforc, trim(LDT_rc%metforc(m))  ! KRA - 08/17/2015
 stop
       enddo
       deallocate(LDT_FORC_State)
#if 0          
       if(LDT_rc%perturb_forcing.ne."none") then 
          deallocate(LDT_FORC_Pert_State)
       endif
#endif
    endif

  end subroutine LDT_metforcing_finalize

!BOP
! 
!  !ROUTINE: get_forcingvar_attributes
!  \label{get_forcingvar_attributes}
! 
! !INTERFACE: 
  subroutine get_forcingvar_attributes(forcConfig, forc_attrib, &
       forcname, varcount, checkentry)

! !USES:     

    implicit none
! !ARGUMENTS:    
    type(ESMF_Config)      :: forcConfig
    type(LDT_forcDataEntry) :: forc_attrib 
    character(len=*)       :: forcname
    integer                :: varcount
    integer                :: checkentry
! 
! !DESCRIPTION: 
! 
!  This subroutine reads the configurable options for each forcing variable. 
!  If the entry for the variable is not present, then the variable is set to be
!  not chosen
! 
!EOP
    character(len=3)    :: fid
    integer             :: k 
    integer             :: status
    integer             :: n 
    
    n = 1
    if(checkentry.eq.0) then ! the forcing field specification is found 
       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%selectOpt, &
            default=0,rc=status)
       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%vlevels,&
            default=0,rc=status)
       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%units, rc=status)

       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%timeAvgOpt, rc=status)
       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%selectProc, rc=status)
       
       if(forc_attrib%selectOpt.eq.1) then 
          allocate(forc_attrib%varname(forc_attrib%vlevels))
          do k=1,forc_attrib%vlevels
             write(fid,fmt='(i3.3)') k
             forc_attrib%varname(k) = trim(forcname)//' Level '//trim(fid)
          enddo
          write(LDT_logunit,*) "FORCING: ",forcname, forc_attrib%vlevels, &
               forc_attrib%units
          varcount = varcount + forc_attrib%vlevels
       endif       
       if(forc_attrib%timeAvgOpt.eq.2) then 
          allocate(forc_attrib%modelOutput(2,LDT_rc%ntiles(n),&
               forc_attrib%vlevels))
       else
          allocate(forc_attrib%modelOutput(1,LDT_rc%ntiles(n),&
               forc_attrib%vlevels))
       endif
       allocate(forc_attrib%count(LDT_rc%max_model_types,forc_attrib%vlevels))
       forc_attrib%modelOutput = 0
       forc_attrib%count = 0
    else
       forc_attrib%selectOpt = 0 
       forc_attrib%vlevels = 1
       allocate(forc_attrib%varname(forc_attrib%vlevels))
       forc_attrib%varname(1) = 'NOT_SPECIFIED'
    endif

  end subroutine get_forcingvar_attributes


!BOP
! 
! !ROUTINE: add_forcing_fields
! \label{add_forcing_fields}
! 
! !INTERFACE: 
  subroutine add_forcing_fields(n, name, FORC_State, FORC_Base_State,forc_attrib)
! !USES:     

    implicit none
! !ARGUMENTS: 
    integer                   :: n 
    character(len=*)          :: name 
    type(ESMF_State)          :: FORC_State
    type(ESMF_State)          :: FORC_Base_State(LDT_rc%nmetforc)
    type(LDT_forcDataEntry)   :: forc_attrib
! 
! !DESCRIPTION: 
!  This subroutine creates fields for the specified forcing variables 
!  and initializes them to be undefined values. 
!
!EOP

    integer                   :: k,m
    type(ESMF_Field)          :: varField_b
    type(ESMF_Field)          :: varField_m
    type(ESMF_ArraySpec)      :: arrspec1
    real,          pointer    :: var(:)
    integer                   :: status

    forc_attrib%short_name = name
    forc_attrib%standard_name = name

    if(forc_attrib%selectOpt.eq.1) then 
       call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LDT_verify(status,&
            'error in ESMF_ArraySpecSet in add_forcing_fields')       

       do k=1,forc_attrib%vlevels
          varField_m = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = trim(forc_attrib%varname(k)),&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')
          
          call ESMF_AttributeSet(varField_m,"Units",trim(forc_attrib%units),&
               rc=status)
          call LDT_verify(status,&
               'error in ESMF_AttributeSet in add_forcing_fields')
          
          call ESMF_FieldGet(varField_m,localDE=0,farrayPtr=var,rc=status)
          call LDT_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = LDT_rc%udef

          call ESMF_StateAdd(FORC_State,(/varField_m/),rc=status)
          call LDT_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

          do m=1,LDT_rc%nmetforc

             varField_b = ESMF_FieldCreate(grid=LDT_vecTile(n), &
               arrayspec=arrspec1, name = trim(forc_attrib%varname(k)),&
               rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldCreate in add_forcing_fields')
             
             call ESMF_AttributeSet(varField_b,"Enabled",0,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in add_forcing_fields')          

             call ESMF_AttributeSet(varField_b,"Units",trim(forc_attrib%units),&
                  rc=status)
             call LDT_verify(status,&
                  'error in ESMF_AttributeSet in add_forcing_fields')
             
             call ESMF_FieldGet(varField_b,localDE=0,farrayPtr=var,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet in add_forcing_fields')
             var = LDT_rc%udef
             
             call ESMF_StateAdd(FORC_Base_State(m),(/varField_b/),rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateAdd in add_forcing_fields')
          enddo
       enddo
    endif

  end subroutine add_forcing_fields

!BOP
! 
!  !ROUTINE: diagnoseForcingOutput
!  \label{diagnoseForcingOutput}
! 
! !INTERFACE: 
  subroutine diagnoseForcingOutput(n)

! !USES: 
    use LDT_constantsMod,  only : LDT_MS2KMDAY

! !ARGUMENTS: 
    integer,  intent(IN) :: n 
! 
! !DESCRIPTION: 
!  This routines issues the diagnose calls to map the forcing variables 
!  to the history output writer. 
! 
!   The routines invoked are: 
!  \begin{description}
!  \item[diagnoseForcingVar](\ref{diagnoseForcingVar}) \newline
!    generic routine to map a single variable to the LDT 
!    history writer
!  \end{description}
!
!EOP
    integer            :: status, t
    integer            :: k 
    type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
    type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField,totprecField
    type(ESMF_Field)   :: swdirField,swdifField,hField,chField,cmField
    type(ESMF_Field)   :: emissField,mixField,coszField,albField
    type(ESMF_Field)   :: tempField
    real, pointer      :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
    real, pointer      :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:),totprec(:)
    real, pointer      :: swdir(:),swdif(:),harray(:),charray(:),cmarray(:)
    real, pointer      :: emiss(:),mix(:),cosz(:),alb(:)
    real, pointer      :: tempPtr(:)
    real               :: windmag, windmag2

    k = 1 ! since forcing is same for all DA instances. 
! Currently this is set to write only the surface forcing variables and not 
! the entire profile

    if(LDT_rc%wout.ne."none") then 
       if(LDT_FORC_Tair%selectOpt.eq.1) then
          do k=1,LDT_FORC_Tair%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Tair%varname(k)),&
                  tmpField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Tair in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Tair in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Tair,&
                     t,k,tmp(t),&
                     unit="K",direction="-",&
                     vmin = 213.0, vmax = 333.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Qair%selectOpt.eq.1) then
          do k=1,LDT_FORC_Qair%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Qair%varname(k)),&
                  q2Field,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Qair in diagnoseForcingOutput')
             
             call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Qair in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Qair,&
                     t,k,q2(t),&
                     unit="kg/kg",direction="-",&
                     vmin = 0.0, vmax = 0.03)
             enddo
          enddo
       endif

       if(LDT_FORC_SWdown%selectOpt.eq.1) then
          do k=1,LDT_FORC_SWdown%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_SWdown%varname(k)),&
                  swdField,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:SWdown in diagnoseForcingOutput')

             call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:SWdown in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_SWdown,&
                     t,k,swd(t),&
                     unit="W/m2",direction="DN",&
                     vmin = 0.0, vmax = 1360.0)
             enddo
          enddo
       endif

       if(LDT_FORC_LWdown%selectOpt.eq.1) then
          do k=1,LDT_FORC_LWdown%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_LWdown%varname(k)),&
                  lwdField,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:LWdown in diagnoseForcingOutput')

             call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:LWdown in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_LWdown,&
                     t,k,lwd(t),&
                     unit="W/m2",direction="DN",&
                     vmin = 0.0, vmax = 750.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Wind_E%selectOpt.eq.1.and.LDT_FORC_Wind_N%selectOpt.eq.1) then  
          do k=1,LDT_FORC_Wind_E%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Wind_E%varname(k)),&
                  uField,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Wind_E in diagnoseForcingOutput')

             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Wind_N%varname(k)),&
                  vField,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_stateGet:Wind_N in diagnoseForcingOutput')
             
             call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Wind_E in diagnoseForcingOutput')
             
             call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Wind_N in diagnoseForcingOutput')             
             
             do t=1,LDT_rc%ntiles(n)
                windmag  = LDT_rc%udef
                windmag2 = LDT_rc%udef
              ! Make sure undefined values are not factored in to wind-mag calc:
                if( uwind(t) == LDT_rc%udef .or. vwind(t) == LDT_rc%udef ) then
                   windmag  = LDT_rc%udef
                   windmag2 = LDT_rc%udef
                else
                   windmag  = sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t))
                   windmag2 = (sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t))*LDT_MS2KMDAY)
                endif
                call diagnoseForcingVar(n, LDT_FORC_WIND, t, k, &
!                     sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t)), &
                     windmag, &
                     unit="m/s",direction="-",&
                     vmin = 0.0, vmax = 0.0)
                call diagnoseForcingVar(n, LDT_FORC_WIND, t, k, &
!                     sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t))*LDT_MS2KMDAY, &
                     windmag2, &
                     unit="km/day",direction="-",&
                     vmin = 0.0, vmax = 0.0)
             enddo

          enddo
       endif

       if(LDT_FORC_Psurf%selectOpt.eq.1) then
          do k=1,LDT_FORC_Psurf%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Psurf%varname(k)),&
                  psurfField,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Psurf in diagnoseForcingOutput')

             call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Psurf in diagnoseForcingOutput')
             
             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_PSURF,&
                     t,k,psurf(t),unit="Pa",direction="DN",&
                     vmin=5000.0,vmax=110000.0)
             enddo
          enddo
       endif

     ! If rain or snow fall are selected forcings:
       if(LDT_FORC_Rainf%selectOpt.eq.1 .or. LDT_FORC_Snowf%selectOpt.eq.1) then
          do k=1,LDT_FORC_Rainf%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_Rainf%varname(k)),&
                  pcpField,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Rainf in diagnoseForcingOutput')

             call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Rainf in diagnoseForcingOutput')

             if(LDT_FORC_Rainf%selectOpt.eq.1) then
                do t=1,LDT_rc%ntiles(n)
                   call diagnoseForcingVar(n, LDT_FORC_Rainf,&
                        t,k,pcp(t),unit="kg/m2s",direction="DN",&
                        vmin=0.0,vmax=0.02)
                   call diagnoseForcingVar(n, LDT_FORC_Rainf,&
                        t,k,pcp(t)*LDT_rc%ts,unit="kg/m2",direction="DN",&
                        vmin=0.0,vmax=0.0)
                enddo
             endif
           ! Snowfall only field:
             if(LDT_FORC_Snowf%selectOpt.eq.1) then
                call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_Snowf%varname(k)),&
                     snowfField,rc=status)
                call LDT_verify(status,&
                     'error in ESMF_StateGet:Snowf in diagnoseForcingOutput')

                call ESMF_FieldGet(snowfField,localDE=0, farrayPtr=snowf,rc=status)
                call LDT_verify(status,&
                     'error in ESMF_FieldGet:Snowf in diagnoseForcingOutput')

                do t=1,LDT_rc%ntiles(n)
                   call diagnoseForcingVar(n, LDT_FORC_Snowf,&
                        t,k,snowf(t),unit="kg/m2s",direction="DN",&
                        vmin=0.0,vmax=0.02)
                   call diagnoseForcingVar(n, LDT_FORC_Snowf,&
                        t,k,snowf(t)*LDT_rc%ts,unit="kg/m2",direction="DN",&
                        vmin=0.0,vmax=0.0)
                enddo

              ! Calculate total precipitation field (rainf+snowf):
                do t=1,LDT_rc%ntiles(n)
                   if( snowf(t) >= 0. ) then
                      pcp(t) = pcp(t) + snowf(t)
                   endif
                enddo
             endif

           ! Write out (diagnose) total precipitation field:
             if(LDT_FORC_TotalPrecip%selectOpt.eq.1) then

!               call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_TotalPrecip%varname(k)),&
!                    totprecField,rc=status)
!               call LDT_verify(status,&
!                    'error in ESMF_StateGet:TotalPrecip in diagnoseForcingOutput')

!               call ESMF_FieldGet(totprecField,localDE=0, farrayPtr=pcp,rc=status)
!               call LDT_verify(status,&
!                    'error in ESMF_FieldGet:TotalPrecip in diagnoseForcingOutput')

               do t=1,LDT_rc%ntiles(n)
                 call diagnoseForcingVar(n, LDT_FORC_TotalPrecip, &
                      t,k,pcp(t),unit="kg/m2s",direction="DN",&
                      vmin = 0.0, vmax=0.02)
                 call diagnoseForcingVar(n, LDT_FORC_TotalPrecip, &
                      t,k,pcp(t)*LDT_rc%ts,unit="kg/m2",direction="DN",&
                      vmin=0.0, vmax=0.0)
               enddo
             endif

          enddo
       endif

     ! Convective rain field:
       if(LDT_FORC_CRainf%selectOpt.eq.1) then
          do k=1,LDT_FORC_CRainf%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_CRainf%varname(k)),&
                  cpcpField,rc=status)
             call LDT_verify(status)

             call ESMF_FieldGet(cpcpField,localDE=0, farrayPtr=cpcp,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_CRainf,&
                     t,k,cpcp(t),unit="kg/m2s",direction="DN",&
                     vmin=0.0,vmax=0.02)
                call diagnoseForcingVar(n, LDT_FORC_CRainf,&
                     t,k,cpcp(t)*LDT_rc%ts,unit="kg/m2",direction="DN",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Wind_N%selectOpt.eq.1) then
          do k=1,LDT_FORC_Wind_N%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Wind_N%varname(k)),&
                  vField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Wind_N in diagnoseForcingOutput')    
             
             call ESMF_FieldGet(vField, localDE=0, farrayPtr=vwind,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Wind_N in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_WIND_N,&
                     t,k,vwind(t),unit="m/s",direction="N",&
                     vmin=-75.0,vmax=75.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Wind_E%selectOpt.eq.1) then
          do k=1,LDT_FORC_Wind_E%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Wind_E%varname(k)),&
                  uField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Wind_E in diagnoseForcingOutput')    
             
             call ESMF_FieldGet(uField, localDE=0, farrayPtr=uwind,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Wind_E in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_WIND_E,&
                     t,k,uwind(t),unit="m/s",direction="E",&
                     vmin=-75.0,vmax=75.0)
             enddo
          enddo
       endif

       if(LDT_FORC_SWdirect%selectOpt.eq.1) then
          do k=1,LDT_FORC_SWdirect%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_SWdirect%varname(k)),&
                  swdirField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:SWdirect in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(swdirField, localDE=0, farrayPtr=swdir,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:SWdirect in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_SWdirect,&
                     t,k,swdir(t),unit="W/m2",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_SWdiffuse%selectOpt.eq.1) then
          do k=1,LDT_FORC_SWdiffuse%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_SWdiffuse%varname(k)),&
                  swdifField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:SWdiffuse in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(swdifField, localDE=0, farrayPtr=swdif,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:SWdiffuse in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_SWdiffuse,&
                     t,k,swdif(t),unit="W/m2",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif


       if(LDT_FORC_Forc_Hgt%selectOpt.eq.1) then
          do k=1,LDT_FORC_Forc_Hgt%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Forc_Hgt%varname(k)),&
                  hField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Forc_Hgt in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(hField, localDE=0, farrayPtr=harray,rc=status)
             call LDT_verify(status,'error in ESMF_FieldGet:Forc_Hgt in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Forc_Hgt,&
                     t,k,harray(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Ch%selectOpt.eq.1) then
          do k=1,LDT_FORC_Ch%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Ch%varname(k)),&
                  chField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Ch in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(chField, localDE=0, farrayPtr=charray,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Ch in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Ch,&
                     t,k,charray(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Cm%selectOpt.eq.1) then
          do k=1,LDT_FORC_Cm%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Cm%varname(k)),&
                  cmField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Cm in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(cmField, localDE=0, farrayPtr=cmarray,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Cm in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Cm,&
                     t,k,cmarray(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Emiss%selectOpt.eq.1) then
          do k=1,LDT_FORC_Emiss%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Emiss%varname(k)),&
                  emissField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Emiss in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(emissField, localDE=0, farrayPtr=emiss,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Emiss in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_EMiss,&
                     t,k,emiss(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Q2sat%selectOpt.eq.1) then
          do k=1,LDT_FORC_Q2sat%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Q2sat%varname(k)),&
                  mixField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Q2sat in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(mixField, localDE=0, farrayPtr=mix,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Q2sat in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Q2sat,&
                     t,k,mix(t),unit="kg/kg",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Cosz%selectOpt.eq.1) then
          do k=1,LDT_FORC_Cosz%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Cosz%varname(k)),&
                  coszField, rc=status)
             call LDT_verify(status,&
                  'error ESMF_StateGet:Cosz in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(coszField, localDE=0, farrayPtr=cosz,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Cosz in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Cosz,&
                     t,k,cosz(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Alb%selectOpt.eq.1) then
          do k=1,LDT_FORC_Alb%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_Alb%varname(k)),&
                  albField, rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:Alb in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(albField, localDE=0, farrayPtr=alb,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:Alb in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_Alb,&
                     t,k,alb(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Pardr%selectOpt.eq.1) then
          do k=1,LDT_FORC_Pardr%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_Pardr%varname(k)),&
                  q2Field,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:PARDR in diagnoseForcingOutput')

             call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:PARDR in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_PARDR,&
                     t,k,q2(t),unit="W/m2",direction="DN",&
                     vmin=0.0,vmax=1360.0)
             enddo
          enddo
       endif

       if(LDT_FORC_Pardf%selectOpt.eq.1) then
          do k=1,LDT_FORC_Pardf%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_Pardf%varname(k)),&
                  q2Field,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_StateGet:PARDF in diagnoseForcingOutput')

             call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
             call LDT_verify(status,&
                  'error in ESMF_FieldGet:PARDF in diagnoseForcingOutput')

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_PARDF,&
                     t,k,q2(t),unit="W/m2",direction="DN",&
                     vmin=0.0,vmax=1360.0)
             enddo
          enddo
       endif

!<for vic>
       if(LDT_FORC_SNOWFLAG%selectOpt.eq.1) then
          do k=1,LDT_FORC_SNOWFLAG%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_SNOWFLAG%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_SNOWFLAG,&
                     t,k,tempPtr(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_DENSITY%selectOpt.eq.1) then
          do k=1,LDT_FORC_DENSITY%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_DENSITY%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_DENSITY,&
                     t,k,tempPtr(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_VAPORPRESS%selectOpt.eq.1) then
          do k=1,LDT_FORC_VAPORPRESS%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_VAPORPRESS%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_VAPORPRESS,&
                     t,k,tempPtr(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_VAPORPRESSDEFICIT%selectOpt.eq.1) then
          do k=1,LDT_FORC_VAPORPRESSDEFICIT%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_VAPORPRESSDEFICIT%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_VAPORPRESSDEFICIT,&
                     t,k,tempPtr(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

       if(LDT_FORC_WIND%selectOpt.eq.1) then
          do k=1,LDT_FORC_WIND%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),trim(LDT_FORC_WIND%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_WIND,&
                     t,k,tempPtr(t),unit="-",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif
!</for vic>
     ! Begin for FEWSNET
       if(LDT_FORC_PET%selectOpt.eq.1) then
          do k=1,LDT_FORC_PET%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_PET%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)

             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_PET,&
                     t,k,tempPtr(t),unit="kg/m2s",direction="-",&
                     vmin=0.0,vmax=0.0)
                call diagnoseForcingVar(n, LDT_FORC_PET,&
                     t,k,tempPtr(t)*LDT_rc%ts,unit="kg/m2",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif
       if(LDT_FORC_RefET%selectOpt.eq.1) then
          do k=1,LDT_FORC_RefET%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_RefET%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)

             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_RefET,&
                     t,k,tempPtr(t),unit="kg/m2",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif
     ! End for FEWSNET

       if(LDT_FORC_CAPE%selectOpt.eq.1) then
          do k=1,LDT_FORC_CAPE%vlevels
             call ESMF_StateGet(LDT_FORC_State(n),(LDT_FORC_CAPE%varname(k)),&
                  tempField, rc=status)
             call LDT_verify(status)

             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LDT_verify(status)

             do t=1,LDT_rc%ntiles(n)
                call diagnoseForcingVar(n, LDT_FORC_CAPE,&
                     t,k,tempPtr(t),unit="J/kg",direction="-",&
                     vmin=0.0,vmax=0.0)
             enddo
          enddo
       endif

    endif

  end subroutine diagnoseForcingOutput



!BOP
! 
! !ROUTINE: diagnoseForcingVar
! \label{diagnoseForcingVar}
! 
! !INTERFACE:
  subroutine diagnoseForcingVar(n, dataentry, t, vlevel,in_value, unit, &
                                direction, vmin, vmax)
! !USES: 

    implicit none
! !ARGUMENTS:     
    integer                 :: n 
    type(LDT_forcdataEntry) :: dataEntry
    integer                 :: t
    integer                 :: vlevel
    real, intent(in)        :: in_value
    character(len=*)        :: unit
    character(len=*)        :: direction
    real                    :: vmin
    real                    :: vmax
! 
! !DESCRIPTION: 
!  This routine maps a single output variable to the appropriate variable 
!  in the generic list of the LDT history writer. 
!EOP
    integer                 :: i
    logical                 :: unit_status
    logical                 :: dir_status
    real                    :: mfactor
    real                    :: value
    integer                 :: sftype
    character(len=20)       :: cfunit

    sftype = LDT_domain(n)%tile(t)%sftype

    unit_status = .false.
!    do i=1,dataEntry%nunits
    if(unit.eq.dataEntry%units) then 
       unit_status = .true. 
    endif
!    enddo
!    
!    dir_status = .false. 
!    do i=1,dataEntry%ndirs
!       if(direction.eq.dataEntry%dirtypes(i)) then 
!          dir_status = .true. 
!          exit
!       endif
!    enddo
!   
    if(unit_status) then 
       call convertToCFunits(unit,cfunit)
!       if(cfunit.eq.dataEntry%units) then 
          ! it is assumed that there will be only two 
          ! directions. 
          ! Determine "multiplication" factor ("mfactor")
!          if(direction.eq.dataEntry%dir) then 
       mfactor = 1.0
!          else
!             mfactor = -1.0
!          endif
          
          ! Correct the direction of value
       value = in_value * mfactor
       
       if(mfactor.eq.1) then 
          dataEntry%valid_min = vmin
          dataEntry%valid_max = vmax
       else
          dataEntry%valid_min = vmax
          dataEntry%valid_max = vmin
       endif
       
       ! accumulate values and record instantaneous values
       if( dataEntry%timeAvgOpt.eq.2 ) then 
          dataEntry%modelOutput(1,t,vlevel) = &   ! Accum.
               dataEntry%modelOutput(1,t,vlevel) + value
          dataEntry%modelOutput(2,t,vlevel) = value
          dataEntry%count(sftype,vlevel) = &      ! Inst.
               dataEntry%count(sftype,vlevel)+1

       ! accumulate values
       elseif(dataEntry%timeAvgOpt.eq.1 .or. &    ! Ave.
            dataEntry%timeAvgOpt.eq.3) then       ! Sum

          dataEntry%modelOutput(1,t,vlevel) = &
               dataEntry%modelOutput(1,t,vlevel) + value

          dataEntry%count(sftype,vlevel) = &
               dataEntry%count(sftype,vlevel)+1

       ! record instantaneous values
       else  ! value = 0
          dataEntry%modelOutput(1,t,vlevel) = value
          dataEntry%count(sftype,vlevel) = 1
       endif

    endif
       !    endif

       !    if(.not.unit_status) then 
       !       write(LDT_logunit,*) 'Error: ',trim(dataEntry%units),&
       !            ' for field ',trim(dataEntry%standard_name),' is not defined '
       !       write(LDT_logunit,*) 'for diagnostic output...'
       !       write(LDT_logunit,*) 'supported unit types: ',dataEntry%unittypes
       !       write(LDT_logunit,*) 'Program stopping ..'
       !       call LDT_endrun()       
       !    endif
       !    if(.not.dir_status) then 
       !       write(LDT_logunit,*) 'Error: ',trim(dataEntry%dir),&
       !            ' for field ',trim(dataEntry%standard_name),' is not defined '
       !       write(LDT_logunit,*) 'for diagnostic output...'
       !!       write(LDT_logunit,*) 'supported direction types: ',&
       !!            dataEntry%dir
       !       write(LDT_logunit,*) 'Program stopping ..'
       !       call LDT_endrun()       
       !    endif
  end subroutine diagnoseForcingVar

!BOP
! 
! !ROUTINE: LDT_output_met_forcing
! \label{LDT_output_met_forcing}
! 
! !INTERFACE:
  subroutine LDT_output_met_forcing(n)

! !USES: 
    use LDT_fileIOMod

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine sets up the final output writing 
!   calls for the meteorological forcing variables.
!
!EOP

    logical         :: alarmCheck
    character(len=LDT_CONST_PATH_LEN)   :: outfile
    character*120   :: source_names
    integer         :: m, i
    
    alarmCheck = .false. 
    if(LDT_rc%wopt.ne."none") then 
       alarmCheck = LDT_isAlarmRinging(LDT_rc,&
            "LDT metforcing output alarm")
    endif
    
    if(alarmCheck) then 

      if(LDT_masterproc) then
         call LDT_create_output_directory('FORCING')
         call LDT_create_output_filename(n,outfile,&
              model_name = 'FORCING',&
              writeint=LDT_rc%metForcOutInterval)
      endif

    ! Write out full list of selected metforcing names:
      i = 1
      source_names = ""
      do m = 1,LDT_rc%nmetforc
         if( m /= LDT_rc%nmetforc ) then
           source_names(i:(i+LEN_TRIM(LDT_rc%metforc(m))+1)) = &
                  TRIM(LDT_rc%metforc(m))//","
         else
           source_names(i:(i+LEN_TRIM(LDT_rc%metforc(m))+1)) = &
                  TRIM(LDT_rc%metforc(m))
         endif
         i = i + LEN_TRIM(LDT_rc%metforc(m))+1
      enddo
      call writeForcingOutput(n,outfile,&
           outInterval=LDT_rc%metForcOutInterval,&
           model_name = source_names )
!           model_name = 'FORCING')
    endif

  end subroutine LDT_output_met_forcing

!BOP
! !ROUTINE: writeForcingOutput
! \label{writeForcingOutput}
! 
! !INTERFACE: 
  subroutine writeForcingOutput(n,outfile,outInterval,model_name)
! !USES:

! !ARGUMENTS: 
    integer,   intent(in)   :: n 
    real,      intent(in)   :: outInterval
    character(len=*), intent(in) :: outfile
    character(len=*), intent(in),optional :: model_name
! 
! !DESCRIPTION: 
!   This subroutine invokes the routine to write forcing output in the selected
!    data format (binary/grib1/netcdf) and using the selected list of variables. 
!    Further, the variables are also written as instantaneous, time averaged,
!    or accumulated based on the user specifications. 
!   
!   The arguments are: 
!   \begin{description}
!    \item[n]  index of the nest \newline
!    \item[lsmoutfile]  name of the LSM history file \newline
!  \end{description}
! 
!   The routines invoked are: 
!   \begin{description}
!   \item[writeBinaryOutput](\ref{writeBinaryOutput}) \newline
!     writes the history files in binary format
!   \item[writeGrib1Output](\ref{writeGrib1Output}) \newline
!     writes the history files in GRIB1 format
!   \item[writeNetcdfOutput](\ref{writeNetcdfOutput}) \newline
!     writes the history files in NETCDF format
!   \item[LDT\_resetOutputVars](\ref{LDT_resetOutputVars}) \newline
!     resets the time averaged varibles for the next output. 
!   \end{description}
!EOP

    integer :: ftn
    integer :: iret
    character(len=LDT_CONST_PATH_LEN) :: mname_temp
    
    if(.NOT.PRESENT(model_name)) then 
       mname_temp = "model_not_specified"
    else
       mname_temp = model_name
    endif

    call rescaleForcingCount(n)

    if(LDT_rc%wout.eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       if(LDT_masterproc) then 
#if (defined USE_NETCDF4)
          iret = nf90_create(path=outfile,cmode=nf90_netcdf4,&
               ncid = ftn)
          call LDT_verify(iret,'creating netcdf file failed')
#endif
#if (defined USE_NETCDF3)
          iret = nf90_create(path=outfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LDT_verify(iret,'creating netcdf file failed')
#endif
       endif
       
     ! Write out meteorological forcing output:
       write(LDT_logunit,*)"... Writing output file: ",trim(outfile) 
       call writeNetcdfOutput(n,ftn,outInterval,mname_temp)
       
       if(LDT_masterproc) then 
          iret = nf90_close(ftn)
       endif
#endif
    endif
    ! After writing reset the variables
    call LDT_resetForcingVars()

  end subroutine writeForcingOutput

!BOP
! 
! !ROUTINE: convertToCFunits
! \label{convertToCFunits}
!
! !INTERFACE: 
  subroutine convertToCFunits(unit,cfunit)
!
! !DESCRIPTION: 
!   This routine converts the LIS specified units to CF compliant
!   unit specifications. 
! 
!EOP

    implicit none

    character(len=*)  :: unit
    character(len=*), intent(out)  :: cfunit
    cfunit = ""
    if(unit.eq."W/m2") then 
       cfunit = "W m-2"
    elseif(unit.eq."J/m2") then 
       cfunit = "J m-2"
    elseif(unit.eq."kg/m2s") then 
       cfunit = "kg m-2 s-1"
    elseif(unit.eq."kg/m2") then 
       cfunit = "kg m-2"
    elseif(unit.eq."kg/m2s2") then 
       cfunit = "kg m-2 s-2"
    elseif(unit.eq."m3/m3") then 
       cfunit = "m^3 m-3"
    elseif(unit.eq."m/s") then 
       cfunit = "m s-1"
    elseif(unit.eq."N/m2") then 
       cfunit = "N m-2"
    elseif(unit.eq."kg/kg") then 
       cfunit = "kg kg-1"
    elseif(unit.eq."g/g") then 
       cfunit = "g g-1"
    elseif(unit.eq."umol/m2s") then 
       cfunit = "umol m-2 s-1"
    elseif(unit.eq."mm/hr") then 
       cfunit = "mm hr-1"
    elseif(unit.eq."-") then
       cfunit="1"
    else
       cfunit = unit
    endif
  end subroutine convertToCFunits


!BOP
! !ROUTINE: rescaleForcingCount
! \label{rescaleForcingCount}
! 
! !INTERFACE: 
  subroutine rescaleForcingCount(n)
! 
    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
     
! !DESCRIPTION: 
!EOP

    call rescaleSingleForcingCount(n,LDT_FORC_Tair)
    call rescaleSingleForcingCount(n,LDT_FORC_Qair)
    call rescaleSingleForcingCount(n,LDT_FORC_SWdown)
    call rescaleSingleForcingCount(n,LDT_FORC_SWdirect)
    call rescaleSingleForcingCount(n,LDT_FORC_SWdiffuse)
    call rescaleSingleForcingCount(n,LDT_FORC_LWdown)
    call rescaleSingleForcingCount(n,LDT_FORC_Wind_E)
    call rescaleSingleForcingCount(n,LDT_FORC_Wind_N)
    call rescaleSingleForcingCount(n,LDT_FORC_Psurf)
    call rescaleSingleForcingCount(n,LDT_FORC_Rainf)
    call rescaleSingleForcingCount(n,LDT_FORC_Snowf)
    call rescaleSingleForcingCount(n,LDT_FORC_CRainf)
    call rescaleSingleForcingCount(n,LDT_FORC_TotalPrecip)

  end subroutine rescaleForcingCount

  subroutine rescaleSingleForcingCount(n,dataEntry)

    use LDT_mpiMod

    integer                  :: n 
    type(LDT_forcDataEntry)  :: dataEntry

    integer :: count
    integer :: k, m
    integer :: ierr
    
    if(dataEntry%selectProc.ne.0) then 
       do k=1,dataEntry%vlevels
#if (defined SPMD)
          call mpi_reduce(sum(dataEntry%count(:,k)),&
               count, 1, MPI_INTEGER, MPI_SUM, 0, &
               MPI_COMM_WORLD,ierr)
#else
          count = sum(dataEntry%count(:,k))
#endif
          if(LDT_masterproc) then 
             if(count.eq.0) then 
                write(LDT_logunit,*) 'Error:',dataEntry%standard_name,&
                     ' field is not defined'
                write(LDT_logunit,*) 'for diagnostic output...'
                write(LDT_logunit,*) 'Please exclude it from the model output attributes table'
                write(LDT_logunit,*) 'Program stopping ..'
                call LDT_endrun()
             endif
          endif
       enddo
       
       ! Rescale the count when the timeAvgOpt indicates to use
       ! time averaging.
       ! If timeAvgOpt indicates instantaneous-only or accumulate
       ! then there is no need to rescale the count
       if ( (dataEntry%timeAvgOpt == 1) .or. &
            (dataEntry%timeAvgOpt == 2) ) then 
          do k=1,dataEntry%vlevels
             do m=1,LDT_rc%max_model_types
                if(dataEntry%count(m,k).gt.0) then 
                   dataEntry%count(m,k) = dataEntry%count(m,k) / &
                        LDT_rc%npatch(n,m)
                endif
             enddo
          enddo
       endif
    endif
  end subroutine rescaleSingleForcingCount

!BOP
! !ROUTINE: LDT_resetForcingVars
! \label{LDT_resetForcingVars}
! 
! !INTERFACE: 
  subroutine LDT_resetForcingVars()
! !ARGUMENTS:

! 
! !DESCRIPTION: 
!   This routine resets the specified variables to enable time averaging 
!    for the next history output step. 
!
!   It also resets the minimum and maximum fields.
!   
!EOP

    call resetSingleForcingVar(LDT_FORC_Tair)
    call resetSingleForcingVar(LDT_FORC_Qair)
    call resetSingleForcingVar(LDT_FORC_SWdown)
    call resetSingleForcingVar(LDT_FORC_SWdirect)
    call resetSingleForcingVar(LDT_FORC_SWdiffuse)
    call resetSingleForcingVar(LDT_FORC_LWdown)
    call resetSingleForcingVar(LDT_FORC_Wind_E)
    call resetSingleForcingVar(LDT_FORC_Wind_N)
    call resetSingleForcingVar(LDT_FORC_Psurf)
    call resetSingleForcingVar(LDT_FORC_Rainf)
    call resetSingleForcingVar(LDT_FORC_Snowf)
    call resetSingleForcingVar(LDT_FORC_CRainf)
    call resetSingleForcingVar(LDT_FORC_TotalPrecip)

  end subroutine LDT_resetForcingVars


!BOP
! !ROUTINE: resetSingleForcingVar
! \label{resetSingleForcingVar}
! 
! !INTERFACE: 
  subroutine resetSingleForcingVar(dataEntry)

  implicit none

! !ARGUMENTS:
    type(LDT_forcDataEntry)  :: dataEntry

! !DESCRIPTION: 
!   This routine resets the specified variable to enable time averaging 
!   for the next history output step. 
!
!   It also resets the minimum and maximum fields.
!
! The arguments are:
!   \begin{description}
!   \item[dataEntry] data structure to reset
!   \end{description}
!EOP
      dataEntry%modelOutput = 0.0
      dataEntry%count = 0 

  end subroutine resetSingleForcingVar


!BOP
! !ROUTINE: writeNetcdfOutput
! \label{writeNetcdfOutput}
! 
! !INTERFACE: writeNetcdfOutput
  subroutine writeNetcdfOutput(n, ftn, outInterval, model_name)
! !USES: 

! !ARGUMENTS: 
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: ftn
    real,      intent(in)   :: outInterval
    character*100,intent(in) :: model_name
! 
! !DESCRIPTION: 
!  This routine writes an output file in the NETCDF format based on the 
!  list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[outInterval]   history output frequency
!    \item[nsoillayers]  Number of soil layers
!    \item[lyrthk]   Thickness of soil layers
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{defineNETCDFheaderVar})
!     writes the required headers for a single variable
!   \item[writeSingleNETCDFvar](\ref{writeSingleNETCDFvar})
!     writes a single variable into a netcdf formatted file. 
!   \item[LDT\_verify](\ref{LDT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    integer               :: dimID(3)
    integer               :: tdimID,xtimeId
    integer               :: t,c,r,index1,i,m
    type(LDT_forcdataEntry), pointer :: xlat, xlong
    character*8           :: xtime_begin_date
    character*6           :: xtime_begin_time
    character*50          :: xtime_units
    character*50          :: xtime_timeInc
    integer               :: iret
! Note that the fix to add lat/lon to the NETCDF output will output
! undefined values for the water points. 
    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    type(LDT_forcdataEntry), pointer :: dataEntry

    character(100)        :: zterp_flag
! __________________________________________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
    call date_and_time(date,time,zone,values)
    
    allocate(xlat)
    allocate(xlong)

    xlat%short_name = "lat"
    xlat%standard_name = "latitude"
    xlat%units = "degree_north"
    xlat%nunits = 1
    xlat%vlevels = 1
    xlat%timeAvgOpt = 0 
    xlat%selectOpt = 1
    allocate(xlat%modelOutput(1,LDT_rc%ntiles(n),xlat%vlevels))
    allocate(xlat%count(1,xlat%vlevels))
    xlat%count = 1
    allocate(xlat%unittypes(1))
    xlat%unittypes(1) = "degree_north"
    xlat%valid_min = 0.0
    xlat%valid_max = 0.0

    xlong%short_name = "lon"
    xlong%standard_name = "longitude"
    xlong%units = "degree_east"
    xlong%nunits = 1
    xlong%vlevels = 1
    xlong%timeAvgOpt = 0 
    xlong%selectOpt = 1
    allocate(xlong%modelOutput(1,LDT_rc%ntiles(n),xlong%vlevels))
    allocate(xlong%count(1,xlong%vlevels))
    xlong%count = 1
    allocate(xlong%unittypes(1))
    xlong%unittypes(1) = "degree_east"
    xlong%valid_min = 0.0
    xlong%valid_max = 0.0

    if(LDT_masterproc) then 
       if(LDT_rc%wopt.eq."1d tilespace") then 
          call LDT_verify(nf90_def_dim(ftn,'ntiles',LDT_rc%glbntiles(n),&
               dimID(1)),&
               'nf90_def_dim for ntiles failed in LDT_metforcingMod')
       elseif(LDT_rc%wopt.eq."2d gridspace") then 
          call LDT_verify(nf90_def_dim(ftn,'east_west',LDT_rc%gnc(n),&
               dimID(1)),&
               'nf90_def_dim for east_west failed in LDT_metforcingMod')
          call LDT_verify(nf90_def_dim(ftn,'north_south',LDT_rc%gnr(n),&
               dimID(2)),&
               'nf90_def_dim for north_south failed in LDT_metforcingMod')
       endif

!LDT output is always writing output for a single time
       call LDT_verify(nf90_def_dim(ftn,'time',1,tdimID),&
            'nf90_def_dim for time failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", LDT_rc%udef),&
            'nf90_put_att for missing_value failed in LDT_metforcingMod')
       
       call defineNETCDFheaderVar(n,ftn,dimID, xlat,&
            non_model_fields = .true. )       
       call defineNETCDFheaderVar(n,ftn,dimID, xlong, &
            non_model_fields = .true. )              
!defining time field

       call LDT_verify(nf90_def_var(ftn,'time',&
            nf90_float,dimids = tdimID, varID=xtimeID),&
            'nf90_def_var for time failed in LDT_metforcingMod')
       
       write(xtime_units,200) LDT_rc%yr, LDT_rc%mo, LDT_rc%da, &
            LDT_rc%hr, LDT_rc%mn, LDT_rc%ss
200    format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
            I2.2,':',I2.2)
       write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
            LDT_rc%yr, LDT_rc%mo, LDT_rc%da
       write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
            LDT_rc%hr, LDT_rc%mn, LDT_rc%ss
       write(xtime_timeInc, fmt='(I20)') &
            nint(outInterval)
       
       call LDT_verify(nf90_put_att(ftn,xtimeID,&
            "units",trim(xtime_units)),&
            'nf90_put_att for units failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,xtimeID,&
            "long_name","time"),&
            'nf90_put_att for long_name failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,xtimeID,&
            "time_increment",trim(adjustl(xtime_timeInc))),&
            'nf90_put_att for time_increment failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,xtimeID,&
            "begin_date",xtime_begin_date),&
            'nf90_put_att for begin_date failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,xtimeID,&
            "begin_time",xtime_begin_time),&
            'nf90_put_att for begin_time failed in LDT_metforcingMod')

       
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_Tair)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_Qair)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_SWdown)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_SWdirect)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_SWdiffuse)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_LWdown)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_Wind_E)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_Wind_N)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_Psurf)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_Rainf)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_Snowf)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_CRainf)
       call defineNETCDFheaderVar(n,ftn,dimId,LDT_FORC_TotalPrecip)

       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
            "LDT metforcing output"),&
            'nf90_put_att for title failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",&
            trim(model_name)),&
            'nf90_put_att for source failed in LDT_metforcingMod')

     ! Write out full string with Zenith interp flags listed:
       zterp_flag = "";  i = 1
       do m = 1,LDT_rc%nmetforc
          if( LDT_rc%met_zterp(m) .eqv. .false. ) then
            zterp_flag(i:(i+6)) = "false,"; i = i + 6
          else
            zterp_flag(i:(i+5)) = "true,";  i = i + 5
          endif
       enddo
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"zenith_interp",&
            zterp_flag),&
            'nf90_put_att for source failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
            trim(LDT_rc%institution)),&
            'nf90_put_att for institution failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
            date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
            'nf90_put_att for history failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
            "Arsenault_etal_GMD_2018, Kumar_etal_EMS_2006"),&
            'nf90_put_att for references failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
            "CF-1.6"),'nf90_put_att for conventions failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
            "website: http://lis.gsfc.nasa.gov/"),&
            'nf90_put_att for comment failed in LDT_metforcingMod')

     ! -- Grid information --
       select case ( LDT_rc%lis_map_proj(n) )

        case( "latlon" )
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "EQUIDISTANT CYLINDRICAL"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LAT", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LON", &
               LDT_rc%gridDesc(n,8)))

          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,9)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,10)))       

        case( "mercator" )
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "MERCATOR"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

        case( "lambert" )   ! Lambert conformal
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "LAMBERT CONFORMAL"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

        case( "polar" )    ! polar stereographic
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "POLAR STEREOGRAPHIC"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

        case default 
       end select
       call LDT_verify(nf90_enddef(ftn))
    endif

    do t=1,LDT_rc%ntiles(n)
       c = LDT_domain(n)%tile(t)%col
       r = LDT_domain(n)%tile(t)%row
       index1 = LDT_domain(n)%gindex(c,r)
       xlat%modelOutput(1,t,1) = LDT_domain(n)%grid(index1)%lat
       xlong%modelOutput(1,t,1) = LDT_domain(n)%grid(index1)%lon
    enddo
    call writeSingleNETCDFvar(ftn,n,xlat, non_model_fields = .true.)
    call writeSingleNETCDFvar(ftn,n,xlong,non_model_fields = .true.)

    if(LDT_masterproc) then 
       call LDT_verify(nf90_put_var(ftn,xtimeID,0.0),&
            'nf90_put_var for xtimeId failed in LDT_metforcingMod')
    endif
    
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_Tair)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_Qair)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_SWdown)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_SWdirect)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_SWdiffuse)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_LWdown)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_Wind_E)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_Wind_N)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_Psurf)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_Rainf)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_Snowf)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_CRainf)
    call writeSingleNETCDFvar(ftn,n,LDT_FORC_TotalPrecip)

    deallocate(xlat%modelOutput)
    deallocate(xlat%count)
    deallocate(xlat%unittypes)
    deallocate(xlat)

    deallocate(xlong%modelOutput)
    deallocate(xlong%count)
    deallocate(xlong%unittypes)
    deallocate(xlong)
#endif
  end subroutine writeNetcdfOutput

!BOP
! !ROUTINE: defineNETCDFheaderVar
! \label{defineNETCDFheaderVar}
! 
! !INTERFACE: 
  subroutine defineNETCDFheaderVar(n,ftn,dimID, dataEntry, non_model_fields)
! !USES: 

! !ARGUMENTS:     
    integer               :: n
    integer               :: ftn
    type(LDT_forcdataEntry) :: dataEntry
    logical,   optional   :: non_model_fields
    integer               :: dimID(3)
! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single variable
! 
!   The arguments are: 
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LDT\_endrun](\ref{LDT_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LDT\_verify](\ref{LDT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    logical       :: nmodel_status
    integer       :: shuffle, deflate, deflate_level
    character*100 :: short_name
    integer       :: fill_value
    logical       :: write_flag

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    nmodel_status = .false.
    if(present(non_model_fields)) then 
       nmodel_status = non_model_fields
    endif

    write_flag = .true. 

    if(nmodel_status) then 
       write_flag = .true.     
    elseif(dataEntry%selectProc.ne.1) then 
       write_flag = .false. 
    endif

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(write_flag) then 
       if(LDT_rc%wopt.eq."1d tilespace") then 
          if(dataEntry%vlevels.gt.1) then 
             call LDT_verify(nf90_def_dim(ftn,&
                  trim(dataEntry%short_name)//'_profiles',&
                  dataEntry%vlevels, dimID(2)),&
                  'nf90_def_dim failed (2d gridspace) in LDT_metforcingMod')
          endif
       elseif(LDT_rc%wopt.eq."2d gridspace") then 
          if(dataEntry%vlevels.gt.1) then 
             call LDT_verify(nf90_def_dim(ftn,&
                  trim(dataEntry%short_name)//'_profiles',&
                  dataEntry%vlevels, dimID(3)),&
                  'nf90_def_dim failed (2d gridspace) in LDT_metforcingMod')
          endif
       endif
       if(LDT_rc%wopt.eq."1d tilespace") then              
          if(dataEntry%timeAvgOpt.eq.2) then 
             if(dataEntry%vlevels.gt.1) then 
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')

                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
                
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#endif
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_opt),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
                
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
#endif
             else
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:1), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)                
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID(1:1), varID=dataEntry%varId_opt),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)                
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
             endif
          else
             if(.not.nmodel_status) then 
                if(dataEntry%timeAvgOpt.eq.0) then
                   short_name = trim(dataEntry%short_name)//'_inst'
                elseif(dataEntry%timeAvgOpt.eq.1) then
                   short_name = trim(dataEntry%short_name)//'_tavg'
                elseif(dataEntry%timeAvgOpt.eq.3) then
                   short_name = trim(dataEntry%short_name)//'_acc'
                endif
             else
                short_name = trim(dataEntry%short_name)
             endif
             
             if(dataEntry%vlevels.gt.1) then 
                call LDT_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
                
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
#endif
             else
                call LDT_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID(1:1), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)                
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
             endif
          endif
       elseif(LDT_rc%wopt.eq."2d gridspace") then 
          if(dataEntry%timeAvgOpt.eq.2) then 
             if(dataEntry%vlevels.gt.1) then 
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID, varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID, varID=dataEntry%varId_opt),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
                
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
             else
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_inst',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
                call LDT_verify(nf90_def_var(ftn,trim(dataEntry%short_name)//'_tavg',&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_opt),&
                     'nf90_def_var for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_opt, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_opt,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
             endif
          else
             if(.not.nmodel_status) then 
                if(dataEntry%timeAvgOpt.eq.0) then
                   short_name = trim(dataEntry%short_name)//'_inst'
                elseif(dataEntry%timeAvgOpt.eq.1) then
                   short_name = trim(dataEntry%short_name)//'_tavg'
                elseif(dataEntry%timeAvgOpt.eq.3) then
                   short_name = trim(dataEntry%short_name)//'_acc'
                endif
             else
                short_name = trim(dataEntry%short_name)
             endif
             
             if(dataEntry%vlevels.gt.1) then 
                call LDT_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID, varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif
             else
                call LDT_verify(nf90_def_var(ftn,trim(short_name),&
                     nf90_float,&
                     dimids = dimID(1:2), varID=dataEntry%varId_def),&
                     'nf90_def_var for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
                call LDT_verify(nf90_def_var_fill(ftn,&
                     dataEntry%varId_def, &
                     1,fill_value), 'nf90_def_var_fill failed for '//&
                     dataEntry%short_name)                
#if(defined USE_NETCDF4)
                call LDT_verify(nf90_def_var_deflate(ftn,&
                     dataEntry%varId_def,&
                     shuffle, deflate, deflate_level),&
                     'nf90_def_var_deflate for '//trim(short_name)//&
                     'failed in defineNETCDFheadervar')                     
#endif                
             endif
          endif
       endif
       
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "units",trim(dataEntry%units)),&
            'nf90_put_att for units failed in defineNETCDFheaderVar')
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "standard_name",trim(dataEntry%standard_name)),&
            'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "scale_factor",1.0),&
            'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "add_offset",0.0),&
            'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "missing_value",LDT_rc%udef),&
            'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "_FillValue",LDT_rc%udef),&
            'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "vmin",dataEntry%valid_min),&
            'nf90_put_att for vmin failed in defineNETCDFheaderVar')
       call LDT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "vmax",dataEntry%valid_max),&
            'nf90_put_att for vmax failed in defineNETCDFheaderVar')
       
       if(dataEntry%timeAvgOpt.eq.2) then
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "units",trim(dataEntry%units)),&
               'nf90_put_att for units failed in defineNETCDFheaderVar')
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "standard_name",trim(dataEntry%standard_name)),&
               'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "scale_factor",1.0),&
               'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "add_offset",0.0),&
               'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "missing_value",LDT_rc%udef),&
               'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "_FillValue",LDT_rc%udef),&
               'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "vmin",dataEntry%valid_min),&
               'nf90_put_att for vmin failed in defineNETCDFheaderVar')
          call LDT_verify(nf90_put_att(ftn,dataEntry%varId_opt,&
               "vmax",dataEntry%valid_max),&
               'nf90_put_att for vmax failed in defineNETCDFheaderVar')
          
       endif       
    endif
#endif
  end subroutine defineNETCDFheaderVar

!BOP
! !ROUTINE: writeSingleNETCDFvar
! \label{writeSingleNETCDFvar}
!
! !INTERFACE: 
  subroutine writeSingleNETCDFvar(ftn, n,dataEntry,&
       non_model_fields)
! !USES: 
    use LDT_coreMod,    only : LDT_rc
    use LDT_historyMod, only : LDT_writevar_netcdf

    implicit none

    integer,   intent(in)   :: n 
    integer,   intent(in)   :: ftn
    type(LDT_forcdataEntry) :: dataEntry
    logical,   optional     :: non_model_fields
! 
! !DESCRIPTION: 
!  This routine writes a single variable to a NETCDF file
!  The arguments are: 
!  \begin{description}
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[n] index of the nest
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LDT\_writevar\_netcdf](\ref{LDT_writevar_netcdf})
!     writes a variable into a netcdf formatted file. 
!   \end{description}
!EOP    
    integer       :: i,k,m,t
    logical       :: nmodel_status
    logical       :: write_flag

    nmodel_status = .false.
    if(present(non_model_fields)) then 
       nmodel_status = non_model_fields
    endif

    write_flag = .true. 

    if(nmodel_status) then 
       write_flag = .true.     
    elseif(dataEntry%selectProc.ne.1) then 
       write_flag = .false. 
    endif

    if(write_flag) then 
          
       if(.not.nmodel_status) then

          do t=1,LDT_rc%ntiles(n)
             if(nmodel_status) then 
                m = 1
             else
                m = LDT_domain(n)%tile(t)%sftype
             endif
             
             do k=1,dataEntry%vlevels
                if(dataEntry%count(m,k).gt.0) then 
                   if(dataEntry%timeAvgOpt.eq.3) then  !do nothing
                      continue   
                   elseif(dataEntry%timeAvgOpt.eq.2.or.dataEntry%timeAvgOpt.eq.1) then 
                      dataEntry%modelOutput(1,t,k) = dataEntry%modelOutput(1,t,k)/&
                           dataEntry%count(m,k)
                   else !do nothing
                      continue   
                   endif
                else
                   dataEntry%modelOutput(1,t,k) = LDT_rc%udef
                endif
             enddo
          enddo
       endif

       do k=1,dataEntry%vlevels

          ! accumulated values
          ! time-averaged values and instantaneous values
          if(dataEntry%timeAvgOpt.eq.2) then 
             call LDT_writevar_netcdf(ftn,n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%varId_def, &
                  dim1=k)
             call LDT_writevar_netcdf(ftn, n,&
                  dataEntry%modelOutput(2,:,k),&
                  dataEntry%varId_opt, &
                  dim1=k)
          ! time-averaged, summed or instantaneous values
          else
             call LDT_writevar_netcdf(ftn, n,&
                  dataEntry%modelOutput(1,:,k),&
                  dataEntry%varId_def, &
                  dim1=k)
          endif
       enddo
    endif

  end subroutine writeSingleNETCDFvar

end module LDT_metforcingMod

