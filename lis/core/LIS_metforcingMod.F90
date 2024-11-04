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
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_metforcingMod
!BOP
!
!  !MODULE: LIS_metforcingMod
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
!                                 from LIS_MOC_RAINFCONV to LIS_MOC_CRAINFFORC.
!                              Added units of [kg/m^2] for PET and CRainf.
! 
  use ESMF
  use LIS_FORC_AttributesMod
  use LIS_spatialDownscalingMod
  use LIS_logMod
  use LIS_coreMod

  implicit none 

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------

  public :: LIS_metforcing_init      ! initialize met forcing setup
  public :: LIS_get_met_forcing      ! retrieve data and interpolate 
                                     !  spatially and temporally
  public :: LIS_perturb_forcing      ! perturbs the met forcing variables
  public :: LIS_metforcing_reset     ! resets required forcing variables
  public :: LIS_metforcing_finalize  ! cleanup allocated structures
!  public :: diagnoseForcingOutput    ! prepares forcing variables for output

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public   :: LIS_forc      
  public   :: LIS_FORC_State
  public   :: LIS_FORC_Base_State
  public   :: LIS_FORC_Pert_State

!EOP

  type, public :: forc_dec_type
     real, allocatable :: modelelev(:)
  end type forc_dec_type

!BOP
! !ROUTINE: LIS_metforcing_init
! \label{LIS_metforcing_init}
!
! !INTERFACE: 
  interface LIS_metforcing_init
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure metforcing_init_offline
     module procedure metforcing_init_coupled
! !DESCRIPTION: 
!  Initialization routine for metforcing structures. There are
!  two private methods metd on whether LIS is run in an offline
!  mode or coupled mode. 
!EOP
  end interface

  type(forc_dec_type), allocatable :: LIS_forc(:,:)
  type(ESMF_State),    allocatable :: LIS_FORC_State(:)
  type(ESMF_State),    allocatable :: LIS_FORC_Base_State(:,:)
  type(ESMF_State),    allocatable :: LIS_FORC_Pert_State(:)

contains
!BOP
! !ROUTINE: metforcing_init_offline
! \label{metforcing_init_offline}
!
! !INTERFACE: 
! !Private name: call using LIS_metforcing_init
  subroutine metforcing_init_offline
! !USES:
    use LIS_metforcing_pluginMod, only : LIS_metforcing_plugin
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
!  \item[LIS\_metforcing\_plugin](\ref{LIS_metforcing_plugin}) \newline
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
    integer :: n,m,j
    integer :: nensem

    TRACE_ENTER("metf_init")
    if(LIS_rc%nmetforc.gt.0) then
       
       allocate(LIS_forc(LIS_rc%nnest,LIS_rc%nmetforc))

       call LIS_metforcing_plugin
       call LIS_init_pcpclimo()

       do n=1,LIS_rc%nnest      
          do m=1,LIS_rc%nmetforc
             if(LIS_rc%met_ecor(m).ne."none") then 
                write(LIS_logunit,*) "[INFO] Elevation correction turned on for:  ",&
                trim(LIS_rc%metforc(m))
                allocate(LIS_forc(n,m)%modelelev(LIS_rc%ngrid(n)))
             endif
          enddo
       enddo
          
     ! Read Forcing Table Attributes:
       call read_forctable_attributes( )

     ! Initialize Metforcing Readers and Parameters:
       do m=1,LIS_rc%nmetforc      
          call initmetforc(trim(LIS_rc%metforc(m))//char(0),m)  
       enddo

       nensem = 0 
       do j=1,LIS_rc%nmetforc
          nensem = nensem + LIS_rc%met_nensem(j)
       enddo
          
       do n=1,LIS_rc%nnest
          if(LIS_rc%metforc_blend_alg.eq."ensemble") then 
             if(LIS_rc%nensem(n).lt.LIS_rc%nmetforc) then 
                write(LIS_logunit,*) '[ERR] The number of ensembles can not be smaller than the'
                write(LIS_logunit,*) '[ERR] number of forcings in the ensemble of forcings mode.'
                write(LIS_logunit,*) '[ERR] Program stopping.... '
                call LIS_endrun()
             endif
             if(mod(LIS_rc%nensem(n), nensem).ne.0) then 
                write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple of the '
                write(LIS_logunit,*) '[ERR] number of ensembles from the forcings.'
                write(LIS_logunit,*) '[ERR] Program stopping.... '
                call LIS_endrun()
             else
                !equally divide the ensembles..
!                LIS_rc%nperforc = LIS_rc%nensem(n)/LIS_rc%nmetforc
                do j=1,LIS_rc%nmetforc
                   LIS_rc%met_nperforc(j) = (LIS_rc%nensem(n)/nensem)*LIS_rc%met_nensem(j)
                enddo
             endif
          endif
       enddo
       
       call create_forcing_structures()
       call forcingPerturbSetup()

    endif
    TRACE_EXIT("metf_init")

  end subroutine metforcing_init_offline

!BOP
! !ROUTINE: metforcing_init_coupled
! \label{metforcing_init_coupled}
!
! !INTERFACE: 
  subroutine metforcing_init_coupled(flag)
! !USES:

    integer, intent(in) :: flag
!
! !DESCRIPTION:
!
! This routine sets up the structures to include meteorological 
! forcing analyses in a coupled runmode. The routine simply allocates
! memory for the forcing variables. 
! 
! The methods invoked are:  
! \begin{description}
!  \item[create\_forcing\_structures](\ref{create_forcing_structures}) \newline
!    create and allocate memory for required structures
!  \item[forcingPerturbSetup](\ref{forcingPerturbSetup}) \newline
!    setup and allocate structures required for forcing 
!    perturbations. 
! \end{description}
!EOP

    TRACE_ENTER("metf_init")
    allocate(LIS_forc(LIS_rc%nnest,LIS_rc%nmetforc))

  ! Read Forcing Table Attributes:
    call read_forctable_attributes()
    call create_forcing_structures()
    call forcingPerturbSetup()
    TRACE_EXIT("metf_init")
  end subroutine metforcing_init_coupled


!BOP
! !ROUTINE: read_forctable_attributes
! \label{read_forctable_attributes}
!
! !INTERFACE:
  subroutine read_forctable_attributes()
! !USES:

!
! !DESCRIPTION:
! This subroutine reads in and creates forcing
! variable attribute structures from the read in  
! forcing table.
! 
!EOP
    integer             :: tnvars
    logical             :: file_exists
    integer             :: status
    type(ESMF_Config)   :: forcConfig
! _______________________________________________________

    inquire(file=LIS_rc%forcvarlistfile,exist=file_exists)

    if(.not.file_exists) then
       write(LIS_logunit,*) "[ERR] Forcing Variables list file, "
       write(LIS_logunit,*) '[ERR] ',trim(LIS_rc%forcvarlistfile),", does not exist ..."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun()
    endif

    write(LIS_logunit,*) "[INFO] Opening Forcing Variables list file: ",&
         trim(LIS_rc%forcvarlistfile)

    forcConfig = ESMF_ConfigCreate(rc=status)

    call ESMF_ConfigLoadFile(forcConfig,LIS_rc%forcvarlistfile,rc=status)

    tnvars = 0

    call ESMF_ConfigFindLabel(forcConfig,"Tair:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Tair,&
         "Near Surface Air Temperature", tnvars, status)

    call ESMF_ConfigFindLabel(forcConfig,"Qair:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Qair,&
         "Near Surface Specific Humidity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWdown:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWdown,&
         "Incident Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWdirect:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWdirect,&
         "Incident Direct Surface Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWdiffuse:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWdiffuse,&
         "Incident Diffuse Surface Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"LWdown:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_LWdown,&
         "Incident Longwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind_E:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Wind_E,&
         "Eastward Wind", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind_N:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Wind_N,&
         "Northward Wind", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Psurf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Psurf,&
         "Surface Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Rainf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Rainf,&
         "Rainfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Snowf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Snowf,&
         "Snowfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CRainf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CRainf,&
         "Convective Rainfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Forc_Hgt:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Forc_Hgt,&
         "Height of Forcing Variables", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Ch:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Ch,&
         "Surface Exchange Coefficient for Heat", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Cm:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Cm,&
         "Surface Exchange Coefficient for Momentum", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Q2sat:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Q2sat,&
         "Saturated Mixing Ratio", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Emiss:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Emiss,&
         "Surface Emissivity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Cosz:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Cosz,&
         "Cosine of Zenith Angle", tnvars,status)


    call ESMF_ConfigFindLabel(forcConfig,"Albedo:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Alb,&
         "Surface Albedo", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"PARDR:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Pardr,&
         "Surface downward PAR direct flux", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"PARDF:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Pardf,&
         "Surface downward PAR diffuse flux", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWnet:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWnet,&
         "Net downward shortwave flux", tnvars,status)

! Begin for FEWSNET
    call ESMF_ConfigFindLabel(forcConfig,"PET:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_PET,&
         "Potential ET", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"RefET:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_RefET,&
         "Reference ET", tnvars,status)
! End for FEWSNET

! CAPE available from NLDAS-2
    call ESMF_ConfigFindLabel(forcConfig,"CAPE:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CAPE,&
         "Convective Available Potential Energy", tnvars,status)

! for CRTM
    call ESMF_ConfigFindLabel(forcConfig,"LPressure:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_lpressure,&
         "Level Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"O3:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_o3,&
         "Ozone Concentration", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Xice:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_XICE,&
         "Sea Ice Mask", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"QSFC:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_QSFC,&
         "Surface Specific Humidity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CHS2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CHS2,&
         "2m Surface Exchange Coefficient for Heat", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CQS2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CQS2,&
         "2m Surface Exchange Coefficient for Moisture", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"T2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_T2,&
         "2m Air Temperature", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Q2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Q2,&
         "2m Specific Humidity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"TH2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_TH2,&
         "2m Potential Temperature", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"TMN:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_TMN,&
         "Soil Temperature at Lower Boundary", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"GVF:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_GVF,&
         "Green Vegetation Fraction", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Z0:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Z0,&
         "Roughness Length", tnvars,status)
!<for vic>
    call ESMF_ConfigFindLabel(forcConfig,"Snowflag:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SNOWFLAG,&
         "Snowflag", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Density:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_DENSITY,&
         "Atmospheric Density", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"VaporPress:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_VAPORPRESS,&
         "Vapor Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"VaporPressDeficit:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_VAPORPRESSDEFICIT,&
         "Vapor Pressure Deficit", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_WIND,&
         "Wind Speed", tnvars,status)
!</for vic>

    if(tnvars.gt.LIS_rc%nf) then
       LIS_rc%nf = tnvars
    endif

  end subroutine read_forctable_attributes


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
    integer             :: n,m,i
    character*100       :: temp
    character*1         :: nestid(2)
    integer             :: tnvars
    logical             :: file_exists
    integer             :: status
    type(ESMF_Config)   :: forcConfig

    LIS_rc%nf = 0
    do i=1, LIS_rc%nmetforc
       if(LIS_rc%nf < LIS_rc%met_nf(i)) then 
          LIS_rc%nf = LIS_rc%met_nf(i)
       endif
    enddo

    allocate(LIS_FORC_State(LIS_rc%nnest))
    allocate(LIS_FORC_Base_State(LIS_rc%nnest,LIS_rc%nmetforc))

#if 0
    inquire(file=LIS_rc%forcvarlistfile,exist=file_exists) 
    
    if(.not.file_exists) then 
       write(LIS_logunit,*) "[ERR] Forcing Variables list file, "
       write(LIS_logunit,*) '[ERR] ',trim(LIS_rc%forcvarlistfile),", does not exist ..."
       write(LIS_logunit,*) "[ERR] Program stopping ..."
       call LIS_endrun()
    endif

    write(LIS_logunit,*) "[INFO] Opening Forcing Variables list file: ",&
         trim(LIS_rc%forcvarlistfile)
    
    forcConfig = ESMF_ConfigCreate(rc=status)
    
    call ESMF_ConfigLoadFile(forcConfig,LIS_rc%forcvarlistfile,rc=status)
  
    tnvars = 0 

    call ESMF_ConfigFindLabel(forcConfig,"Tair:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Tair,&
         "Near Surface Air Temperature", tnvars, status)

    call ESMF_ConfigFindLabel(forcConfig,"Qair:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Qair,&
         "Near Surface Specific Humidity", tnvars,status)
   
    call ESMF_ConfigFindLabel(forcConfig,"SWdown:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWdown,&
         "Incident Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWdirect:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWdirect,&
         "Incident Direct Surface Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWdiffuse:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWdiffuse,&
         "Incident Diffuse Surface Shortwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"LWdown:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_LWdown,&
         "Incident Longwave Radiation", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind_E:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Wind_E,&
         "Eastward Wind", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind_N:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Wind_N,&
         "Northward Wind", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Psurf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Psurf,&
         "Surface Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Rainf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Rainf,&
         "Rainfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Snowf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Snowf,&
         "Snowfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CRainf:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CRainf,&
         "Convective Rainfall Rate", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Forc_Hgt:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Forc_Hgt,&
         "Height of Forcing Variables", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Ch:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Ch,&
         "Surface Exchange Coefficient for Heat", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Cm:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Cm,&
         "Surface Exchange Coefficient for Momentum", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Q2sat:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Q2sat,&
         "Saturated Mixing Ratio", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Emiss:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Emiss,&
         "Surface Emissivity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Cosz:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Cosz,&
         "Cosine of Zenith Angle", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Albedo:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Alb,&
         "Surface Albedo", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"PARDR:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Pardr,&
         "Surface downward PAR direct flux", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"PARDF:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Pardf,&
         "Surface downward PAR diffuse flux", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"SWnet:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SWnet,&
         "Net downward shortwave flux", tnvars,status)

! Begin for FEWSNET
    call ESMF_ConfigFindLabel(forcConfig,"PET:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_PET,&
         "Potential ET", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"RefET:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_RefET,&
         "Reference ET", tnvars,status)
! End for FEWSNET

! CAPE available from NLDAS-2
    call ESMF_ConfigFindLabel(forcConfig,"CAPE:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CAPE,&
         "Convective Available Potential Energy", tnvars,status)

! for CRTM
    call ESMF_ConfigFindLabel(forcConfig,"LPressure:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_lpressure,&
         "Level Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"O3:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_o3,&
         "Ozone Concentration", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Xice:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_XICE,&
         "Sea Ice Mask", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"QSFC:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_QSFC,&
         "Surface Specific Humidity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CHS2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CHS2,&
         "2m Surface Exchange Coefficient for Heat", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"CQS2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_CQS2,&
         "2m Surface Exchange Coefficient for Moisture", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"T2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_T2,&
         "2m Air Temperature", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Q2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Q2,&
         "2m Specific Humidity", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"TH2:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_TH2,&
         "2m Potential Temperature", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"TMN:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_TMN,&
         "Soil Temperature at Lower Boundary", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"GVF:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_GVF,&
         "Green Vegetation Fraction", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Z0:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_Z0,&
         "Roughness Length", tnvars,status)
!<for vic>
    call ESMF_ConfigFindLabel(forcConfig,"Snowflag:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_SNOWFLAG,&
         "Snowflag", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Density:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_DENSITY,&
         "Atmospheric Density", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"VaporPress:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_VAPORPRESS,&
         "Vapor Pressure", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"VaporPressDeficit:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_VAPORPRESSDEFICIT,&
         "Vapor Pressure Deficit", tnvars,status)

    call ESMF_ConfigFindLabel(forcConfig,"Wind:",rc=status)
    call get_forcingvar_attributes(forcConfig,LIS_FORC_WIND,&
         "Wind Speed", tnvars,status)
!</for vic>

    if(tnvars.gt.LIS_rc%nf) then 
       LIS_rc%nf = tnvars
    endif
#endif

    do n=1,LIS_rc%nnest       
       write(unit=temp,fmt='(i2.2)') n
       read(unit=temp,fmt='(2a1)') nestid
       
       LIS_FORC_State(n) = ESMF_StateCreate(name=&
            "Forcing State"//nestid(1)//nestid(2),&
            rc=status)
       call LIS_verify(status, &
            'error in ESMF_StateCreate:LIS_FORC_State in create_forcing_structures')       

       do m=1,LIS_rc%nmetforc
          LIS_FORC_Base_State(n,m) = ESMF_StateCreate(name=&
               "Forcing State"//nestid(1)//nestid(2),&
               rc=status)
          call LIS_verify(status,&
               'error in ESMF_StateCreate:LIS_FORC_BaseState in create_forcing_structures')       
       enddo

       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Tair)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Qair)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_SWdown)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_SWdirect)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_SWdiffuse)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_LWdown)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Wind_E)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Wind_N)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Psurf)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Rainf)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Snowf)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_CRainf)

       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Forc_Hgt)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Ch)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Cm)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Q2sat)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Emiss)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Cosz)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Alb)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Pardr)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Pardf)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_SWnet)

       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Xice)

       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_QSFC)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_CHS2)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_CQS2)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_T2)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Q2)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_TH2)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_TMN)

       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_GVF)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_Z0)

       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_lpressure)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_forc_o3)

!<for vic>
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_SNOWFLAG)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_DENSITY)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_VAPORPRESS)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_VAPORPRESSDEFICIT)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_WIND)
!</for vic>
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_PET)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_RefET)
       call add_forcing_fields(n,LIS_FORC_State(n),&
            LIS_FORC_Base_State(n,:),LIS_FORC_CAPE)
    enddo
    
    call ESMF_ConfigDestroy(forcConfig)

  end subroutine create_forcing_structures

!
!BOP
! 
! !ROUTINE: forcingPerturbSetup
! \label{forcingPerturbSetup}
! 
! !INTERFACE: 
  subroutine forcingPerturbSetup
! !USES: 
    use LIS_perturbMod 

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
    type(ESMF_Field)        :: varField(LIS_rc%nf)
    type(ESMF_Field)        :: varpertField(LIS_rc%nf)
    type(ESMF_ArraySpec)    :: arrspec1
    integer                 :: ftn 

    if(LIS_rc%perturb_forcing.ne."none") then 
       allocate(LIS_FORC_Pert_State(LIS_rc%nnest))

       do n=1,LIS_rc%nnest
          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid
          
          LIS_FORC_Pert_State(n) = ESMF_StateCreate(name=&
               "Forcing Perturbations"//&
               nestid(1)//nestid(2),  rc=status)
          call LIS_verify(status, 'force pert state create')
          write(LIS_logunit,*) '[INFO] Successfully created forcing perturbations states ..'
       enddo

       do n=1,LIS_rc%nnest
          write(LIS_logunit,*) '[INFO] Opening Forcing Attributes file: ',&
               trim(LIS_rc%forcattribfile)
          
          ftn = LIS_getNextUnitNumber()
          open(ftn,file=(LIS_rc%forcattribfile), status='old')
          read(ftn,*)
          read(ftn,*) LIS_rc%nforcepert
          read(ftn,*)
          allocate(varname(LIS_rc%nforcepert))
          allocate(varmax(LIS_rc%nforcepert))
          allocate(varmin(LIS_rc%nforcepert))
          
          allocate(forc_pert%vname(LIS_rc%nforcepert))
          allocate(forc_pert%perttype(LIS_rc%nforcepert))
          allocate(forc_pert%ssdev(LIS_rc%nforcepert))
          allocate(forc_pert%stdmax(LIS_rc%nforcepert))
          allocate(forc_pert%zeromean(LIS_rc%nforcepert))
          allocate(forc_pert%tcorr(LIS_rc%nforcepert))
          allocate(forc_pert%xcorr(LIS_rc%nforcepert))
          allocate(forc_pert%ycorr(LIS_rc%nforcepert))
          allocate(forc_pert%ccorr(LIS_rc%nforcepert,LIS_rc%nforcepert))
          
          allocate(ssdev(LIS_rc%ngrid(n)))

          do i=1,LIS_rc%nforcepert
             read(ftn,fmt='(a40)') varname(i)
             read(ftn,*) varmin(i),varmax(i)
             write(LIS_logunit,*) '[INFO] ',varname(i),varmin(i),varmax(i)
          enddo
          call LIS_releaseUnitNumber(ftn)

          call LIS_readPertAttributes(LIS_rc%nforcepert, &
               LIS_rc%forcpertattribfile,forc_pert)
      
          call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
               rc=status)
          call LIS_verify(status,'arrspecset in perturbForcing')
             
          do i=1,LIS_rc%nforcepert
             
             call ESMF_StateGet(LIS_FORC_state(n),&
                  trim(varname(i)), varField(i),rc=status)
             call LIS_verify(status,'stateget in perturbSetup')
             call ESMF_AttributeSet(varField(i),"Max Value",&
                  varmax(i),rc=status)
             call LIS_verify(status,'attributeset(max) in perturbSetup')
             call ESMF_AttributeSet(varField(i),"Min Value",&
                  varmin(i),rc=status)
             call LIS_verify(status,'attributeset(min) in perturbSetup')    
                
             varpertField(i) = ESMF_FieldCreate(grid=LIS_vecTile(n),&
                  arrayspec=arrspec1,&
                  name=trim(varname(i)),rc=status)
             call LIS_verify(status,'fieldcreate in perturbSetup')
                          
             call ESMF_StateAdd(LIS_FORC_Pert_State(n),(/varpertField(i)/),&
                  rc=status)
             call LIS_verify(status, &
                  'error in ESMF_StateAdd in perturbSetup')
          end do
          deallocate(varname)
          deallocate(varmax)
          deallocate(varmin)

          allocate(pertobjs(LIS_rc%nforcepert))
          allocate(order(LIS_rc%nforcepert))
          allocate(ccorr(LIS_rc%nforcepert,LIS_rc%nforcepert))

          call ESMF_StateGet(LIS_FORC_Pert_State(n), &
               itemNameList=pertobjs, rc=status)
          call LIS_verify(status, &
               'error in ESMF_StateGet in perturbSetup')

          order = -1
          do i=1,LIS_rc%nforcepert
             do j=1,LIS_rc%nforcepert
                if(forc_pert%vname(j).eq.pertobjs(i)) then 
                   order(i) = j
                   exit;
                endif
             enddo             
          enddo

          do i=1,LIS_rc%nforcepert             
             do j=1,LIS_rc%nforcepert
                ccorr(i,j) = forc_pert%ccorr(order(i),order(j))
             enddo
          enddo

          do i=1,LIS_rc%nforcepert
                          
             call ESMF_StateGet(LIS_FORC_Pert_State(n), &
                  pertobjs(i),pertField,&
                  rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,"Perturbation Type",&
                  forc_pert%perttype(order(i)), rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')

             if(LIS_rc%ngrid(n).gt.0) then 
                ssdev = forc_pert%ssdev(order(i))

                call ESMF_AttributeSet(pertField,"Standard Deviation",&
                     ssdev,itemCount=LIS_rc%ngrid(n),rc=status)
                call LIS_verify(status,&
                     'error in ESMF_AttributeSet in forcingPerturbSetup')
             endif


             call ESMF_AttributeSet(pertField,"Std Normal Max",&
                  forc_pert%stdmax(order(i)),rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')

             call ESMF_AttributeSet(pertField,"Ensure Zero Mean",&
                  forc_pert%zeromean(order(i)),rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,&
                  "Temporal Correlation Scale",&
                  forc_pert%tcorr(order(i)),rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,"X Correlation Scale",&
                  forc_pert%xcorr(order(i)),rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,"Y Correlation Scale",&
                  forc_pert%ycorr(order(i)),rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
             
             call ESMF_AttributeSet(pertField,&
                  "Cross Correlation Strength",&
                  ccorr(i,:),&
                  itemCount=LIS_rc%nforcepert,&
                  rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in forcingPerturbSetup')
          enddo
          deallocate(pertobjs)
          deallocate(order)
          deallocate(ccorr)
          deallocate(ssdev)
       enddo
       call perturbinit(trim(LIS_rc%perturb_forcing)//char(0), 1)
       call perturbsetup(trim(LIS_rc%perturb_forcing)//char(0),&
            1, 1, LIS_FORC_State,&
            LIS_FORC_Pert_State)

    endif
  end subroutine forcingPerturbSetup

!BOP
! !ROUTINE: LIS_get_met_forcing
! \label{LIS_get_met_forcing}
!
! !INTERFACE:
  subroutine LIS_get_met_forcing(n)

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
!  \item[retrievemetforc](\ref{retrievemetforc}) \newline
!    invokes the generic method in the registry to retrieve the 
!    met forcing data
!  \item[timeinterpmetforc](\ref{timeinterpmetforc}) \newline
!    invokes the generic method in the registry to perform
!    temporal interpolation
!  \item[LIS\_lapseRateCorrection](\ref{LIS_lapseRateCorrection}) \newline
!    method to apply lapse-rate based topographical corrections
!  \item[LIS\_slopeAspectCorrection](\ref{LIS_slopeAspectCorrection}) \newline
!    method to apply slope aspect based topographical corrections
!  \item[LIS\_microMetCorrection](\ref{LIS_microMetCorrection}) \newline
!    method to apply topographical corrections (updated with latest code)
! \end{description}
!EOP

    integer :: m

    TRACE_ENTER("metf_get")
    if(LIS_rc%nmetforc.gt.0) then 
       do m=1,LIS_rc%nmetforc

          call retrievemetforc(trim(LIS_rc%metforc(m))//char(0),n,m)
          call timeinterpmetforc(trim(LIS_rc%metforc(m))//char(0),n,m)

        ! Perform different topography-related downscaling:
          if(LIS_rc%met_ecor(m) .eq. "lapse-rate") then 

             if( LIS_rc%useelevationmap(n) == "none" ) then
                write(LIS_logunit,*) "[ERR] 'lapse-rate' correction turned on for"
                write(LIS_logunit,*) "[ERR] the forcing dataset, ",trim(LIS_rc%metforc(m)),","
                write(LIS_logunit,*) "[ERR] ... Though NO LDT-generated elevation file read in ... "
                write(LIS_logunit,*) "[ERR] -- This LIS run is ending ..."
                call LIS_endrun
             endif

             call LIS_lapseRateCorrection(n,LIS_forc(n,m)%modelelev,&
                  LIS_FORC_Base_State(n,m))

          elseif(LIS_rc%met_ecor(m).eq."slope-aspect") then 

             if( LIS_rc%useslopemap(n)  == "none" .or. & 
                 LIS_rc%useaspectmap(n) == "none" ) then
                write(LIS_logunit,*) "[ERR] 'slope-aspect' correction turned on for"
                write(LIS_logunit,*) "[ERR] the forcing dataset, ",trim(LIS_rc%metforc(m)),","
                write(LIS_logunit,*) "[ERR] ... Though NO LDT-generated slope/aspect files read in ... "
                write(LIS_logunit,*) "[ERR] -- This LIS run is ending ..."
                call LIS_endrun
             endif 
            
             call LIS_slopeAspectCorrection(n,LIS_FORC_Base_State(n,m))

          elseif(LIS_rc%met_ecor(m).eq."lapse-rate and slope-aspect") then 

             if( LIS_rc%useelevationmap(n)  == "none" .or. &
                 LIS_rc%useslopemap(n)  == "none" .or. &
                 LIS_rc%useaspectmap(n) == "none" ) then
                write(LIS_logunit,*) "[ERR] 'lapse-rate and slope-aspect' turned on for"
                write(LIS_logunit,*) "[ERR] the forcing dataset, ",trim(LIS_rc%metforc(m)),","
                write(LIS_logunit,*) "[ERR] ... Though NO LDT-generated elev/slope/aspect files read in ... "
                write(LIS_logunit,*) "[ERR] -- This LIS run is ending ..."
                call LIS_endrun
             endif

             call LIS_lapseRateCorrection(n, LIS_forc(n,m)%modelelev,&
                  LIS_FORC_Base_State(n,m))
             call LIS_slopeAspectCorrection(n, LIS_FORC_Base_State(n,m))

          ! New MicroMet option:
          elseif(LIS_rc%met_ecor(m).eq."micromet") then

             if( LIS_rc%useelevationmap(n)  == "none" .or. &
                 LIS_rc%useslopemap(n)      == "none" .or. &
                 LIS_rc%useaspectmap(n)     == "none" .or. &
                 LIS_rc%usecurvaturemap(n)  == "none" ) then

                write(LIS_logunit,*) "[ERR] 'micromet' turned on for"
                write(LIS_logunit,*) "[ERR] the forcing dataset, ",trim(LIS_rc%metforc(m)),","
                write(LIS_logunit,*) "[ERR] ... Though NO LDT-generated topo fields read in ... "
                write(LIS_logunit,*) "[ERR] This LIS run is ending ..."
                call LIS_endrun
             endif

             call LIS_MicroMetCorrection(n, LIS_forc(n,m)%modelelev,&
                  LIS_FORC_Base_State(n,m))

          elseif( LIS_rc%met_ecor(m) .eq. "micromet and slope-aspect" ) then
             write(LIS_logunit,*) "[ERR] Slope-aspect correction option not supported "
             write(LIS_logunit,*) "[ERR]  with micromet, since the micromet option accounts"
             write(LIS_logunit,*) "[ERR]  for slope, aspect and curvature corrections."
             write(LIS_logunit,*) "[ERR]  Please check lis.config.adoc for options."
             call LIS_endrun
          end if
       enddo
       
       ! Blending algorithms (overlay, forcing ensembles, bias correction..)
       if(LIS_rc%metforc_blend_alg.eq."overlay") then ! simple overlays
          call overlayForcings(n)
       elseif(LIS_rc%metforc_blend_alg.eq."ensemble") then !forcing ensembles
          call ensembleForcings(n)
       else
          write(LIS_logunit,*) "[ERR] incorrect 'Blending method for forcings'"
          write(LIS_logunit,*) "[ERR] -- This LIS run is ending ..."
          call LIS_endrun
       endif

    endif
    TRACE_EXIT("metf_get")

  end subroutine LIS_get_met_forcing

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
    character*100, allocatable :: forcobjs(:)
    type(ESMF_Field)       :: mrgField, baseField
    real, pointer          :: forcdata_base(:), forcdata_mrg(:)
    integer                :: status, status1,status2
    
    call ESMF_StateGet(LIS_FORC_State(n),itemCount=fobjcount,rc=status)
    call LIS_verify(status,'ESMF_StateGet failed for in overlayForcings')
    
    allocate(forcobjs(fobjcount))
    
    call ESMF_StateGet(LIS_FORC_State(n),itemNameList=forcobjs,rc=status)
    call LIS_verify(status,'ESMF_StateGet failed in overlayForcings')
    
    do i=1,fobjcount
       call ESMF_StateGet(LIS_FORC_State(n),forcobjs(i),mrgField,&
            rc=status)
       call LIS_verify(status, 'ESMF_StateGet failed for '//trim(forcobjs(i)))
       call ESMF_FieldGet(mrgField,localDE=0,farrayPtr= forcdata_mrg, &
            rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for in overlayForcings')
       
       do m=1,LIS_rc%nmetforc
          call ESMF_StateGet(LIS_FORC_Base_State(n,m),forcobjs(i),baseField,&
               rc=status1)
          if(status1.eq.0) then 
             call ESMF_FieldGet(baseField,localDE=0,farrayPtr= forcdata_base, &
                  rc=status2)
             call LIS_verify(status2,'ESMF_FieldGet failed in LIS_metForcingMod')
             
             do t=1,LIS_rc%ntiles(n)
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
    integer, intent(in)  :: n 
!
! !DESCRIPTION:
!  This subroutine blends different forcing inputs by assigning them 
!  to different ensemble members. 
!
!EOP
    integer                :: fobjcount
    integer                :: i,t,m,tid,tid1,tid2,index1
    character*100, allocatable :: forcobjs(:)
    type(ESMF_Field)       :: mrgField, baseField
    real, pointer          :: forcdata_base(:), forcdata_mrg(:)
    integer                :: status, status1,status2
    
    call ESMF_StateGet(LIS_FORC_State(n),itemCount=fobjcount,rc=status)
    call LIS_verify(status,'ESMF_StateGet failed for in ensembleForcings')
    
    allocate(forcobjs(fobjcount))

    call ESMF_StateGet(LIS_FORC_State(n),itemNameList=forcobjs,rc=status)
    call LIS_verify(status,'ESMF_StateGet failed in ensembleForcings')
    
    do i=1,fobjcount
       call ESMF_StateGet(LIS_FORC_State(n),forcobjs(i),mrgField,&
            rc=status)
       call LIS_verify(status, 'ESMF_StateGet failed for '//trim(forcobjs(i)))
       call ESMF_FieldGet(mrgField,localDE=0,farrayPtr= forcdata_mrg, &
            rc=status)
       call LIS_verify(status,'ESMF_FieldGet failed for in ensembleForcings')
       
       do m=1,LIS_rc%nmetforc
          call ESMF_StateGet(LIS_FORC_Base_State(n,m),forcobjs(i),baseField,&
               rc=status1)
          if(status1.eq.0) then 
             call ESMF_FieldGet(baseField,localDE=0,farrayPtr= forcdata_base, &
                  rc=status2)
             call LIS_verify(status2,&
                  'ESMF_FieldGet failed in LIS_metForcingMod')
             if(m.ne.1) then      
                do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                   tid1 =(t-1)*LIS_rc%nensem(n)+&
                        sum(LIS_rc%met_nperforc(1:m-1))+ 1
                   tid2 = tid1 + ( LIS_rc%met_nperforc(m)-1 )
                   do tid=tid1,tid2
                      if(forcdata_base(tid).ne.LIS_rc%udef) then 
                         forcdata_mrg(tid)=forcdata_base(tid)	
                      endif
                   enddo
                enddo
             else
                forcdata_mrg = forcdata_base
             endif

#if 0 
             if(m.ne.1) then      
                do t=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                   tid1 =(t-1)*LIS_rc%nensem(n)+(m-1)*LIS_rc%nperforc+ 1
                   tid2 = tid1 + ( LIS_rc%nperforc-1 )
                   do tid=tid1,tid2
                      if(forcdata_base(tid).ne.LIS_rc%udef) then 
                         forcdata_mrg(tid)=forcdata_base(tid)	
                      endif
                   enddo
                enddo
             else
                forcdata_mrg = forcdata_base
             endif
#endif
          endif
       enddo
       
    enddo

    deallocate(forcobjs)
 !   stop
  end subroutine ensembleForcings

!BOP
! !ROUTINE: LIS_perturb_forcing
! \label{LIS_perturb_forcing}
! 
! !INTERFACE: 
  subroutine LIS_perturb_forcing(n)
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
!    routine to map a forcing fields to the LIS history writer
!  \end{description}
!EOP
    integer                   :: k
    real                      :: curr_time
    integer                   :: status
    integer                   :: i, f,j,t,m,t_unpert
    integer                   :: objcount
    integer                   :: fobjcount
    integer                   :: offset1
    character*100             :: fname
    character*100,    allocatable :: forcobjs(:)
    real,             pointer :: forcvar(:)    
    type(ESMF_Field), allocatable :: forcField(:)
    character*100,    allocatable :: stateobjs(:)
    type(ESMF_Field), allocatable :: stateField(:)
    real,             pointer :: pertdata(:)
    integer                   :: perttype
    real                      :: maxval
    real                      :: minval
    real                      :: delta
    integer                   :: nxx

    TRACE_ENTER("metf_perturb")
    k = 1
    if(LIS_rc%perturb_forcing.ne."none") then 
       curr_time = float(LIS_rc%hr)*3600+60*float(LIS_rc%mn)+float(LIS_rc%ss)
       if(mod(curr_time,real(LIS_rc%pertforcinterval)).eq.0)then
          
          call perturbmethod(trim(LIS_rc%perturb_forcing)//char(0),1,n,k, &
               LIS_FORC_State(n), LIS_FORC_Pert_State(n))
          
          call ESMF_StateGet(LIS_FORC_Pert_State(n),itemCount=objcount,rc=status)
          call LIS_verify(status, &
               'error in ESMF_StateGet in LIS_perturb_forcing')
          
          call ESMF_StateGet(LIS_FORC_State(n),itemCount=fobjcount,rc=status)
          call LIS_verify(status,&
               'error in ESMF_StateGet in LIS_perturb_forcing')
          
          allocate(stateobjs(objcount))
          allocate(stateField(objcount))
          allocate(forcobjs(fobjcount))
          allocate(forcField(objcount))
          
          call ESMF_StateGet(LIS_FORC_Pert_State(n),itemNameList=stateobjs,rc=status)
          call LIS_verify(status,&
               'error in ESMF_StateGet in LIS_perturb_forcing')
          
          call ESMF_StateGet(LIS_FORC_State(n),itemNameList=forcobjs,rc=status)
          call LIS_verify(status,&
               'error in ESMF_StateGet in LIS_perturb_forcing')
          
          do i=1,objcount

             call ESMF_StateGet(LIS_FORC_Pert_State(n),stateobjs(i),stateField(i),&
                  rc=status)
             call ESMF_FieldGet(stateField(i),localDE=0,farrayPtr= pertdata, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet in LIS_perturb_forcing')
             
             call ESMF_AttributeGet(stateField(i),"Perturbation Type",&
                  perttype,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeGet in LIS_perturb_forcing')
             
             fname   = stateobjs(i)

             call ESMF_StateGet(LIS_FORC_State(n),trim(fname),forcField(i),rc=status)
             if(status.eq.0) then 
                
                call ESMF_FieldGet(forcField(i),localDE=0,farrayPtr= forcvar,rc=status)
                call LIS_verify(status,&
                     'error in ESMF_FieldGet in LIS_perturb_forcing')
                call ESMF_AttributeGet(forcField(i), "Max Value",maxval,rc=status)
                call LIS_verify(status,&
                     'error in ESMF_AttributeGet in LIS_perturb_forcing')
                
                call ESMF_AttributeGet(forcField(i), "Min Value",minval,rc=status)
                call LIS_verify(status,&
                     'error in ESMF_AttributeGet in LIS_perturb_forcing')
                
                do j=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)

                   do m=1,LIS_rc%nensem(n)
                      t = (j-1)*LIS_rc%nensem(n)+m                      
                      if(perttype.eq.0) then            
                         if(LIS_rc%pert_bias_corr.eq.1) then 
                            if(m.ne.LIS_rc%nensem(n)) then 
                               forcvar(t) = forcvar(t)+ pertdata(t)
                            endif
                         else
                            forcvar(t) = forcvar(t)+ pertdata(t)
                         endif
                      elseif(perttype.eq.1) then 
                         if(LIS_rc%pert_bias_corr.eq.1) then 
                            if(m.ne.LIS_rc%nensem(n)) then 
                               forcvar(t) = forcvar(t)*pertdata(t)
                            endif
                         else
                            forcvar(t) = forcvar(t)*pertdata(t)
                         endif
                      endif
                   enddo
                enddo

                if(LIS_rc%pert_bias_corr.eq.1) then 
                   nxx = LIS_rc%nensem(n)/LIS_rc%nmetforc
                   
                   do j=1, LIS_rc%ntiles(n)/LIS_rc%nensem(n)
                      do f=1, LIS_rc%nmetforc
                         delta    = 0.0
                         t_unpert = (j-1)*LIS_rc%nensem(n)+f*nxx
                         do m = 1, nxx-1
                            t = (j-1)*LIS_rc%nensem(n)+&
                                 (k-1)*nxx + m
                            delta = delta + (forcvar(t)-forcvar(t_unpert))
                         enddo
                         delta = delta/(nxx-1)
                         do m=1, nxx-1
                            t = (j-1)*LIS_rc%nensem(n)+&
                                 (f-1)*nxx + m
                            forcvar(t) = forcvar(t) - delta
                         enddo
                      enddo
                   enddo
                endif

                do t=1,LIS_rc%ntiles(n)
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
    call diagnoseForcingOutput(n)
    TRACE_EXIT("metf_perturb")

 end subroutine LIS_perturb_forcing

!BOP
! !ROUTINE: LIS_metforcing_reset
!  \label{LIS_metforcing_reset}
! 
! !INTERFACE: 
  subroutine LIS_metforcing_reset
! !USES: 
! 
! !DESCRIPTION: 
!  This subroutine issues the invocation to reset variables for 
!  specific forcing scheme. 
!
!EOP
    integer       :: m

    TRACE_ENTER("metf_reset")
    if(LIS_rc%nmetforc.gt.0) then 
       do m=1,LIS_rc%nmetforc
          call resetmetforc(trim(LIS_rc%metforc(m))//char(0),m)
       enddo
    endif
    TRACE_EXIT("metf_reset")

  end subroutine LIS_metforcing_reset
!BOP
! !ROUTINE: LIS_metforcing_finalize
! \label{LIS_metforcing_finalize}
! 
! !INTERFACE: 
  subroutine LIS_metforcing_finalize
! !USES:

!
! !DESCRIPTION:
!  This routine issues the invocation to deallocate and cleanup
!  any allocated data structures in the specific instance of the 
!  forcing scheme. 
! 
! The methods invoked are:  
! \begin{description}
!  \item[finalmetforc](\ref{finalmetforc}) \newline
!    invokes the generic method in the registry to cleanup 
!    the allocated structures for the met forcing scheme. 
! \end{description}
!EOP
    integer       :: m

    if(LIS_rc%nmetforc.gt.0) then 
       do m=1,LIS_rc%nmetforc
          call finalmetforc(trim(LIS_rc%metforc(m))//char(0),m)
       enddo
       deallocate(LIS_FORC_State)
          
       if(LIS_rc%perturb_forcing.ne."none") then 
          deallocate(LIS_FORC_Pert_State)
       endif
    endif

  end subroutine LIS_metforcing_finalize

!BOP
! 
!  !ROUTINE: get_forcingvar_attributes
!  \label{get_forcingvar_attributes}
! 
! !INTERFACE: 
  subroutine get_forcingvar_attributes(forcConfig, forc_attrib, &
       forcname, varcount, &
       checkentry)
! !USES:     

    implicit none
! !ARGUMENTS:    
    type(ESMF_Config)      :: forcConfig
    type(forc_attrib_type) :: forc_attrib 
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
    character(len=3)           :: fid
    integer                    :: k 
    integer                    :: status
    
    if(checkentry.eq.0) then ! the forcing field specification is found 
       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%selectOpt, &
            default=0,rc=status)
       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%vlevels,&
            default=0,rc=status)
       call ESMF_ConfigGetAttribute(forcConfig,forc_attrib%units, rc=status)
       
       if(forc_attrib%selectOpt.eq.1) then 
          allocate(forc_attrib%varname(forc_attrib%vlevels))
          do k=1,forc_attrib%vlevels
             write(fid,fmt='(i3.3)') k
             forc_attrib%varname(k) = trim(forcname)//' Level '//trim(fid)
          enddo
          write(LIS_logunit,*) '[INFO] FORCING: ',forcname, forc_attrib%vlevels, &
               forc_attrib%units
          varcount = varcount + forc_attrib%vlevels
       endif       
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
  subroutine add_forcing_fields(n, FORC_State, FORC_Base_State,forc_attrib)
! !USES:     

    implicit none
! !ARGUMENTS: 
    integer                   :: n 
    type(ESMF_State)          :: FORC_State
    type(ESMF_State)          :: FORC_Base_State(LIS_rc%nmetforc)
    type(forc_attrib_type)    :: forc_attrib
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
    real,           pointer   :: var(:)
    integer                   :: status

    if(forc_attrib%selectOpt.eq.1) then 
       call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
            rc=status)
       call LIS_verify(status,&
            'error in ESMF_ArraySpecSet in add_forcing_fields')       

       do k=1,forc_attrib%vlevels
          varField_m = ESMF_FieldCreate(grid=LIS_vecTile(n), &
               arrayspec=arrspec1, name = trim(forc_attrib%varname(k)),&
               rc=status)
          call LIS_verify(status,&
               'error in ESMF_FieldCreate in add_forcing_fields')
          
          call ESMF_AttributeSet(varField_m,"Units",trim(forc_attrib%units),&
               rc=status)
          call LIS_verify(status,&
               'error in ESMF_AttributeSet in add_forcing_fields')
          
          call ESMF_FieldGet(varField_m,localDE=0,farrayPtr=var,rc=status)
          call LIS_verify(status,&
               'error in ESMF_FieldGet in add_forcing_fields')
          var = LIS_rc%udef

          call ESMF_StateAdd(FORC_State,(/varField_m/),rc=status)
          call LIS_verify(status,&
               'error in ESMF_StateAdd in add_forcing_fields')

          do m=1,LIS_rc%nmetforc

             varField_b = ESMF_FieldCreate(grid=LIS_vecTile(n), &
               arrayspec=arrspec1, name = trim(forc_attrib%varname(k)),&
               rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldCreate in add_forcing_fields')
             
             call ESMF_AttributeSet(varField_b,"Units",trim(forc_attrib%units),&
                  rc=status)
             call LIS_verify(status,&
                  'error in ESMF_AttributeSet in add_forcing_fields')
             
             call ESMF_FieldGet(varField_b,localDE=0,farrayPtr=var,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet in add_forcing_fields')
             var = LIS_rc%udef
             
             call ESMF_StateAdd(FORC_Base_State(m),(/varField_b/),rc=status)
             call LIS_verify(status,&
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
    use LIS_histDataMod
    use LIS_constantsMod,   only : LIS_MS2KMDAY
! !ARGUMENTS: 
    integer,   intent(IN)  :: n 
! 
! !DESCRIPTION: 
!  This routines issues the diagnose calls to map the forcing variables 
!  to the history output writer. 
! 
!   The routines invoked are: 
!  \begin{description}
!  \item[LIS\_diagnoseOutputVar](\ref{LIS_diagnoseSurfaceOutputVar}) \newline
!    generic routine to map a single variable to the LIS 
!    history writer
!  \end{description}
!
!EOP
    integer            :: status, t
    type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
    type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
    real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
    real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
    integer            :: k 
    type(ESMF_Field)   :: swdirField,swdifField,hField,chField,cmField
    type(ESMF_Field)   :: emissField,mixField,coszField,albField
    type(ESMF_Field)   :: tempField
    real,pointer       :: swdir(:),swdif(:),harray(:),charray(:),cmarray(:)
    real,pointer       :: emiss(:),mix(:),cosz(:),alb(:)
    real,pointer       :: tempPtr(:)
    real               :: windmag, windmag2

    real, allocatable  :: totpcp(:)
! ___________________________________________________________________________

    k = 1 ! since forcing is same for all DA instances. 

! Currently this is set to write only the surface forcing variables and not 
! the entire profile

    if(LIS_rc%wout.ne."none") then 
       if(LIS_FORC_Tair%selectOpt.eq.1) then
          do k=1,LIS_FORC_Tair%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Tair%varname(k)),&
                  tmpField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Tair in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(tmpField, localDE=0, farrayPtr= tmp,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Tair in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_TAIRFORC,value=&
                     tmp(t),vlevel=k,unit="K",direction="-",&
                     valid_min = 213.0, valid_max=333.0)
             enddo
          enddo
       endif

       if(LIS_FORC_Qair%selectOpt.eq.1) then
          do k=1,LIS_FORC_Qair%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Qair%varname(k)),&
                  q2Field,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Qair in diagnoseForcingOutput')
             
             call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Qair in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_QAIRFORC,value=&
                     q2(t),vlevel=k,unit="kg kg-1",direction="-",&
                     valid_min = 0.0, valid_max=0.03)
             enddo
          enddo
       endif

       if(LIS_FORC_SWdown%selectOpt.eq.1) then
          do k=1,LIS_FORC_SWdown%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SWdown%varname(k)),&
                  swdField,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:SWdown in diagnoseForcingOutput')

             call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:SWdown in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_SWDOWNFORC,value=&
                     swd(t),vlevel=k,unit="W m-2",direction="DN",&
                     valid_min = 0.0, valid_max=1360.0)
             enddo
          enddo
       endif

       if(LIS_FORC_LWdown%selectOpt.eq.1) then
          do k=1,LIS_FORC_LWdown%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_LWdown%varname(k)),&
                  lwdField,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:LWdown in diagnoseForcingOutput')

             call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:LWdown in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_LWDOWNFORC,value=&
                     lwd(t),vlevel=k,unit="W m-2",direction="DN",&
                     valid_min = 0.0, valid_max = 750.0)
             enddo
          enddo
       endif

       if(LIS_FORC_Wind_E%selectOpt.eq.1.and.LIS_FORC_Wind_N%selectOpt.eq.1) then  
          do k=1,LIS_FORC_Wind_E%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_E%varname(k)),&
                  uField,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Wind_E in diagnoseForcingOutput')

             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_N%varname(k)),&
                  vField,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_stateGet:Wind_N in diagnoseForcingOutput')
             
             call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Wind_E in diagnoseForcingOutput')
             
             call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Wind_N in diagnoseForcingOutput')             
             
             windmag = LIS_rc%udef
             windmag2 = LIS_rc%udef
             do t=1,LIS_rc%ntiles(n)

              ! Make sure undefined values are not factored in to wind-mag calc:
                if( uwind(t) == LIS_rc%udef .or. vwind(t) == LIS_rc%udef ) then
                   windmag  = LIS_rc%udef
                   windmag2 = LIS_rc%udef
                else
                   windmag  = sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t))
                   windmag2 = (sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t))*LIS_MS2KMDAY)
                endif

                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_WINDFORC,value=&
!                     sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t)),vlevel=k,&
                     windmag, vlevel=k,&
                     unit="m s-1",direction="-")

                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_WINDFORC,value=&
!                     (sqrt(uwind(t)*uwind(t)+vwind(t)*vwind(t))*LIS_MS2KMDAY),&
                     windmag2, &
                     vlevel=k,unit="km day-1",direction="-")
             enddo

          enddo
       endif

       if(LIS_FORC_Psurf%selectOpt.eq.1) then
          do k=1,LIS_FORC_Psurf%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Psurf%varname(k)),&
                  psurfField,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Psurf in diagnoseForcingOutput')

             call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Psurf in diagnoseForcingOutput')
             
             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_PSURFFORC,value=&
                     psurf(t),vlevel=k,unit="Pa",direction="-", &
                     valid_min = 5000.0, valid_max=110000.0)
             enddo
          enddo
       endif

       if(LIS_FORC_Rainf%selectOpt.eq.1 .or. LIS_FORC_Snowf%selectOpt.eq.1) then

          allocate(totpcp(LIS_rc%ntiles(n)))

          do k=1,LIS_FORC_Rainf%vlevels

             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Rainf%varname(k)),&
                  pcpField,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Rainf in diagnoseForcingOutput')
             
             call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Rainf in diagnoseForcingOutput')

             do t=1,lis_rc%ntiles(n)
                totpcp(t) = pcp(t)
             enddo

             if(LIS_FORC_Rainf%selectOpt.eq.1) then
                do t=1,LIS_rc%ntiles(n)
                   call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINFFORC, &
                          value=pcp(t),vlevel=k,unit="kg m-2 s-1",direction="DN", &
                          valid_min=0.0,valid_max=0.02)
                   call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_RAINFFORC, &
                          value=pcp(t)*LIS_rc%ts,vlevel=k,unit="kg m-2",direction="DN")
                enddo
             endif

             if(LIS_FORC_Snowf%selectOpt.eq.1) then
                call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Snowf%varname(k)),&
                     snowfField,rc=status)
                call LIS_verify(status,&
                     'error in ESMF_StateGet:Snowf in diagnoseForcingOutput')
               
                call ESMF_FieldGet(snowfField,localDE=0, farrayPtr=snowf,rc=status)
                call LIS_verify(status,&
                     'error in ESMF_FieldGet:Snowf in diagnoseForcingOutput')

                do t=1,lis_rc%ntiles(n)
                   if( snowf(t) >= 0. ) then
                      totpcp(t) = totpcp(t) + snowf(t)
                   endif
                enddo

                do t=1,LIS_rc%ntiles(n)
                   call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_SNOWFFORC, &
                          value=snowf(t),vlevel=k,unit="kg m-2 s-1",direction="DN", &
                          valid_min=0.0,valid_max=0.02)
                   call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_SNOWFFORC, &
                          value=snowf(t)*LIS_rc%ts,vlevel=k,unit="kg m-2",direction="DN")
                enddo

             endif

             do t=1,LIS_rc%ntiles(n)

                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_TOTALPRECIP,value=&
                     totpcp(t),vlevel=k,unit="kg m-2 s-1",direction="DN",&
                     valid_min = 0.0, valid_max=0.02)
                call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_TOTALPRECIP,               &
                     value=totpcp(t)*LIS_rc%ts,vlevel=1,unit="kg m-2",direction="DN")
             enddo
             
          enddo

          deallocate(totpcp)
       endif


       if(LIS_FORC_CRainf%selectOpt.eq.1) then
          do k=1,LIS_FORC_CRainf%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_CRainf%varname(k)),&
                  cpcpField,rc=status)
             call LIS_verify(status)

             call ESMF_FieldGet(cpcpField,localDE=0, farrayPtr=cpcp,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_CRAINFFORC,value=&
                     cpcp(t),vlevel=k,unit="kg m-2 s-1",direction="DN",&
                     valid_min = 0.0, valid_max = 0.02)
                call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_CRAINFFORC, &
                       value=cpcp(t)*LIS_rc%ts,vlevel=k,unit="kg m-2",direction="DN")
             enddo
          enddo
       endif

       if(LIS_FORC_Wind_N%selectOpt.eq.1) then
          do k=1,LIS_FORC_Wind_N%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_N%varname(k)),&
                  vField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Wind_N in diagnoseForcingOutput')    
             
             call ESMF_FieldGet(vField, localDE=0, farrayPtr=vwind,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Wind_N in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_NWINDFORC,value=&
                     vwind(t),vlevel=k,unit="m s-1",direction="N",&
                     valid_min = -75.0, valid_max = 75.0)
             enddo
          enddo
       endif

       if(LIS_FORC_Wind_E%selectOpt.eq.1) then
          do k=1,LIS_FORC_Wind_E%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Wind_E%varname(k)),&
                  uField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Wind_E in diagnoseForcingOutput')    
             
             call ESMF_FieldGet(uField, localDE=0, farrayPtr=uwind,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Wind_E in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_EWINDFORC,value=&
                     uwind(t),vlevel=k,unit="m s-1",direction="E",&
                     valid_min = -75.0, valid_max = 75.0)
             enddo
          enddo
       endif

       if(LIS_FORC_SWdirect%selectOpt.eq.1) then
          do k=1,LIS_FORC_SWdirect%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SWdirect%varname(k)),&
                  swdirField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:SWdirect in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(swdirField, localDE=0, farrayPtr=swdir,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:SWdirect in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_DIRECTSWFORC,value=&
                                    swdir(t),vlevel=k,unit="W m-2",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_SWdiffuse%selectOpt.eq.1) then
          do k=1,LIS_FORC_SWdiffuse%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SWdiffuse%varname(k)),&
                  swdifField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:SWdiffuse in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(swdifField, localDE=0, farrayPtr=swdif,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:SWdiffuse in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_DIFFUSESWFORC,value=&
                                    swdif(t),vlevel=k,unit="W m-2",direction="-")
             enddo
          enddo
       endif


       if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
          do k=1,LIS_FORC_Forc_Hgt%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Forc_Hgt%varname(k)),&
                  hField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Forc_Hgt in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(hField, localDE=0, farrayPtr=harray,rc=status)
             call LIS_verify(status,'error in ESMF_FieldGet:Forc_Hgt in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_FHEIGHTFORC,value=&
                                      harray(t),vlevel=k,unit="m",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_Ch%selectOpt.eq.1) then
          do k=1,LIS_FORC_Ch%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Ch%varname(k)),&
                  chField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Ch in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(chField, localDE=0, farrayPtr=charray,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Ch in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_CHFORC,value=&
                                     charray(t),vlevel=k,unit="m s-1",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_Cm%selectOpt.eq.1) then
          do k=1,LIS_FORC_Cm%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Cm%varname(k)),&
                  cmField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Cm in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(cmField, localDE=0, farrayPtr=cmarray,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Cm in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_CMFORC,value=&
                                     cmarray(t),vlevel=k,unit="m s-1",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_Emiss%selectOpt.eq.1) then
          do k=1,LIS_FORC_Emiss%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Emiss%varname(k)),&
                  emissField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Emiss in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(emissField, localDE=0, farrayPtr=emiss,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Emiss in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_EMISSFORC,value=&
                                       emiss(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_Q2sat%selectOpt.eq.1) then
          do k=1,LIS_FORC_Q2sat%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Q2sat%varname(k)),&
                  mixField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Q2sat in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(mixField, localDE=0, farrayPtr=mix,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Q2sat in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_MIXRATIOFORC,value=&
                                     mix(t),vlevel=k,unit="kg kg-1",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_Cosz%selectOpt.eq.1) then
          do k=1,LIS_FORC_Cosz%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Cosz%varname(k)),&
                  coszField, rc=status)
             call LIS_verify(status,&
                  'error ESMF_StateGet:Cosz in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(coszField, localDE=0, farrayPtr=cosz,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Cosz in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_COSZENFORC,value=&
                                        cosz(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_Alb%selectOpt.eq.1) then
          do k=1,LIS_FORC_Alb%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Alb%varname(k)),&
                  albField, rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:Alb in diagnoseForcingOutput')         
             
             call ESMF_FieldGet(albField, localDE=0, farrayPtr=alb,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:Alb in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n,t,LIS_MOC_ALBEDOFORC,value=&
                                         alb(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_Pardr%selectOpt.eq.1) then
          do k=1,LIS_FORC_Pardr%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Pardr%varname(k)),&
                  q2Field,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:PARDR in diagnoseForcingOutput')
             
             call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:PARDR in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_PARDRFORC,value=&
                     q2(t),vlevel=k,unit="W m-2",direction="DN",&
                     valid_min = 0.0, valid_max=1360.0)
             enddo
          enddo
       endif

       if(LIS_FORC_Pardf%selectOpt.eq.1) then
          do k=1,LIS_FORC_Pardf%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_Pardf%varname(k)),&
                  q2Field,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_StateGet:PARDF in diagnoseForcingOutput')
             
             call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
             call LIS_verify(status,&
                  'error in ESMF_FieldGet:PARDF in diagnoseForcingOutput')

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_PARDFFORC,value=&
                     q2(t),vlevel=k,unit="W m-2",direction="DN",&
                     valid_min = 0.0, valid_max=1360.0)
             enddo
          enddo
       endif

!<for vic>
       if(LIS_FORC_SNOWFLAG%selectOpt.eq.1) then
          do k=1,LIS_FORC_SNOWFLAG%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_SNOWFLAG%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_SNOWFLAGFORC,value=&
                                     tempPtr(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_DENSITY%selectOpt.eq.1) then
          do k=1,LIS_FORC_DENSITY%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_DENSITY%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_DENSITYFORC,value=&
                                     tempPtr(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_VAPORPRESS%selectOpt.eq.1) then
          do k=1,LIS_FORC_VAPORPRESS%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_VAPORPRESS%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_VAPORPRESSFORC,value=&
                                     tempPtr(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_VAPORPRESSDEFICIT%selectOpt.eq.1) then
          do k=1,LIS_FORC_VAPORPRESSDEFICIT%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_VAPORPRESSDEFICIT%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_VAPORPRESSDEFICITFORC,value=&
                                     tempPtr(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_WIND%selectOpt.eq.1) then
          do k=1,LIS_FORC_WIND%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_WIND%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_WINDFORC,value=&
                                     tempPtr(t),vlevel=k,unit="-",direction="-")
             enddo
          enddo
       endif
!</for vic>
       ! SY: Begin for FEWSNET
       if(LIS_FORC_PET%selectOpt.eq.1) then
          do k=1,LIS_FORC_PET%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_PET%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_PETFORC,value=&
                     tempPtr(t),vlevel=k,unit="kg m-2 s-1",direction="-")
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_PETFORC,value=&
                     (tempPtr(t)*LIS_rc%ts),vlevel=k,unit="kg m-2",direction="-")
             enddo
          enddo
       endif

       if(LIS_FORC_RefET%selectOpt.eq.1) then
          do k=1,LIS_FORC_RefET%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_RefET%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_REFETFORC,value=&
                     tempPtr(t),vlevel=k,unit="kg m-2",direction="-")
             enddo
          enddo
       endif
       ! SY: End for FEWSNET

       if(LIS_FORC_CAPE%selectOpt.eq.1) then
          do k=1,LIS_FORC_CAPE%vlevels
             call ESMF_StateGet(LIS_FORC_State(n),(LIS_FORC_CAPE%varname(k)),&
                  tempField, rc=status)
             call LIS_verify(status)         
             
             call ESMF_FieldGet(tempField, localDE=0, farrayPtr=tempPtr,rc=status)
             call LIS_verify(status)

             do t=1,LIS_rc%ntiles(n)
                call LIS_diagnoseSurfaceOutputVar(n, t,LIS_MOC_CAPEFORC,value=&
                     tempPtr(t),vlevel=k,unit="J kg-1",direction="-")
             enddo
          enddo
       endif

    endif
  end subroutine diagnoseForcingOutput

end module LIS_metforcingMod

