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
module CRTM2_EMMod
!BOP
!
! !MODULE: CRTM2_EMMod
!
! !DESCRIPTION:
!    This module provides the routines to control the execution of 
!    CRTM 2.x from within LIS. The subroutines provide mapping of the surface
!    and atmospheric profiles to CRTM2, which then runs the forward 
!    model to generate the radiative quantities. 
!
! !REVISION HISTORY:
! 20 Oct 2010: Yudong Tian; Modifed from CRTM2_handlerMod.F90 to support 
!  CRTM2EM (Emissivity only). 
!
! !USES:        

#if (defined RTMS)

  USE CRTM_Module
  USE CRTM_GeometryInfo_Define,  ONLY: CRTM_GeometryInfo_type
  USE CRTM_SfcOptics
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: CRTM2_EMonly_initialize
  public :: CRTM2_EMonly_f2t
  public :: CRTM2_EMonly_geometry
  public :: CRTM2_EMonly_run
  public :: CRTM2_EMonly_output
  public :: CRTM_landmatch
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: crtm_struc, sm_correction_struc
!EOP
  type, public ::  crtm_type_dec 
     integer                               :: nsensors
     integer                               :: nprofiles
     integer                               :: nlayers
     integer                               :: nabsorbers
     integer                               :: nclouds
     integer                               :: naerosols
     character*256                         :: sensor_id
     character*100                         :: coeff_data
     real                                  :: zenith_angle
     integer				   :: SensorIndex
     integer				   :: ChannelIndex
     type(CRTM_ChannelInfo_type),  allocatable :: ChannelInfo(:)
     type(CRTM_Geometry_type), allocatable :: Geometry(:)
     type(CRTM_SfcOptics_type),    allocatable :: SfcOptics(:, :)
     type(CRTM_SOVariables_type),  allocatable :: SOV(:)
! forward variables
     type(CRTM_Surface_type),      allocatable :: Sfc(:)
     real, allocatable                         :: emissivity_ave(:,:)
     real, allocatable                         :: emissivity(:,:)
     real, allocatable                         :: emissivity_bg(:,:)  !bare ground

     real, allocatable                         :: bgf_fixed(:)     ! bare ground fraction
     real, allocatable                         :: bgf_total(:)         ! bare ground fraction total for grid cell
     real, allocatable                         :: bgf_v(:)         ! bare ground fraction total for veg portion of grid cell
     real, allocatable                         :: k_lai2vgf(:)  ! gvf of vegetated portion = 1*exp(-k_lai2vgf*LAI) 

  end type crtm_type_dec

  type, public ::  sm_correction_dec 
     integer                               :: c_type 
     character(len=LIS_CONST_PATH_LEN)     :: src_mean_file
     character(len=LIS_CONST_PATH_LEN)     :: src_sigma_file
     character(len=LIS_CONST_PATH_LEN)     :: dst_mean_file
     character(len=LIS_CONST_PATH_LEN)     :: dst_sigma_file
     real                                  :: gridDesc(8)
     real, allocatable                         :: src_mean(:, :), src_sigma(:, :)  
     real, allocatable                         :: dst_mean(:, :), dst_sigma(:, :)  
  end type sm_correction_dec

  type(crtm_type_dec), allocatable :: crtm_struc(:)
  type(sm_correction_dec), allocatable :: sm_correction_struc(:) 
  SAVE

contains
!BOP
! 
! !ROUTINE: CRTM2_EMonly_initialize
! \label{CRTM2_EMonly_initialize}
! 
! !INTERFACE:
  subroutine CRTM2_EMonly_initialize()
! !USES:
!    use CRTM_Module
    use LIS_coreMod,    only : LIS_rc, LIS_config
    use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_getNextUnitNumber, &
                               LIS_releaseUnitNumber
    use LIS_RTMMod,     only : LIS_sfcState
    use LIS_fileIOMod,  only : LIS_readData, LIS_readDomainConfigSpecs

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for noah-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noah from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readCRTM2EMcrd](\ref{readCRTM2EMcrd}) \newline
!    reads the runtime options for CRTM2 EMonly
!  \end{description}
!EOP
    implicit none

    integer :: n , ftn, rc
    integer :: m 
    integer :: n_Sensors, n_Channels, n_Profiles, nc, np
    integer :: status
    real, allocatable ::  tmp_gridDesc(:, :)

    allocate(tmp_gridDesc(LIS_rc%nnest, 8))
    allocate(crtm_struc(LIS_rc%nnest))
    allocate(sm_correction_struc(LIS_rc%nnest))

    call ESMF_ConfigFindLabel(LIS_config,"RTM input soil moisture correction:",rc=rc)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,sm_correction_struc(n)%c_type,rc=rc)
    end do 

    do n=1,LIS_rc%nnest
       
       if (sm_correction_struc(n)%c_type .eq. 1) then 
          call ESMF_ConfigFindLabel(LIS_config,"RTM input soil moisture correction src mean file:",rc=rc)
          call ESMF_ConfigGetAttribute(LIS_config,sm_correction_struc(n)%src_mean_file,rc=rc)
          
          call ESMF_ConfigFindLabel(LIS_config,"RTM input soil moisture correction src sigma file:",rc=rc)
          call ESMF_ConfigGetAttribute(LIS_config,sm_correction_struc(n)%src_sigma_file,rc=rc)
          
          call ESMF_ConfigFindLabel(LIS_config,"RTM input soil moisture correction dst mean file:",rc=rc)
          call ESMF_ConfigGetAttribute(LIS_config,sm_correction_struc(n)%dst_mean_file,rc=rc)
          
          call ESMF_ConfigFindLabel(LIS_config,"RTM input soil moisture correction dst sigma file:",rc=rc)
          call ESMF_ConfigGetAttribute(LIS_config,sm_correction_struc(n)%dst_sigma_file,rc=rc)
          
          call LIS_readDomainConfigSpecs("RTM input soil moisture correction", tmp_gridDesc)
       endif
    enddo

    do n=1,LIS_rc%nnest
       if (sm_correction_struc(n)%c_type .eq. 1) then 
          sm_correction_struc(n)%gridDesc(:) = tmp_gridDesc(n, :)
          allocate(sm_correction_struc(n)%src_mean(LIS_rc%lnc(n),LIS_rc%lnr(n))) 
          allocate(sm_correction_struc(n)%src_sigma(LIS_rc%lnc(n),LIS_rc%lnr(n))) 
          allocate(sm_correction_struc(n)%dst_mean(LIS_rc%lnc(n),LIS_rc%lnr(n))) 
          allocate(sm_correction_struc(n)%dst_sigma(LIS_rc%lnc(n),LIS_rc%lnr(n))) 
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%src_mean_file, access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc, sm_correction_struc(n)%src_mean) 
          call LIS_releaseUnitNumber(ftn)
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%src_sigma_file, access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc, sm_correction_struc(n)%src_sigma) 
          call LIS_releaseUnitNumber(ftn)
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%dst_mean_file, access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc, sm_correction_struc(n)%dst_mean)
          call LIS_releaseUnitNumber(ftn)
          
          ftn = LIS_getNextUnitNumber()
          open(ftn, file=sm_correction_struc(n)%dst_sigma_file, access="direct",  &
               form="unformatted", recl=4)
          call LIS_readData(n, ftn, sm_correction_struc(n)%gridDesc, sm_correction_struc(n)%dst_sigma)
          call LIS_releaseUnitNumber(ftn)
          
       end if
    end do 

    call readCRTM2EMcrd()

    do n=1,LIS_rc%nnest
       crtm_struc(n)%nprofiles = LIS_rc%ntiles(n) ! all grid points in LIS

       allocate(crtm_struc(n)%ChannelInfo(crtm_struc(n)%nsensors))
       allocate(crtm_struc(n)%Geometry(crtm_struc(n)%nprofiles))
       allocate(crtm_struc(n)%Sfc(crtm_struc(n)%nprofiles))
       allocate(crtm_struc(n)%SOV(crtm_struc(n)%nprofiles))
!---------------------------------------------------------------------
!  This initializes the CRTM for the sensors defined through 
!  the configuration file. 
!---------------------------------------------------------------------

       write(LIS_logunit,*) ' Initializing CRTM .... '
       status = CRTM_Init((/crtm_struc(n)%Sensor_Id/), & 
            crtm_struc(n)%ChannelInfo, &
            File_Path = trim(crtm_struc(n)%coeff_data))

       call LIS_verify(status, 'Error Initializing CRTM')
!---------------------------------------------------------------------
! Determine the total number of channels for which CRTM was
! initialized
!---------------------------------------------------------------------
       n_Sensors  = SIZE(crtm_struc(n)%ChannelInfo)
       n_Channels = SUM(crtm_struc(n)%ChannelInfo%n_Channels)
       n_Profiles = crtm_struc(n)%nprofiles

      allocate(crtm_struc(n)%SfcOptics(n_Channels, n_Profiles))
      allocate(crtm_struc(n)%emissivity_ave(n_Channels, n_Profiles))
      allocate(crtm_struc(n)%emissivity(n_Channels, n_Profiles))
      allocate(crtm_struc(n)%emissivity_bg(n_Channels, n_Profiles))

      Do np=1, n_Profiles
       Do nc=1, n_Channels 
        status=CRTM_Allocate_SfcOptics( MAX_N_ANGLES, &  ! Input
                                              MAX_N_STOKES, &  ! Input
                                      crtm_struc(n)%SfcOptics(nc, np)   )  ! Output

        call LIS_verify(status, 'Error in CRTM_Allocate_SfcOptics')

       End Do
      End Do

      IF ( n_Sensors == 0 .OR. n_Channels == 0 ) then
         write(*, *)"Error: n_Sensors == 0 .OR. n_Channels == 0!"
         RETURN
      End If

!---------------------------------------------------------------------
! Allocate 'tiled' approach of bg and lai fixes
!---------------------------------------------------------------------
      allocate( crtm_struc(n)%bgf_fixed     (crtm_struc(n)%nprofiles) ) 
      allocate( crtm_struc(n)%bgf_total     (crtm_struc(n)%nprofiles) ) 
      allocate( crtm_struc(n)%bgf_v         (crtm_struc(n)%nprofiles) ) 
      allocate( crtm_struc(n)%k_lai2vgf     (crtm_struc(n)%nprofiles) ) 

      crtm_struc(n)%k_lai2vgf = 0.52  ! redefined now as exponent for lai-gvf relationship (formerly bare ground fraction) 
      crtm_struc(n)%bgf_fixed = 0.01  ! default is zero bare ground fraction
      
!---------------------------------------------------------------------
! Allocate ESMF State for mapping surface properties
!---------------------------------------------------------------------

       call add_sfc_fields(n,LIS_sfcState(n), "Wind Speed")
       call add_sfc_fields(n,LIS_sfcState(n), "Land Coverage")
       call add_sfc_fields(n,LIS_sfcState(n), "Land Type")
       call add_sfc_fields(n,LIS_sfcState(n), "Land Temperature")
       call add_sfc_fields(n,LIS_sfcState(n), "Snow Temperature")
       call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Content")
       call add_sfc_fields(n,LIS_sfcState(n), "Soil Temperature")
       call add_sfc_fields(n,LIS_sfcState(n), "Canopy Water Content")
       call add_sfc_fields(n,LIS_sfcState(n), "Vegetation Fraction")
       call add_sfc_fields(n,LIS_sfcState(n), "Snow Coverage")
       call add_sfc_fields(n,LIS_sfcState(n), "Snow Depth")
       call add_sfc_fields(n,LIS_sfcState(n), "Snow Density")
       call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")

    enddo


  end subroutine CRTM2_EMonly_initialize 


  subroutine add_sfc_fields(n, sfcState,varname)

    use LIS_logMod,   only : LIS_verify
    use LIS_coreMod,  only : LIS_vecTile

    implicit none 

    integer            :: n 
    type(ESMF_State)   :: sfcState
    character(len=*)   :: varname

    type(ESMF_Field)     :: varField
    type(ESMF_ArraySpec) :: arrspec
    integer              :: status
    
    call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    varField = ESMF_FieldCreate(arrayspec=arrSpec, & 
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
    call LIS_verify(status, 'Error in field_create of '//trim(varname))
    
    call ESMF_StateAdd(sfcState, (/varField/), rc=status)
    call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

  end subroutine add_sfc_fields

  subroutine CRTM2_EMonly_f2t(n)

    use LIS_coreMod,         only : LIS_rc
    use LIS_metforcingMod,  only : LIS_FORC_State
    use LIS_logMod,          only : LIS_verify
!    USE Profile_Utility, only: sa_to_mr, mr_to_ppmv
    USE Units_Conversion
    use LIS_FORC_AttributesMod

    implicit none

    integer, intent(in)    :: n 

    integer                :: k, t, kk
    integer                :: status
    type(ESMF_Field)       :: lpressField, pressField, tmpField, q2Field,o3Field
    real,   allocatable        :: lpress(:), press(:), tmp(:), q2(:),o3(:)
    real                   :: qsmall, o3small 
    
    qsmall = 1E-6
    o3small = 1E-10

  end subroutine CRTM2_EMonly_f2t


  subroutine CRTM2_EMonly_geometry(n)

    implicit none

    integer,   intent(in)    :: n
    real, parameter :: PI = 3.14159265
    real                     :: r  !distance ratio (see CRTM User Guide; Geometry Structure)
    real                     :: h  !altititude of satellite
    real                     :: EarthRadius !can vary
    integer :: np


  crtm_struc(n)%Geometry(:)%Sensor_Zenith_Angle = crtm_struc(n)%zenith_angle 

! See CRTM documentation on Geometry Structure
! Sensor Scan angle is CRTM input; here computed by following assumptions
! (If we work with satellite files, may be supplied by them)
! Here Earth Radius at surface pt  below satellite assumed equal to that of  Field of View (FOV)
! Height of satellite assumed fixed, but would vary with satellite
! Testing with AMSRE 10.7GHz shows scan angle *not* impacting any Tb calculation to any decimal pt.
! Yudong: converted radians to degrees 

    EarthRadius = 6371  !km
    h           =  705  !km
    r = EarthRadius/(EarthRadius+h)
    Do np = 1, crtm_struc(n)%nprofiles 
      crtm_struc(n)%Geometry(np)%Sensor_Scan_Angle &
         = DASIN(r*DSIN(crtm_struc(n)%Geometry(np)%Sensor_Zenith_Angle*PI/180.0)) & 
             *180.0/PI 
    End Do


  end subroutine CRTM2_EMonly_geometry


  subroutine CRTM2_EMonly_run(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_Domain
    use LIS_RTMMod, only : LIS_sfcState
    use LIS_logMod, only : LIS_verify, LIS_getNextUnitNumber 
  USE CRTM_Module
  USE CRTM_SfcOptics
  USE CRTM_GeometryInfo
  USE CRTM_GeometryInfo_Define
  use LIS_histDataMod

    
    implicit none

    integer, intent(in) :: n 
    integer :: n_Sensors, n_Channels, n_Profiles
    integer :: iFOV
    real(fp)    :: Source_ZA
    integer :: AllocStatus
    TYPE(CRTM_GeometryInfo_type) :: GeometryInfo
    TYPE(CRTM_SfcOptics_type)    :: SfcOptics

    integer             :: status
    integer             :: t, ns, nc, ln, i, j, k 
    integer             :: SensorIndex, ChannelIndex, Index_Sat_Ang
    real                :: emiss_ave !average emissivity
    real,    pointer    :: wind_speed(:), land_coverage(:), land_temperature(:), snow_temperature(:), &
       soil_moisture_content(:), soil_temperature(:), canopy_water_content(:), &
       vegetation_fraction(:), snow_depth(:), snow_density(:), land_type(:), & 
       snow_coverage(:), leaf_area_index(:)
!    real, parameter     :: bgf=0.80 ! bare ground fraction of pixel
    real, parameter     :: bg_lai=0.005 ! lai of bare ground

    n_Sensors  = SIZE(crtm_struc(n)%ChannelInfo)
    n_Channels = SUM(crtm_struc(n)%ChannelInfo%n_Channels)
    n_Profiles = crtm_struc(n)%nprofiles

!   map surface properties to SFC
    call getsfcvar(LIS_sfcState(n), "Wind Speed", wind_speed)
    call getsfcvar(LIS_sfcState(n), "Land Coverage", land_coverage)
    call getsfcvar(LIS_sfcState(n), "Land Type", land_type)
    call getsfcvar(LIS_sfcState(n), "Land Temperature", land_temperature)
    call getsfcvar(LIS_sfcState(n), "Snow Temperature", snow_temperature)
    call getsfcvar(LIS_sfcState(n), "Soil Moisture Content",soil_moisture_content)
    call getsfcvar(LIS_sfcState(n), "Soil Temperature", soil_temperature)
    call getsfcvar(LIS_sfcState(n), "Canopy Water Content", canopy_water_content)
    call getsfcvar(LIS_sfcState(n), "Vegetation Fraction", vegetation_fraction)
    call getsfcvar(LIS_sfcState(n), "Snow Coverage", snow_coverage)
    call getsfcvar(LIS_sfcState(n), "Snow Depth", snow_depth)
    call getsfcvar(LIS_sfcState(n), "Snow Density", snow_density)
    call getsfcvar(LIS_sfcState(n), "Leaf Area Index", leaf_area_index)
    
    !---------------------------------------------
    ! Profile (tile) loop 
    !--------------------------------------------
    AllocStatus=CRTM_Allocate_SfcOptics( MAX_N_ANGLES, &  ! Input
         MAX_N_STOKES, &  ! Input
         SfcOptics     )  ! Output
    
    call LIS_verify(AllocStatus, 'Error in CRTM_Allocate_SfcOptics')
    ! write(*, *)"2--SfcOptics allocated"
    
    Profile_Loop: do t=1, n_Profiles 
       do k=1, 2  ! 1=basic; 2=bare ground
          crtm_struc(n)%SFC(t)%wind_speed            = wind_speed(t)
          crtm_struc(n)%SFC(t)%land_coverage         = land_coverage(t)
          crtm_struc(n)%SFC(t)%land_type             = CRTM_landmatch(nint(land_type(t)),&
               LIS_rc%lcscheme)
          crtm_struc(n)%SFC(t)%land_temperature      = land_temperature(t)
          crtm_struc(n)%SFC(t)%snow_temperature      = snow_temperature(t)

          
          i = LIS_domain(n)%tile(t)%col
          j = LIS_domain(n)%tile(t)%row 
          ! Soil moisture correction: rescaling method  
          If (sm_correction_struc(n)%c_type .eq. 1) then 
             crtm_struc(n)%SFC(t)%soil_moisture_content =   &
                  sm_correction_struc(n)%dst_mean(i, j) + &
                  (soil_moisture_content(t)-sm_correction_struc(n)%src_mean(i, j) ) * &
                  sm_correction_struc(n)%dst_sigma(i, j)/sm_correction_struc(n)%src_sigma(i, j)
             if ( crtm_struc(n)%SFC(t)%soil_moisture_content .LT. 0 ) & 
                  crtm_struc(n)%SFC(t)%soil_moisture_content = 0.0
          End if
          
          !      YDT 11/3/11: turn off canopy interception as it is too much relative 
          !        to CRTM's vegetation optical depth computation. This change is applied 
          !        to runs coded as EXP7xx. 
          !       crtm_struc(n)%SFC(t)%Canopy_Water_Content  = canopy_water_content(t)
          crtm_struc(n)%SFC(t)%Canopy_Water_Content  = 0.0 
          crtm_struc(n)%SFC(t)%soil_temperature      = soil_temperature(t)
          crtm_struc(n)%SFC(t)%vegetation_fraction   = vegetation_fraction(t)
          crtm_struc(n)%SFC(t)%snow_coverage         = snow_coverage(t)
          crtm_struc(n)%SFC(t)%snow_depth            = snow_depth(t)
          crtm_struc(n)%SFC(t)%snow_density          = snow_density(t)      
          
          if (k.eq.1) then
             crtm_struc(n)%bgf_v(t)      =  1.0*exp(-1.0*crtm_struc(n)%k_lai2vgf(t)*(leaf_area_index(t)/(1.0-crtm_struc(n)%bgf_fixed(t))))
             crtm_struc(n)%bgf_total(t)  =  crtm_struc(n)%bgf_fixed(t) + crtm_struc(n)%bgf_v(t)*(1.0-crtm_struc(n)%bgf_fixed(t))
             crtm_struc(n)%SFC(t)%lai       =  leaf_area_index(t)/(1.0-crtm_struc(n)%bgf_total(t))

             crtm_struc(n)%SFC(t)%soil_moisture_content = soil_moisture_content(t)
          else !k=2  bare ground
             crtm_struc(n)%SFC(t)%soil_moisture_content = soil_moisture_content(t)
             crtm_struc(n)%SFC(t)%lai                   = bg_lai
          end if
          
          ! Process geometry
          ! ...Compute derived geometry
          CALL CRTM_GeometryInfo_SetValue( GeometryInfo, &
               Geometry=crtm_struc(n)%Geometry(t),  &
               Sensor_Zenith_Angle = crtm_struc(n)%Geometry(t)%Sensor_Zenith_Angle) 
          CALL CRTM_GeometryInfo_Compute( GeometryInfo )
          ! ...Retrieve components into local variable
          CALL CRTM_GeometryInfo_GetValue( &
               GeometryInfo, &
               iFOV = iFOV, &
               Source_Zenith_Angle = Source_ZA )
          
          ! write(*, *)"1--Geometry processed, Sensor_Zenith_Angle=", & 
          !              GeometryInfo%Sensor_Zenith_Radian*180.0/3.1415926
          
          ! Average surface skin temperature for multi-surface types
          
          CALL CRTM_Compute_SurfaceT( crtm_struc(n)%sfc(t), SfcOptics )
          
          ! write(*, *)"3--SurfaceT computed "
          
          ! -----------
          ! Sensor loop
          ! -----------
          ! Initialise channel counter for channel(l)/sensor(n) count
          ln = 0
          
          Sensor_Loop: DO ns = 1, n_Sensors
             
             ! Shorter name
             SensorIndex = crtm_struc(n)%ChannelInfo(ns)%Sensor_Index
             
             ! ------------
             ! Channel loop
             ! ------------
             Channel_Loop: DO nc = 1, n_Channels 
                
                ! Shorter name
                ChannelIndex = crtm_struc(n)%ChannelInfo(ns)%Channel_Index(nc)
                
                ! Increment channel counter
                ln = ln + 1
                
                !YDT:  needs careful exam later. The following settings can reproduce
                ! the emissivity values from CRTM2_Forward runs (and Ken's CRTM v1 runs)
                ! But may not be robust. Need to implement other checks as in CRTM_RTSolution
                ! Indicate SfcOptics ARE to be computed
                SfcOptics%Compute_Switch = SET
                SfcOptics%Index_Sat_Ang = 1        
                SfcOptics%Angle(1) = crtm_struc(n)%Geometry(t)%Sensor_Zenith_Angle
                SfcOptics%Weight(1) = 1.0 
                
                ! ---------------------------------------------------------------------
                ! Assign the number of Stokes parameters. Currently, this is set to 1
                ! for decoupled polarization between surface and atmosphere.
                ! Remember for polarised microwave instruments, each frequency's Stokes
                ! parameter is treated as a separate channel.
                ! ---------------------------------------------------------------------
                SfcOptics%n_Stokes = 1
                SfcOptics%n_Angles = 1
                
                !write(*, *)"4--Before Compute_SfcOptics, Index_Sat_Ang= ", SfcOptics%Index_Sat_Ang
                status = CRTM_Compute_SfcOptics(      &
                     crtm_struc(n)%Sfc(t),                   &
                     GeometryInfo,          &
                     SensorIndex,           &
                     ChannelIndex,           &
                     SfcOptics,             & 
                     crtm_struc(n)%SOV(t)) 
                
                call LIS_verify(status, 'Error in CRTM2_Compute_SfcOptics Model')
                
                Index_Sat_Ang = SfcOptics%Index_Sat_Ang
                !write(*, *)"5--After Compute_SfcOptics, Index_Sat_Ang= ", Index_Sat_Ang, &
                !         "ln=", ln
                
                crtm_struc(n)%SfcOptics(ln, t)%Emissivity( Index_Sat_Ang, 1 ) = & 
                     SfcOptics%Emissivity( Index_Sat_Ang, 1 )
                if (k.eq.1) then
                   crtm_struc(n)%emissivity(ln, t)=real(crtm_struc(n)%SfcOptics(ln, t)%Emissivity( Index_Sat_Ang, 1 ))
                else !k=2 assign to bare ground emissivity
                   crtm_struc(n)%emissivity_bg(ln, t)=real(crtm_struc(n)%SfcOptics(ln, t)%Emissivity( Index_Sat_Ang, 1 ))
                endif
!                if (ln.eq.2) then 
!                   print *, crtm_struc(n)%SfcOptics(ln, t)%Emissivity( Index_Sat_Ang, 1 )
!                end if
!                call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_EMISSIVITY, value=      &
!                     real(crtm_struc(n)%SfcOptics(ln, t)%Emissivity(Index_Sat_Ang, 1 )), & 
!                     vlevel=ln, unit="-",direction="-")
!                
                !   call LIS_diagnoseOutputVar(n, t, LIS_MOC_RTM_MPDI, value=       &
                !         real(crtm_struc(n)%SfcOptics(ln, t)%Emissivity(Index_Sat_Ang, 1 )), & 
                !         vlevel=ln, unit="-",direction="-")
                
                !write(*, *)"6--After saving emissivity "
             END DO Channel_Loop
          END DO Sensor_Loop
       end do
          
          !    call LIS_diagnoseOutputVar(n, t, LIS_MOC_RTM_SM, value=       &
          !         real(crtm_struc(n)%SFC(t)%soil_moisture_content),             &
          !         vlevel=1, unit="m3/m3",direction="-")
          
          
!          call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_SM, value=       &
!               real(crtm_struc(n)%SFC(t)%soil_moisture_content),             &
!               vlevel=1, unit="m3/m3",direction="-")
          
          
          ! ---------------------------------------------------
          ! Deallocate local sensor independent data structures
          ! ---------------------------------------------------
          
             ! Initialise channel counter for channel(l)/sensor(n) count
          ln = 0
          Nuther_Sensor_Loop: DO ns = 1, n_Sensors
             Nuther_Channel_Loop: DO nc = 1, n_Channels 
                ! Increment channel counter
                ln = ln + 1
             !   print *,  crtm_struc(n)%emissivity(ln, t), crtm_struc(n)%emissivity_bg(ln, t)
                emiss_ave=(1.0-crtm_struc(n)%bgf_total(t))*crtm_struc(n)%emissivity(ln, t)+(crtm_struc(n)%bgf_total(t))*crtm_struc(n)%emissivity_bg(ln, t)
                crtm_struc(n)%emissivity_ave(ln, t)=emiss_ave
!                emiss_ave=crtm_struc(n)%emissivity(ln, t)
                
                call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_EMISSIVITY, value=      &
                     emiss_ave, & 
                     vlevel=ln, unit="-",direction="-")
                
                !if (ln.eq.2) then 
                !   print *, crtm_struc(n)%SfcOptics(ln, t)%Emissivity( Index_Sat_Ang, 1 )
                !end if
                
                END DO Nuther_Channel_Loop
             END DO Nuther_Sensor_Loop
             call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_SM, value=       &
                  real(crtm_struc(n)%SFC(t)%soil_moisture_content),             &
                  vlevel=1, unit="m3/m3",direction="-")
       END DO Profile_Loop
       
       AllocStatus = CRTM_Destroy_SfcOptics( SfcOptics )
       call LIS_verify(AllocStatus, 'Error in CRTM_Destroy_SfcOptics')
       
       !write(*, *)"7--After CRTM_Destroy_SfcOptics "
       
       !write(*, *)"8--End of Forward run "
       
     end subroutine CRTM2_EMonly_run
        
!BOP
! !ROUTINE: CRTM2_EMonly_output
! \label{CRTM2_EMonly_output}
!
! !INTERFACE: 
  subroutine CRTM2_EMonly_output(n)
! !USES: 
!!!!!$    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_masterproc
!!!!!$    use LIS_logMod, only : LIS_verify, LIS_getNextUnitNumber, &
!!!!!$         LIS_releaseUnitNumber
!!!!!$    use LIS_fileIOMod,  only : LIS_create_output_directory
!!!!!$    use LIS_historyMod, only : LIS_writevar_gridded

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine writes the CRTM solution fields to a time stamped output
!  file
!EOP
    
!!!!!$    integer           :: ftn 
!!!!!$    integer           :: c,r,k,gid, n_channels, Index_Sat_Ang
!!!!!$    character*100     :: crtm_filename
!!!!!$    real              :: datafield(LIS_rc%ngrid(n))

!!!!!$     if(LIS_masterproc) then 
!!!!!$        ftn = LIS_getNextUnitNumber()
!!!!!$        call create_crtm_output_filename(crtm_filename)

!!!!!$        call LIS_create_output_directory('CRTM',style=1)
!!!!!$        open(ftn,file=trim(crtm_filename), form='unformatted')        
!!!!!$     endif
     
!!!!!$     n_Channels = SUM(crtm_struc(n)%ChannelInfo%n_Channels)
     
! svk
! Currently we are assuming that tiling is not used. If tiling is used, the 
! data needs to be averaged upto the grid scale first. 

!!!!!$     do k=1,n_Channels
!!!!!$        do r=1,LIS_rc%lnr(n)
!!!!!$           do c=1,LIS_rc%lnc(n)
!!!!!$              gid = LIS_domain(n)%gindex(c,r)
!!!!!$              if(gid.ne.-1) then 
!!!!!$                 Index_Sat_Ang = crtm_struc(n)%SfcOptics(k, gid)%Index_Sat_Ang
!!!!!$                 datafield(gid) =  &
!!!!!$                     crtm_struc(n)%SfcOptics(k, gid)%Emissivity(Index_Sat_Ang, 1)
!!!!!$              endif
!!!!!$           enddo
!!!!!$        enddo
        
!!!!!$        call LIS_writevar_gridded(ftn,n,datafield)
!!!!!$     enddo

!!!!!$     if(LIS_masterproc) then 
!!!!!$        call LIS_releaseUnitNumber(ftn)
!!!!!$     endif


  end subroutine CRTM2_EMonly_output



  subroutine create_crtm_output_filename(obsname)
! !USES: 
    use LIS_coreMod, only : LIS_rc
    
! !ARGUMENTS: 
    character(len=*)      :: obsname
! 
! !DESCRIPTION: 
! 
!EOP

    character(len=12) :: cdate1
    
    write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
         LIS_rc%yr, LIS_rc%mo, &
         LIS_rc%da, LIS_rc%hr,LIS_rc%mn
    
!!$    obsname = trim(LIS_rc%odir)//'/EXP'//trim(adjustl(LIS_rc%expcode))//&
!!$         '/CRTM/'//cdate1(1:4)// '/' // cdate1(1:8) // '/' //cdate1//'.gs4r'
! for consistency with other rtm output writing, keep same as the others for now.
      obsname = trim(LIS_rc%odir)//'/'//&
         'CRTM/'//cdate1//'.1gs4r'
end subroutine create_crtm_output_filename


  subroutine getsfcvar(sfcState, varname, var)
! !USES: 
    use LIS_logMod,  only : LIS_verify
    
    implicit none
    
    type(ESMF_State)      :: sfcState
    type(ESMF_Field)      :: varField
    character(len=*)      :: varname
    real, pointer         :: var(:)
    integer               :: status

    call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
    call LIS_verify(status, 'Error in StateGet: CRTM_handlerMod '//trim(varname))
    call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
    call LIS_verify(status, 'Error in FieldGet: CRTM_handlerMod '//trim(varname))

  end subroutine getsfcvar


  integer function CRTM_landmatch(i, classification)

  !USES:      

    use CRTM_Module

    implicit none
  ! !ARGUMENTS: 
    integer, intent(in) :: i               ! index into CRTM lookup table
    character(len=*), intent(in) :: classification  ! classification scheme, e.g., UMD
  !                                          (1), USGS (2), UMD (3), IGBP (4)
  ! Description
  ! Output 
    integer :: j                           ! temp variable to store returned value

  ! UMD Land Cover Classification from:  http://www.geog.umd.edu/landcover/1km-map/meta-data.html
  !kwh: ARBITRARY VALUES CURRENTLY ASSIGNED FOR TESTING.

    integer, dimension (0:13), parameter :: umd_crtm_match = &
      (/  &
           INVALID_LAND              , &        ! 0  - WATER (and Goode's interrupted space)   
           PINE_FOREST               , &        ! 1  - EVERGREEN NEEDLELEAF FOREST             
           BROADLEAF_FOREST          , &        ! 2  - EVERGREEN BROADLEAF FOREST              
           PINE_FOREST               , &        ! 3  - DECIDUOUS NEEDLELEAF FOREST             
           BROADLEAF_FOREST          , &        ! 4  - DECIDUOUS BROADLEAF FOREST              
           BROADLEAF_PINE_FOREST     , &        ! 5  - MIXED FOREST                            
           BROADLEAF_BRUSH           , &        ! 6  - WOODLAND                                
           BROADLEAF_BRUSH           , &        ! 7  - WOODED GRASSLAND                        
           SCRUB                     , &        ! 8  - CLOSED SHRUBLAND                        
           SCRUB_SOIL                , &        ! 9  - OPEN SHRUBLAND                          
           GRASS_SCRUB               , &        ! 10 - GRASSLAND                               
           TILLED_SOIL               , &        ! 11 - CROPLAND                                
           COMPACTED_SOIL            , &        ! 12 - BARE GROUND                             
           URBAN_CONCRETE            &          ! 13 - URBAN AND BUILT-UP                      
       /)                                                         
               
  !                
  !USGS Land Use/Land Cover System Legend (Modified Level 2) (as taken from USGSEDC)
  !kwh: NEED TO PUT IN CRTM MATCH HERE, AS WITH GFS.  CURRENLT ARBITRARY VALUES ASSIGNED

    integer, dimension (24), parameter :: usgs_crtm_match = &
    (/  &                                         !Value 	Code 	Description
       URBAN_CONCRETE            , &          	! 1  	100 	Urban and Built-Up Land
       TILLED_SOIL               , &          	! 2  	211 	Dryland Cropland and Pasture
       TILLED_SOIL               , &          	! 3  	212 	Irrigated Cropland and Pasture
       TILLED_SOIL               , &          	! 4  	213 	Mixed Dryland/Irrigated Cropland and Pasture
       TILLED_SOIL               , &          	! 5  	280 	Cropland/Grassland Mosaic
       TILLED_SOIL               , &          	! 6  	290 	Cropland/Woodland Mosaic
       MEADOW_GRASS              , &          	! 7  	311 	Grassland
       SCRUB                     , &          	! 8  	321 	Shrubland
       GRASS_SCRUB               , &          	! 9  	330 	Mixed Shrubland/Grassland
       BROADLEAF_BRUSH           , &          	! 10 	332 	Savanna
       BROADLEAF_FOREST          , &          	! 11 	411 	Deciduous Broadleaf Forest
       PINE_FOREST               , &          	! 12 	412 	Deciduous Needleleaf Forest
       BROADLEAF_FOREST          , &          	! 13 	421 	Evergreen Broadleaf Forest
       PINE_FOREST               , &         	! 14 	422 	Evergreen Needleleaf Forest
       BROADLEAF_PINE_FOREST     , &         	! 15 	430 	Mixed Forest
       INVALID_LAND              , &		! 16 	500 	Water Bodies
       BROADLEAF_BRUSH           , &		! 17 	620 	Herbaceous Wetland
       BROADLEAF_BRUSH           , &		! 18 	610 	Wooded Wetland
       COMPACTED_SOIL 	       , &		! 19 	770 	Barren or Sparsely Vegetated
       TUNDRA		       , &		! 20 	820 	Herbaceous Tundra
       TUNDRA  		       , &		! 21 	810 	Wooded Tundra
       TUNDRA  		       , &		! 22 	850 	Mixed Tundra
       TUNDRA  		       , &		! 23 	830 	Bare Ground Tundra
       INVALID_LAND   	       &		! 24 	900 	Snow or Ice  
     /)  
                                  
  !  integer, parameter :: USGS_  =99 	  	! Interrupted Areas (Goodes Homolosine Projection)
  !  integer, parameter :: USGS_  =100 	Missing Data
  !
  !   GFS/GDAS vegtype taken from Weizhong Zeng driver

    integer, dimension (13), parameter :: gfs_crtm_match = &
    (/  &                              
         BROADLEAF_FOREST         , &   ! 1  - BROADLEAF-EVERGREEN TREES  (TROPICAL FOREST)
         BROADLEAF_FOREST         , &   ! 2  - BROADLEAF-DECIDUOUS TREES
         BROADLEAF_PINE_FOREST    , &   ! 3  - BROADLEAF AND NEEDLELEAF TREES (MIXED FOREST)
         PINE_FOREST              , &   ! 4  - NEEDLELEAF-EVERGREEN TREES
         PINE_FOREST              , &   ! 5  - NEEDLELEAF-DECIDUOUS TREES (LARCH)
         BROADLEAF_BRUSH          , &   ! 6  - BROADLEAF TREES WITH GROUNDCOVER (SAVANNA)
         SCRUB                    , &   ! 7  - GROUNDCOVER ONLY (PERENNIAL)  
         SCRUB                    , &   ! 8  - BROADLEAF SHRUBS WITH PERENNIAL GROUNDCOVER
         SCRUB_SOIL               , &   ! 9  - BROADLEAF SHRUBS WITH BARE SOIL
         TUNDRA                   , &   ! 10 - DWARF TREES AND SHRUBS WITH GROUNDCOVER (TUNDRA)
         COMPACTED_SOIL           , &   ! 11 - BARE SOIL
         TILLED_SOIL              , &   ! 12 - CULTIVATIONS (THE SAME PARAMETERS AS FOR TYPE 7)
         COMPACTED_SOIL           &     ! 13 - GLACIAL (THE SAME PARAMETERS AS FOR TYPE 11) 
    /)                                


  ! IGBP Land Cover Categories as taken from USGS EDC
  !kwh: NEED TO PUT IN CRTM MATCH HERE, AS WITH GFS.  CURRENTLY ARBITRARY VALUES ASSIGNED
    integer, dimension (17), parameter :: igbp_crtm_match = &
    (/  &                              
         PINE_FOREST                , &               ! 1  Evergreen Needleleaf Forest
         BROADLEAF_FOREST           , &               ! 2  Evergreen Broadleaf Forest
         PINE_FOREST                , &               ! 3  Deciduous Needleleaf Forest
         BROADLEAF_FOREST           , &               ! 4  Deciduous Broadleaf Forest
         BROADLEAF_PINE_FOREST      , &               ! 5  Mixed Forest
         SCRUB                      , &               ! 6  Closed Shrublands
         SCRUB_SOIL                 , &               ! 7  Open Shrublands
         BROADLEAF_BRUSH            , &               ! 8  Woody Savannas
         BROADLEAF_BRUSH            , &               ! 9  Savannas
         MEADOW_GRASS               , &               ! 10 Grasslands
         BROADLEAF_BRUSH            , &               ! 11 Permanent Wetlands
         TILLED_SOIL                , &               ! 12 Croplands
         URBAN_CONCRETE             , &               ! 13 Urban and Built-Up
         TILLED_SOIL                , &               ! 14 Cropland/Natural Vegetation Mosaic
         INVALID_LAND               , &               ! 15 Snow and Ice
         COMPACTED_SOIL             , &               ! 16 Barren or Sparsely Vegetated
         INVALID_LAND               &                 ! 17 Water Bodies                                     
  !                                                   ! 99 Interrupted Areas (Goodes Homolosine Projection)
  !                                                   !100 Missing Data           
    /)                                
                                    
  ! 
  ! !DESCRIPTION
  ! scheme:
  ! #1 use the UMD landcover
  ! #2 use the USGS landcover data
  ! #3 use the GFS landcover data
  ! #4 use the IGBP landcover data
  !
  ! This subroutine matches the land (veg type) classification
  ! to that of CRTM

  !EOP

  if     (classification .eq. "UMD") then !UMD
    !need to bounce matchings off someone; similar matchings to weizhong driver
    j = umd_crtm_match(i)
  elseif (classification .eq. "USGS") then !USGS
    !insert error code for 'not yet implemented'
    !j = usgs_crtm_match(i)
  elseif (classification .eq. "GFS") then !GFS
    j = gfs_crtm_match(i)
  elseif (classification .eq. "MODIS") then !IGBP
    !insert error code for 'not yet implemented'
    !j = igbp_crtm_match(i)
  else
    !kwh: insert error code
  end if
  CRTM_landmatch = j
  return
end function  CRTM_landmatch

#endif
end module CRTM2_EMMod



