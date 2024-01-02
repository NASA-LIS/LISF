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
module CRTM2_handlerMod
!BOP
!
! !MODULE: CRTM2_handlerMod
!
! !DESCRIPTION:
!    This module provides the routines to control the execution of 
!    CRTM 2.x from within LIS. The subroutines provide mapping of the surface
!    and atmospheric profiles to CRTM2, which then runs the forward 
!    model to generate the radiative quantities. 
!
! !REVISION HISTORY:
! 14 Sep 2010: Yudong Tian; Modifed from CRTM_handlerMod.F90 to support CRTM2. 
!
! !USES:        

#if (defined RTMS)

  USE CRTM_Module
  use ESMF

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: CRTM2_Forward_initialize
  public :: CRTM2_Forward_f2t
  public :: CRTM2_Forward_geometry
  public :: CRTM2_Forward_run
  public :: CRTM2_Forward_output
  public :: CRTM_landmatch
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: crtm_struc
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
     type(CRTM_ChannelInfo_type),  allocatable :: ChannelInfo(:)
     type(CRTM_Geometry_type), allocatable :: Geometry(:)
! forward variables
     type(CRTM_Atmosphere_type),   allocatable :: Atm(:)
     type(CRTM_Surface_type),      allocatable :: Sfc(:)
     type(CRTM_RTSolution_type),   allocatable :: RTSolution(:,:)
  end type crtm_type_dec

  type(crtm_type_dec), allocatable :: crtm_struc(:)
  SAVE

contains
!BOP
! 
! !ROUTINE: CRTM2_Forward_initialize
! \label{CRTM2_Forward_initialize}
! 
! !INTERFACE:
  subroutine CRTM2_Forward_initialize()
! !USES:
!    use CRTM_Module
    use LIS_coreMod,    only : LIS_rc
    use LIS_logMod,     only : LIS_logunit, LIS_verify
    use LIS_RTMMod,     only : LIS_sfcState

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for noah-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noah from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readCRTM2crd](\ref{readCRTM2crd}) \newline
!    reads the runtime options for Noah LSM
!  \end{description}
!EOP
    implicit none

    integer :: n 
    integer :: m 
    integer :: n_Channels, n_profiles
    integer :: status

    allocate(crtm_struc(LIS_rc%nnest))

    call readCRTM2crd()

    do n=1,LIS_rc%nnest
       crtm_struc(n)%nprofiles = LIS_rc%ntiles(n) ! all grid points in LIS

       allocate(crtm_struc(n)%ChannelInfo(crtm_struc(n)%nsensors))
       allocate(crtm_struc(n)%Geometry(crtm_struc(n)%nprofiles))
       allocate(crtm_struc(n)%Atm(crtm_struc(n)%nprofiles))
       allocate(crtm_struc(n)%Sfc(crtm_struc(n)%nprofiles))
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
       n_channels = SUM(crtm_struc(n)%ChannelInfo%n_Channels)
       n_Profiles = crtm_struc(n)%nprofiles
!---------------------------------------------------------------------
! Allocate structure arrays
!---------------------------------------------------------------------
 
       allocate(crtm_struc(n)%RTSolution(n_Channels, n_Profiles))

!---------------------------------------------------------------------
! Allocate FORWARD structure
!---------------------------------------------------------------------
       call CRTM_Atmosphere_Create(crtm_struc(n)%Atm, & 
            crtm_struc(n)%nlayers, &
            crtm_struc(n)%nabsorbers, crtm_struc(n)%nclouds, &
            crtm_struc(n)%naerosols) 
       if ( ANY(.NOT. CRTM_Atmosphere_Associated(crtm_struc(n)%Atm) ) ) Then
          Write(*, *) 'Error allocating CRTM Atmosphere structure'
          Stop 
       End If

       call CRTM_RTSolution_Create(crtm_struc(n)%RTSolution, & 
            crtm_struc(n)%nlayers) 
       if ( ANY(.NOT. CRTM_RTSolution_Associated(crtm_struc(n)%RTSolution) ) ) Then
          Write(*, *) 'Error allocating CRTM RTSolution structure'
          Stop
       End If

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

    enddo


  end subroutine CRTM2_Forward_initialize 


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

  subroutine CRTM2_Forward_f2t(n)

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
    real,   pointer        :: lpress(:), press(:), tmp(:), q2(:),o3(:)
    real                   :: qsmall, o3small 
    real                   :: qmax
    
    qsmall = 1E-6
    o3small = 1E-10
    do t=1,LIS_rc%ntiles(n)
! These are currently hardcoded. Modify later
       crtm_struc(n)%Atm(t)%Climatology = TROPICAL
       crtm_struc(n)%Atm(t)%Absorber_Id = (/H2O_ID, O3_ID/)
       crtm_struc(n)%Atm(t)%Absorber_Units = &
            (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    enddo
!assigning forcing****  k=1 reserved for surface value  ****

    do k=2,crtm_struc(n)%nlayers + 2

! The order (top to bottom) might need to be reordered
!    Level Pressure       
       call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_LPressure%varname(k)), &
            lpressField, rc=status)
       call LIS_verify(status,'ESMF_StateGet with LPressure in CRTM2_Forward_f2t')
       
       call ESMF_FieldGet(lpressField, localDE=0,farrayPtr=lpress,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet with LPressure in CRTM2_Forward_f2t')
       kk = (crtm_struc(n)%nlayers+2)-k
       do t=1,LIS_rc%ntiles(n)
          crtm_struc(n)%Atm(t)%Level_Pressure(kk) = lpress(t)/100.0
       enddo
    enddo

    do k=2,crtm_struc(n)%nlayers+1
! Pressure
       call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Psurf%varname(k)), &
            pressField, rc=status)
       call LIS_verify(status,'ESMF_StateGet with Pressure in CRTM2_Forward_f2t')
       
       call ESMF_FieldGet(pressField, localDE=0,farrayPtr=press,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet with Pressure in CRTM2_Forward_f2t')
       kk = crtm_struc(n)%nlayers+2-k
       do t=1,LIS_rc%ntiles(n)
          crtm_struc(n)%Atm(t)%Pressure(kk) = press(t)/100.0
       enddo
! Temperature
       call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(k)), &
            tmpField, rc=status)
       call LIS_verify(status,'ESMF_StateGet with Temperature in CRTM2_Forward_f2t')
       
       call ESMF_FieldGet(tmpField, localDE=0,farrayPtr=tmp,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet with Temperature in CRTM2_Forward_f2t')
       kk = crtm_struc(n)%nlayers+2-k
       do t=1,LIS_rc%ntiles(n)
          crtm_struc(n)%Atm(t)%Temperature(kk) = tmp(t)
       enddo
! Absorber 1: H2O
       call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Qair%varname(k)), &
            q2Field, rc=status)
       call LIS_verify(status,'ESMF_StateGet with Q2 in CRTM2_Forward_f2t')
       
       call ESMF_FieldGet(q2Field, localDE=0,farrayPtr=q2,rc=status)
       call LIS_verify(status, 'ESMF_FieldGet with Q2 in CRTM2_Forward_f2t')

       kk = crtm_struc(n)%nlayers+2-k
       do t=1,LIS_rc%ntiles(n)
          !LIS unit is Specific Humidity (SH) in kg h20/kg air with water
          !CRTM unit is mixing ratio (MR) in units of  g h20/kg dry air
          !CRTM Profile Utility has conversion that assumes g/kg SH
!          crtm_struc(n)%Atm(t)%Absorber(kk,1) = max(q2(t),qsmall)*1000.0_fp
          qmax = max(q2(t)*1000.0,qsmall)
          crtm_struc(n)%Atm(t)%Absorber(kk,1) = SA_to_MR(qmax*1.0_fp)
       enddo 

! Absorber 2: O3
       if(LIS_FORC_O3%selectOpt.eq.1) then 
          call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_O3%varname(k)), &
               o3Field, rc=status)
          call LIS_verify(status,'ESMF_StateGet with O3 in CRTM2_Forward_f2t')
          
          call ESMF_FieldGet(o3Field, localDE=0,farrayPtr=o3,rc=status)
          call LIS_verify(status, 'ESMF_FieldGet with O3 in CRTM2_Forward_f2t')
       endif
       kk = crtm_struc(n)%nlayers+2-k
       do t=1,LIS_rc%ntiles(n)
          if(LIS_FORC_O3%selectOpt.eq.1) then 
             crtm_struc(n)%Atm(t)%Absorber(kk,2) = &
                  MR_to_PPMV(max(o3(t),o3small)*1000.0_fp)
          else
!hack
             crtm_struc(n)%Atm(t)%Absorber(kk,2) = 0.0
          endif
       enddo
    enddo

  end subroutine CRTM2_Forward_f2t


  subroutine CRTM2_Forward_geometry(n)

    implicit none

    integer,   intent(in)    :: n 
    real, parameter :: PI = 3.14159265
    real                     :: r  !distance ratio (see CRTM User Guide; Geometry Structure)
    real                     :: h  !altititude of satellite
    real                     :: EarthRadius !can vary


  crtm_struc(n)%Geometry%Sensor_Zenith_Angle = crtm_struc(n)%zenith_angle

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
    crtm_struc(n)%Geometry%Sensor_Scan_Angle &
         = DASIN(r*DSIN(crtm_struc(n)%Geometry%Sensor_Zenith_Angle*PI/180.0)) &
           *180.0/PI 


  end subroutine CRTM2_Forward_geometry


  subroutine CRTM2_Forward_run(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc
    use LIS_RTMMod, only : LIS_sfcState
    use LIS_logMod, only : LIS_verify
    
    implicit none

    integer, intent(in) :: n 
    integer             :: status
    integer             :: t
    real,    pointer    :: wind_speed(:), land_coverage(:), land_temperature(:), snow_temperature(:),&
       soil_moisture_content(:), soil_temperature(:), canopy_water_content(:), &
       vegetation_fraction(:), snow_depth(:), snow_density(:), land_type(:), snow_coverage(:)

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
    
    do t=1, LIS_rc%ntiles(n)
       crtm_struc(n)%SFC(t)%wind_speed            = wind_speed(t)
       crtm_struc(n)%SFC(t)%land_coverage         = land_coverage(t)
       crtm_struc(n)%SFC(t)%land_type             = CRTM_landmatch(nint(land_type(t)),&
            LIS_rc%lcscheme)
       crtm_struc(n)%SFC(t)%land_temperature      = land_temperature(t)
       crtm_struc(n)%SFC(t)%snow_temperature      = snow_temperature(t)
       crtm_struc(n)%SFC(t)%soil_moisture_content = soil_moisture_content(t)
       crtm_struc(n)%SFC(t)%soil_temperature      = soil_temperature(t)
       crtm_struc(n)%SFC(t)%vegetation_fraction   = vegetation_fraction(t)
       crtm_struc(n)%SFC(t)%snow_coverage         = snow_coverage(t)
       crtm_struc(n)%SFC(t)%snow_depth            = snow_depth(t)
       crtm_struc(n)%SFC(t)%snow_density          = snow_density(t)      
    enddo

! 9/20/2010 Yudong: trying CRTM_Forward() instead CRTM_K_Matrix()
! for computational efficiency. No other parts
! were changed. It turns out it runs twice faster now, with identical 
! emissivity and Tb results.  

    status = CRTM_Forward(crtm_struc(n)%Atm, &
         crtm_struc(n)%Sfc,                   &
         crtm_struc(n)%Geometry,          &
         crtm_struc(n)%ChannelInfo,           &
         crtm_struc(n)%RTSolution)

    call LIS_verify(status, 'Error in CRTM2 Forward Model')


  end subroutine CRTM2_Forward_run

!BOP
! !ROUTINE: CRTM2_Forward_output
! \label{CRTM2_Forward_output}
!
! !INTERFACE: 
  subroutine CRTM2_Forward_output(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_masterproc
    use LIS_logMod, only : LIS_verify, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_historyMod, only : LIS_writevar_gridded
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine writes the CRTM solution fields to a time stamped output
!  file
!EOP
    
    integer           :: ftn 
    integer           :: c,r,k,gid, n_channels
    character(len=LIS_CONST_PATH_LEN) :: crtm_filename
    real              :: datafield(LIS_rc%ngrid(n))

     if(LIS_masterproc) then 
        ftn = LIS_getNextUnitNumber()
        call create_crtm_output_filename(crtm_filename)

        call LIS_create_output_directory('CRTM')
        open(ftn,file=trim(crtm_filename), form='unformatted')        
     endif
     
     n_channels = SUM(crtm_struc(n)%ChannelInfo%n_Channels)
     
! svk
! Currently we are assuming that tiling is not used. If tiling is used, the 
! data needs to be averaged upto the grid scale first. 

     do k=1,n_Channels
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              gid = LIS_domain(n)%gindex(c,r)
              if(gid.ne.-1) then 
                 datafield(gid) = crtm_struc(n)%RTSolution(k,gid)%Brightness_Temperature
              endif
           enddo
        enddo
        
        call LIS_writevar_gridded(ftn,n,datafield)
     enddo

     do k=1,n_Channels
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              gid = LIS_domain(n)%gindex(c,r)
              if(gid.ne.-1) then 
                 datafield(gid) = crtm_struc(n)%RTSolution(k,gid)%Surface_Emissivity
              endif
           enddo
        enddo
        
        call LIS_writevar_gridded(ftn,n,datafield)
     enddo

     if(LIS_masterproc) then 
        call LIS_releaseUnitNumber(ftn)
     endif


  end subroutine CRTM2_Forward_output



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
end module CRTM2_handlerMod



