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
! !ROUTINE: mos_main
!  \label{mos_main}
!
! !REVISION HISTORY:
! 11  Apr 2000: Brian Cosgrove; Initial Code
! 12  May 2000: Brian Cosgrove/Jared Entin; Added code to prevent 
!               numerical instablility due to small values of wind
!               and humidity
! 22 Aug. 2000: Brian Cosgrove; Modified code for output of
!               standard LDAS output variables.  Added LDAS and MOS
!               modules into call for MOSTILE and also added many
!               new variables in the MOS%RETURN section 
! 25 Aug. 2000: Brian Cosgrove; Fixed snowcover fraction output variable
! 21 Sep. 2000: Brian Cosgrove; Fixed code so that rain output field
!               is set to zero when snow output field is greater than zero
! 27 Sep. 2000: Brian Cosgrove; Fixed code so that 1 meter soil moisture
!               is correctly calculated (FACTOR variable had been
!               incorrectly computed before this fix)
! 23 Mar. 2001: Jon Radakovich; Updated for PSAS temperature assimilation
! 19 Apr. 2001  Updated scheme for calculating albedo, using sibalb module
! 05 Sep. 2001: Brian Cosgrove; Added in volumetric variables soilwm,soilww
!               for use in LDAS output.  Zthick changed from 50 to 10 sometime
!               by someone else.  Added volumetric output at 4 more levels
!               corresponding to OK mesonet levels.  Changed calculation
!               of soil wetness output variables
! 05 Feb. 2002: Brian Cosgrove; Changed Zthick back to 50 from 10 after talk
!               with Jon Radokovich.  10 is correct for observation height, but 
!               was causing fluxes that were too high, so changed back to
!               'incorrect' value of 50
! 22 Apr. 2002: Urszula Jambor; Added conditional to suppress calculation
!               of dewpt temp. if using Aaron Berg's reanalysis ECMWF data.
! 13 Sep. 2002: Urszula Jambor; Reversed suppression of dew point temp.
!               calculation if using Aaron Berg's reanalysis ECMWF data.
! 27 Nov. 2002: Urszula Jambor; Restricted mos%snow to max. of 100000.0
!               to prevent problems in GRIB output snow fields.
! 12 Dec. 2002: Brian Cosgrove; Fixed usage of Wiltpoint variable.  Before,
!               Wiltpoint1 and Wiltpoint2 were used in calculation of 
!               root zone soil moisture availability...now, only Wiltpoint2
!               is used since wiltpoint1 is not the correct wilting point
!               needed for the calculation.
!  25 Sep 2007: Sujay Kumar; Upgraded for LIS5.0
! 
! !INTERFACE: 
subroutine mos_main(n)
! !USES:   
  use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
  use LIS_constantsMod,  only : LIS_CONST_RHOFW
  use LIS_timeMgrMod, only : LIS_isAlarmRinging
  use LIS_logMod, only : LIS_logunit
  use LIS_albedoMod, only : LIS_alb
  use LIS_histDataMod
  use mos_lsmMod     
  use sibalb_module  

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
!
! !DESCRIPTION:
!  This is the entry point for calling the Mosaic LSM physics. This routine
!  calls the {\tt MOSTILE} routine that performs the land surface
!  computations, to solve for water and energy equations. For documentation
!  of the {\tt MOSTILE} and other Mosaic routines, please see the reference: 
!  Koster, R. D., and M. J. Suarez, Energy and Water Balance 
!  Calculations in the MOSAIC LSM. NASA Technical Memorandum 
!  104606, 9, 76 pp., 1996.
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

  integer index1
  INTEGER FLAG            !Flag for top 1 meter soil mst derivation
  !1=have value, 0=don't have value
  INTEGER SEC           !Current number of seconds into day
  INTEGER RESET         !Flag for TILE SUBROUTINE
  !1=convert ETA volumetric to %saturation
  !  soil moisture because MOSAIC soil
  !  moisture is being reset to ETA values
  !  and unit conversion needs to occur
  REAL DEPTH         !Depth value used in top 1meter soil mst calcs
  REAL FACTOR        !Factor by which to multiply soil moisture
  !in a layer to arrive at its contribution
  !to soil mst in top 1 meter
  REAL TEMP1,TEMP2        !Factors used in top 1 meter soil mst calcs
  REAL ALBEDO1,ALBEDO2
  REAL ALBEDO3,ALBEDO4
  REAL M_ESAT             !Function which uses DPT2M to
  !  compute EM.
  REAL CTT                ! 
  REAL M_QSAT           !Function 
  REAL SNWMID           !
  REAL ALHE
  REAL ALHS
  REAL ALHM
  REAL ESATTC
  REAL GRAV,RGAS,CP     !Constants
  REAL VKRMN,EPSI       !Constants
  REAL B,C,D,ZTHICK     !Constants
  REAL RHOAIR           !Density of air
  REAL QM               !Parameter from EM,EPSI,PSUR
  REAL QC                 !Parameter from EA,EPSI,PSUR
  REAL DTV                !
  REAL RI                 !
  REAL DRIDTC             !
  REAL DRIDQC             !
  REAL DUMMY              !
  REAL CN                 !
  REAL T                  !
  REAL FU                 !
  REAL FT                 !
  REAL DFTDRI             !
  REAL CS                 !
  REAL R                  !
  REAL S                  !
  REAL DFTDQC             !
  REAL DFTDTC             !
  REAL HSTURB             !Sensible heat flux computed by GCM
  REAL ETURB              !Evaporation rate computed by GCM
  REAL DEDQA              !Derivative of evaporation rate
  !  w.r.t. specific humidity
  REAL DEDTC              !Derivative of evaporation rate
  !  w.r.t. surface-canopy system
  !  temperature
  REAL DHSDQA           !
  REAL DHSDTC           !Derivative of sensible heat flux
  !  w.r.t. surface-canopy system
  !  temperature
  REAL DEDEA            !
  REAL DHSDEA           !
  REAL AVISDR           !Visible direct snow albedo
  REAL ANIRDR           !Near infrared direct snow albedo
  REAL AVISDF           !Visible diffuse snow albedo
  REAL ANIRDF           !Near infrared diffuse snow albedo
  REAL ALBAVE           !Average snow albedo
  REAL SWNET            !Net shortwave radiation absorbed
  !  by the surface
  REAL PAR              !Set to 0.5*SWDN
  REAL PDIR             !Set to a value of 0.5
  REAL TOTALB           !Set equal to ALBAVES
  REAL ZLAI             !Set equal to LAI
  REAL QSATTC           !The saturated specific humidity
  !  based on the surface-canopy 
  !  system temperature
  REAL DQSDTC           !The derivative of the saturated
  !  specific humidity with respect
  !  to the surface-canopy system
  !  temperature
  REAL PARDIR           !The direct component of the 
  !  photosynthetically active
  !  radiation (W m-2)
  REAL PARDIF           !The diffuse component of the 
  !  photosynthetically active
  !  radiation (W m-2)
  REAL CD               !Drag coefficient
  REAL LWATTC           !
  REAL DLWDTC           !
  REAL ALWRAD           !First term in longwave radiation
  !  linearization (W m-2)
  REAL BLWRAD           !Second term in longwave radiation
  !  linearization (W m-2K)
  REAL EVAP             !Evaporation rate (W m-2)
  REAL SHFLUX           !Sensible heat flux (W m-2)
  REAL RUNOFF           !Total runoff generated (Sum of
  !  surface runoff and moisture
  !  diffusion flux at the bottom of
  !  the lowest soil layer (kg m-2 s-1ec)
  REAL BOMB             !Flag set to 0 in SUBROUTINE TILE
  REAL ESOI             !The evaporation rate from bare soil (W m-2)
  REAL EINT             !The interception loss (W m-2)
  REAL ESNO             !The evaporation rate from snowpack (W m-2)
  REAL EVEG             !The transpiration rate (W m-2)
  REAL SMELT            !The rate of snowmelt (kg m-2 s-1ec)
  REAL HLATN            !The latent heat flux (W m-2)
  REAL HLWUP            !The outgoing longwave radiation flux (W m-2)
  REAL GDRAIN           !The diffusion of moisture across the
  !  bottom of the root zone (middle soil layer)
  REAL RUNSRF           !The overland flow (kg m-2 s-1ec)
  REAL FWSOIL           !The infiltration of rainwater into
  !  the top soil layer (kg m-2 s-1ec)
  REAL STRDG1           !Diagnostic not currently used
  REAL STRDG2             !Diagnostic not currently used
  REAL STRDG3             !Diagnostic not currently used
  REAL STRDG4             !Diagnostic not currently used
  REAL WSMAX1             !Parameter derived from VGWMAX
  REAL WSMAX2             !Parameter derived from VGWMAX
  REAL WSMAX3             !Parameter derived from VGWMAX
  
  REAL CONDRY             !
  
  REAL GHFLUX             !Ground heat flux (W m-2) ?
  
  REAL POR1               !Porosity of first soil layer?
  REAL POR2               !Porosity of second soil layer?
  REAL POR3               !Porosity of third soil layer?
  REAL WATER1             !SOWET(1)*WSMAX1
  REAL WATER2             !SOWET(2)*WSMAX2
  REAL WATER3             !SOWET(3)*WSMAX3
  REAL  DUMBSNOW     !for the Call to SIBALB
  INTEGER DUMBIRUN    !for the Call to SIBALB
  
  REAL  XV2WIND     !for calculating the wind speed from vector amounts
  REAL  YV2WIND     !for calculating the wind speed from vector amounts
  
  Real  Wiltpoint1,Wiltpoint2,Wiltpoint3
  REAL EPSILON            !
  PARAMETER (ALHE = 2.4548e6)
  PARAMETER (ALHS = 2.8368e6)
  PARAMETER (ALHM = ALHS-ALHE)
  PARAMETER (EPSILON = 18.01/28.97)


!   MOSAIC MONTHLY VEGETATION PARAMETERS
!   HOWEVER, these have already been interpolated
!   So the value passed in is singular for each variable
  !
  REAL GREEN          ! greenness fraction of vegetation
                      !  Holds 12 months of data for interpolation
  REAL LAI            ! leaf area index of vegetation
                      !  Table 5 in Mosaic manual
  REAL VGZ0           !Monthly roughness length
                      !  Table 7 in Mosaic manual
                      !  Used to derive Z0 and U2FAC
  REAL VGROTL         !Monthly variation of root length density
                      !  Table 9 in Mosaic manual
                      !  Used to derive RSOIL1 and RSOIL2
  REAL VGRDC          !Monthly variation of vegetation specific
                      !  constant used to determin subcanopy
                      !  aerodynamic resistance.
                      !  Table 8 in Mosaic manual
                      !  Used to derive RDC
  REAL VGDD           !Monthly variation of zero plane
                                !  displacement height
                                !  Table 10 in Mosaic manual
                                !  Used to derive U2FAC




!   MOSAIC STATIC PARAMETERS (Passed into mos_middle in MOS%VEGP array)
!  THESE ARE THE STATIC VEGETATION ONES
  REAL VGRF11             !Vegetation parameter used to derive SQSCAT
  REAL VGRF12             !Vegetation parameter used to derive SQSCAT
  REAL VGTR11             !Vegetation parameter used to derive SQSCAT
  REAL VGTR12             !Vegetation parameter used to derive SQSCAT
  REAL VGZ2               !Height of the canopy, Z2 is set to this,
                          !  and it's used to derive U2FAC
  REAL VGROTD             !Rooting depth parameter used to derive 
                          !  RSOIL1 and RSOIL2
  REAL VGRDRS             !Resistance to moisture transport per unit
                          !  root length used to derive RSOIL1
  REAL VGROCA             !Average cross sectional area of root,
                          !  Used to derive RSOIL2
  REAL VKC                !Parameter used to derive U2FAC (same for
                                !  all vegetation types (.35))
  REAL VGPH1X             !Soil moisture potential above which
                                !  vegetation is not moisture-stressed (m)
  REAL VGPH2X             !Soil moisture potential below which
                                !  transpiration ceases due to wilting (m)
  REAL VGRPLX             !Average resistance to moisture transport
                                !  within the vegetation itself (s)
  REAL VGCHIL             !Parameter describing departure of leaf
                                !  angles from a spherical distribution
  REAL VGZMEW             !
  REAL VGRST1             !Stomatal resistance parameter a
  REAL VGRST2             !Stomatal resistance parameter b
  REAL VGRST3             !Stomatal resistance parameter c
  REAL VGDFAC             !Parameter controlling vapor pressure
                                !  deficit stress (1/mb)
  REAL VGTLL              !Temperature below which temperature
                                !  stress prevents transpiration (K)
  REAL VGTU               !Temperature above which temperature
                                !  stress prevents transpiration (K)
  REAL VGTCF1             !Coefficient in temperature stress
                                !  equation (K-4)
  REAL VGTCF2             !Coefficient in temperature stress
                                !  equation (K-3)
  REAL VGTCF3             !Coefficient in temperature stress
                                !  equation (K-2)

! MOSAIC STATIC SOIL PARAMETERS

  REAL VGBEEX             !Soil paramter B related to pore size
                                !  distribution index
  REAL VGPSAX             !Soil moisture potential of a saturated
                                !  soil (m)
  REAL VGCSAX             !Hydraulic conductivity of the soil
                                !  at saturation (m s-1)
  REAL VGZDEX(3)          !Soil layer thickness (m) of layer i
  REAL VGSLOX             !Cosine of theta
        
  REAL VGWMAX(3)          !Moisture holding capacity of soil
                                !  layer i (kg m-2)



!   MOSAIC TEMPORALLY INTERPOLATED OR OTHERWISE DERIVED PARAMETERS
!  THESE VALUES ARE GENERALLY BASED ON THE VALUES CONTAINED IN THE
!   VEGPARAM ARRAY OR ON FORCING VARIABLES

!        REAL GREEN              !Current greenness fraction of vegetation
!                                !  Used among other things to derive
!                                !  the SQSCAT parameter
!        REAL LAI                !Current leaf area index of vegetation
  REAL Z0                 !Current roughness parameter (based on VGZ0)
  REAL SQSCAT             !Vegetation parameter derived from
                          !  VGRF11,VGRF12,VGTR11,VGTR12,GREEN
  REAL Z2                 !Height of the canopy, set to VGZ2
  REAL DZM                !Parameter that is set to 600 inside 'PMONTH'
                          !  also used to derive U2FAC
  REAL RSOIL1             !Parameter that is derived from
                          !  VGRDRS, VGROTD, VGROTL
  REAL RSOIL2             !Parameter that is derived from
                          !  VGROCA, VGROTD, VGROTL
  REAL SATCAP             !Parameter that is derived from
                          !  LAI
  REAL RDC                !Parameter derived from VGRDC
  REAL U2FAC              !Parameter derived from VGZZ, VGZ0,
                                !  DZM
  REAL EM                 !Value computed from M_ESAT
  REAL EA                 !Value computed from QA, PSUR
                                !  and EPSILON
  REAL PSUR               !Surface Pressure (Mb), from SFCPRS/100
  REAL DPT2M              !Dew Point (K), derived from HUMID
  REAL EAI                !Variable used to derive DPT2M
  REAL ESA                !Variable used to derive DPT2M
  REAL SUNANG             !Cosine of the Solar Zenith angle computed
                                !  in the ASTRO subroutine
  REAL RA                 !Earth-Sun distance in units of the
                                !  orbits semi-major axis
  REAL TRAINC             !Convective rainfall rate (kg m-2 s-1ec)
                                !  derived from CPCP and TPCP
  REAL TRAINL             !Large-scale rainfall rate (kg m-2 s-1ec)
                                !  derived from TPCP and TRAINC
  REAL TSNOW              !The snowfall rate (kg m-2 s-1ec)
  REAL RAINF
  REAL SNWFRC             !

! MOSAIC FORCING VARIABLES (IN TILE SPACE)
  REAL WND              !Wind Speed (m s-1)
  REAL TMP2M            !2 Meter Temperature (K)
  REAL HUMID            !2 Meter Humidity (kg kg-1) (is corrected for
                        !numerical instability by setting low 
                        !values to .00001
  REAL HUMIDORIG        !2 Meter Humidity (uncorrected version 
                        !of HUMID)
  REAL CPCP             !Convective Precipitation (kg m-2 s-1ec)
  REAL TPCP             !Total Precipitation (kg m-2 s-1ec)
  REAL LWDWN            !Downwelling longwave radiation (W m-2)
  REAL SWDN             !Downwelling shortwave radiation (W m-2)
  REAL SFCPRS           !Surface Pressure (Pa)

  REAL SOILM_PREV
  REAL SOILWM             !Soil wetness variable, max avail volumetric mst
  REAL SOILWW             !Soil wetness variable, actual avail vol. mst
  REAL ROOTTEMP
  integer M
  logical             :: alarmCheck
  integer             :: iret,tid
!********************************        
  real tempac,tempcc
!*********************************        
  character*3         :: fnest
  write(fnest,'(i3.3)') n    

!=== End Variable Definition =============================================

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "Mosaic model alarm "//trim(fnest))
  if(alarmCheck) then 
     do m = 1, LIS_rc%npatch(n,LIS_rc%lsm_index)

        !  RESET was originially used to read data from ETA files
        !  since the LDAS Subdriver has its own restart files
        !  and doesn't have the capability to read ETA soil moisture
        !  data then RESET SHOULD ALWAYS BE ZERO
        !    note that if ETA reading ability is installed
        !   then reset is only equal to one for the very first call
        !   and then should be set to zero again!!!!
        RESET=0

        !  THE FOLLOWING SECTIONS BREAKS DOWN THE THREE PARAMETER ARRAYS
        !   INTO THEIR INDIVIDUAL VARIABLE NAMES FOR LATER USAGE
        !
        !  STATIC VEGETATION PARAMETERS
        VGRF11=MOS_STRUC(N)%MOS(M)%VEGP(1)
        VGRF12=MOS_STRUC(N)%MOS(M)%VEGP(2)
        VGTR11=MOS_STRUC(N)%MOS(M)%VEGP(3)
        VGTR12=MOS_STRUC(N)%MOS(M)%VEGP(4)
        VGZ2  =MOS_STRUC(N)%MOS(M)%VEGP(5)
        VGROTD=MOS_STRUC(N)%MOS(M)%VEGP(6)
        VGRDRS=MOS_STRUC(N)%MOS(M)%VEGP(7)
        VGROCA=MOS_STRUC(N)%MOS(M)%VEGP(8)
        VKC   =MOS_STRUC(N)%MOS(M)%VEGP(9)
        VGPH1X=MOS_STRUC(N)%MOS(M)%VEGP(10)
        VGPH2X=MOS_STRUC(N)%MOS(M)%VEGP(11)
        VGRPLX=MOS_STRUC(N)%MOS(M)%VEGP(12)
        VGCHIL=MOS_STRUC(N)%MOS(M)%VEGP(13)
        VGZMEW=MOS_STRUC(N)%MOS(M)%VEGP(14)
        VGRST1=MOS_STRUC(N)%MOS(M)%VEGP(15)
        VGRST2=MOS_STRUC(N)%MOS(M)%VEGP(16)
        VGRST3=MOS_STRUC(N)%MOS(M)%VEGP(17)
        VGDFAC=MOS_STRUC(N)%MOS(M)%VEGP(18)
        VGTLL =MOS_STRUC(N)%MOS(M)%VEGP(19)
        VGTU  =MOS_STRUC(N)%MOS(M)%VEGP(20)
        VGTCF1=MOS_STRUC(N)%MOS(M)%VEGP(21)
        VGTCF2=MOS_STRUC(N)%MOS(M)%VEGP(22)
        VGTCF3=MOS_STRUC(N)%MOS(M)%VEGP(23)
        SNWMID=MOS_STRUC(N)%MOS(M)%VEGP(24)

        !  VARIABLE VEGETATION PARAMETERS
        GREEN =MOS_STRUC(N)%MOS(M)%VEGIP(1)
        LAI   =MOS_STRUC(N)%MOS(M)%VEGIP(2)
        VGZ0  =MOS_STRUC(N)%MOS(M)%VEGIP(3)
        VGROTL=MOS_STRUC(N)%MOS(M)%VEGIP(4)
        VGRDC =MOS_STRUC(N)%MOS(M)%VEGIP(5)
        VGDD  =MOS_STRUC(N)%MOS(M)%VEGIP(6)

        !  STATIC SOIL PARAMETERS
        VGBEEX   =  MOS_STRUC(N)%MOS(M)%SOILP(1)
        VGPSAX   =  MOS_STRUC(N)%MOS(M)%SOILP(2)
        VGCSAX   =  MOS_STRUC(N)%MOS(M)%SOILP(3)
        VGZDEX(1)=  MOS_STRUC(N)%MOS(M)%SOILP(4)
        VGZDEX(2)=  MOS_STRUC(N)%MOS(M)%SOILP(5)
        VGZDEX(3)=  MOS_STRUC(N)%MOS(M)%SOILP(6)
        VGSLOX   =  MOS_STRUC(N)%MOS(M)%SOILP(7)
        VGWMAX(1)=  MOS_STRUC(N)%MOS(M)%SOILP(8)
        VGWMAX(2)=  MOS_STRUC(N)%MOS(M)%SOILP(9)
        VGWMAX(3)=  MOS_STRUC(N)%MOS(M)%SOILP(10)      

        !  THE FOLLOWING BREAKS DOWN THE FORCING VARIABLES

        TMP2M = MOS_STRUC(N)%MOS(M)%tair/mos_struc(n)%forc_count
        HUMID = MOS_STRUC(N)%MOS(M)%qair/mos_struc(n)%forc_count
        SWDN  = MOS_STRUC(N)%MOS(M)%swdown/mos_struc(n)%forc_count
        LWDWN = MOS_STRUC(N)%MOS(M)%lwdown/mos_struc(n)%forc_count
        XV2WIND=(MOS_STRUC(N)%MOS(M)%uwind)*(MOS_STRUC(N)%MOS(M)%uwind)
        YV2WIND=(MOS_STRUC(N)%MOS(M)%vwind)*(MOS_STRUC(N)%MOS(M)%vwind)
        WND   = SQRT( XV2WIND + YV2WIND )/mos_struc(n)%forc_count
        SFCPRS= MOS_STRUC(N)%MOS(M)%psurf/mos_struc(n)%forc_count
        TPCP  = MOS_STRUC(N)%MOS(M)%rainf/mos_struc(n)%forc_count
        CPCP  = MOS_STRUC(N)%MOS(M)%rainf_c/mos_struc(n)%forc_count
        
        !=== Prevent Numerical Instability
        if(WND.le.0.01) WND=0.01

        !=== Prevent Numerical Instability with HUMID
        !  Store original variable for output
        HUMIDORIG=HUMID
        if(HUMID.le.0.00001)  HUMID=0.00001    ! 1.0E-5  or 0.1E-4

        ! Compute number of seconds into day
        SEC=(LIS_rc%MN*60)+(LIS_rc%HR*60*60)  !total sec into the day

        ! Compute Zenith angle of sun, SUNANG
        !      CALL ASTRO(LIS_rc%YR,LIS_rc%MO,LIS_rc%DA,SEC, &
        !       LIS_DOMAIN(N)%TILE(M)%LAT,LIS_DOMAIN(N)%TILE(M)%LON,1,SUNANG,RA)
        !       index = gindex(LIS_domain(n)%tile(m)%col,LIS_domain(n)%tile(m)%row)
        index1 = LIS_surface(n,LIS_rc%lsm_index)%tile(m)%index

        CALL ASTRO(LIS_rc%YR,LIS_rc%MO,LIS_rc%DA,SEC, &
             LIS_domain(n)%grid(index1)%lat,LIS_domain(n)%grid(index1)%lon,1,SUNANG,RA)

        !      if (m .eq. 100000) then
        !          write(LIS_logunit,*) LIS_rc%YR,LIS_rc%MO,LIS_rc%DA,SEC
        !	  write(LIS_logunit,*) LIS_domain(n)%grid(index1)%lat,LIS_domain(n)%grid(index1)%lon
        !	  write(LIS_logunit,*) SUNANG,RA
        !      endif				  
        !  Call the calculation of the albedo for the given tile, Sunang, veg type
        !  This is in the fortran program calc_albedo.f
        DUMBSNOW=1.0
        DUMBIRUN=1                  !this should stay equal to one
        !      CALL oldSIBALB( avisdr,anirdr, !albedo Visible/Near Ir both direct
        !     & avisdf,anirdf,             !albedo Vis/Near IR both Diffuse
        !     & LAI,GREEN,SUNANG,          !LAI and Green, Solar Zenith Angle
        !     & DUMBSNOW,TILE%VEGT,        !dumby for snow (not used), Vegetation type
        !     & DUMBIRUN,MOS_STRUC(N)%MOS(M)%CT)           !dumby run variable, Temp of Canopy (not used?)

        !      do i=1,13
        !        write(LIS_logunit,*) i, sib%ALVDR(1,1,i)
        !      enddo

        !      stop

        !      write(*,63) m, LIS_domain(n)%tile(m)%row,LIS_domain(n)%tile(m)%col,LIS_DOMAIN(N)%TILE(M)%VEGT,avisdr
        ! 63   format(i6,1x,i4,1x,i4,1x,i2,1x,f8.5)
        call umd_sibalb   &   
             ( avisdr,anirdr,  &       !albedo Visible/Near Ir both direct
             avisdf,anirdf,  &         !albedo Vis/Near IR both Diffuse
             LAI,GREEN,SUNANG, &        !LAI and Green, Solar Zenith Angle
             DUMBSNOW,LIS_surface(N,LIS_rc%lsm_index)%TILE(M)%VEGT, &      !dumby for snow(not used), Vegetation type
             DUMBIRUN,MOS_STRUC(N)%MOS(M)%CT)           !dumby run variable, Temp of Canopy
        !                      !sibalb look up table of coefficients
        !       write(LIS_logunit,*) "M = ",M

        !       if (m .eq. 100000) then
        !        write(LIS_logunit,*) LIS_DOMAIN(N)%TILE(M)%VEGT,LAI,GREEN,SUNANG
        ! write(LIS_logunit,*) avisdr,anirdr,avisdf,anirdf
        !       endif


        !      write(*,63) m, LIS_domain(n)%tile(m)%row,LIS_domain(n)%tile(m)%col,LIS_DOMAIN(N)%TILE(M)%VEGT,avisdr
        ! 64   format(i6,1x,i4,1x,i4,1x,i2,1x,f8.5)

        !      write(LIS_logunit,*) "hereM3 ",M,MOS_STRUC(N)%MOS(M)%CT

        ALBEDO1=avisdr
        ALBEDO2=anirdr
        ALBEDO3=avisdf
        ALBEDO4=anirdf
        !to 0 to 1 scale.
        !=== Process Section
        call m_pmonth(LIS_rc%doy,LIS_domain(n)%grid(index1)%lat,green,lai, &
             z0,sqscat,z2,dzm,rsoil1, rsoil2, &
             satcap,rdc,u2fac,vgrdrs,vgz2,vgrotd, &
             vgroca,vgrdc,vgrotl,vgz0,vgdd, &
             vgtr12,vgtr11,vgrf11,vgrf12) 

        !       if (m .eq. 100000) then
        !       write(LIS_logunit,*) MOS_STRUC(N)%MOS(M)%LAT,MOS_STRUC(N)%MOS(M)%QA
        !       endif

        TSNOW=0.
        TRAINC=AMIN1(CPCP,TPCP)
        TRAINL=TPCP-TRAINC

        !=== Determine if Precip is snow
        if(TMP2M.lt.(273.15))THEN
           tsnow=tpcp
           trainl=0.0
           trainc=0.0     
        endif

        psur=sfcprs/100.0
        !=== Compute ESI and ESA to get the Dew Point, DPT2M
        EAI=HUMID/0.622*SFCPRS/1000.0
        ESA=0.6108*EXP((17.27*(TMP2M-273.15))/ &
             ((237.3+(TMP2M-273.15))))
        DPT2M=TMP2M/(1-(LOG(EAI/ESA)/17.27))
        em=m_esat(dpt2m)             !m_esat is a function
        ea=MOS_STRUC(N)%MOS(M)%qa*psur/epsilon
        sunang=amax1(sunang,0.01)

        !      if (isnan(MOS_STRUC(N)%MOS(M)%snow)) then
        !        MOS_STRUC(N)%MOS(M)%snow = 0.0
        !      end if

        !=== account for fact that forcing can lead the land
        snwfrc=(MOS_STRUC(N)%MOS(M)%snow/(MOS_STRUC(N)%MOS(M)%snow+snwmid))
        esattc=(snwfrc*m_qsat(MOS_STRUC(N)%MOS(M)%ct,psur,alhs))+ &
             (1.-snwfrc)*m_qsat(MOS_STRUC(N)%MOS(M)%ct,psur,alhe)*psur/epsilon

        if(ea.gt.esattc.and.ea.gt.em) ea=amax1(em,esattc)
        if(ea.lt.esattc.and.ea.lt.em) ea=amin1(em,esattc)

        !=== Duplicate what happens in the Call of Turb {Sub-process is turb}    
        data grav/9.81/,rgas/287./,cp/1010./
        data vkrmn/0.41/,epsi/0.611/,b/5./,c/5./,d/5./

        !=== Special parameters for this section
        grav=9.81
        rgas=287.0
        cp=1010.0
        vkrmn=0.41
        epsi=0.611
        b=5.0
        c=5.0
        d=5.0
        ! Height from forcing or default of 50.0, see mos_f2t.F90
        zthick=mos_struc(n)%mos(m)%obsz/mos_struc(n)%forc_count
        rhoair=psur*100./(rgas*MOS_STRUC(N)%MOS(M)%ct)
        qm=em*epsi/psur
        qc=ea*epsi/psur
        dtv=TMP2M*(1.+epsi*qm)-MOS_STRUC(N)%MOS(M)%ct*(1.+epsi*qc)
        ri=grav*zthick*dtv/(TMP2M*WND*WND)
        dridtc=-grav*zthick*(1.+epsi*qc)/(TMP2M*WND*WND)
        dridqc=-grav*zthick*epsi*MOS_STRUC(N)%MOS(M)%ct/(TMP2M*WND*WND)
        dummy=alog(zthick/z0+1.)
        cn=vkrmn*vkrmn/(dummy*dummy)

        if(ri.ge.0.) then
           t=sqrt(1.+d*ri)
           fu=1./(1.+2.*b*ri/t)
           ft=1./(1.+3.*b*ri*t)
           dftdri=-3.*b*ft*ft*( (d*ri)/(2.*t) + t )
        endif

        if(ri.lt.0.) then
           cs=cn*sqrt((zthick/z0)+1.)
           r=3.*b*cs*sqrt(-ri)    !   so that it can be used in the sqrt function
           s=1./(1.+c*r)
           t=b*ri*s 
           ft=1.-3.*t        
           dftdri=-1.5*b*s*s*(2.+c*r)
        endif

        dftdqc=dftdri*dridqc
        dftdtc=dftdri*dridtc
        ctt=cn*ft
        hsturb=rhoair*cp*ctt*WND*(MOS_STRUC(N)%MOS(M)%ct-TMP2M)
        eturb=rhoair*ctt*WND*(qc-qm)
        dedqa=rhoair*cn*WND*( dftdqc*(qc-qm) + ft )
        dedtc=rhoair*cn*WND*dftdtc*(qc-qm)
        dhsdqa=rhoair*cp*cn*WND*dftdqc*(MOS_STRUC(N)%MOS(M)%ct-TMP2M)
        dhsdtc=rhoair*cp*cn*WND*( dftdtc*(MOS_STRUC(N)%MOS(M)%ct-TMP2M) + ft )
        dedea=dedqa*epsi/psur
        dhsdea=dhsdqa*epsi/psur
        ! Calculate Offline RA now
        ra=1/(ctt*WND)

        if(dhsdtc.lt.0.) dhsdtc=0.
        if(dhsdea.lt.0.) dhsdea=0.
        if(dedtc.lt.0.) dedtc=0.
        if(dedea.lt.0.) dedea=0.

        !=== Adjustment to eturb and hsturb with dtcanal based on PSAS and BC
        !      if(LIS_rc%rpsas.eq.1)then
        !       eturb=eturb+dedtc*MOS_STRUC(N)%MOS(M)%dtcanal
        !       hsturb=hsturb+dhsdtc*MOS_STRUC(N)%MOS(M)%dtcanal
        !      endif
        if(LIS_rc%usealbedomap(n).ne."none") then 
           tid = LIS_surface(n,LIS_rc%lsm_index)%tile(m)%tile_id           
           albave = LIS_alb(n)%albsf(tid)
        else
           if(MOS_STRUC(N)%MOS(M)%snow .gt. 0.) then
              snwfrc=MOS_STRUC(N)%MOS(M)%snow/(MOS_STRUC(N)%MOS(M)%snow+snwmid)
              avisdr=avisdr*(1.-snwfrc) + 0.85*snwfrc
              anirdr=anirdr*(1.-snwfrc) + 0.50*snwfrc
              avisdf=avisdf*(1.-snwfrc) + 0.85*snwfrc
              anirdf=anirdf*(1.-snwfrc) + 0.50*snwfrc
           endif
           albave=0.25*(avisdr+anirdr+avisdf+anirdf)
        endif

        !      write(LIS_logunit,*) mos_struc(n)%mos(m)%snow, avisdr, anirdr, avisdf, anirdf

        !      write(LIS_logunit,*) "ALB check1 ", swdn, MOS_STRUC(N)%MOS(M)%snow, snwfrc
        !      write(LIS_logunit,*) "ALB check2 ", avisdr, anirdr, avisdf, anirdf
        swnet=(1.-albave)*swdn
        par=0.5*swdn
        pdir=0.5
        totalb=albave

        zlai=LAI
        qm = em * epsilon / psur
        MOS_STRUC(N)%MOS(M)%qa = ea * epsilon / psur
        snwfrc = MOS_STRUC(N)%MOS(M)%snow / ( MOS_STRUC(N)%MOS(M)%snow + snwmid)
        qsattc = (snwfrc*m_qsat(MOS_STRUC(N)%MOS(M)%ct,psur,alhs)+ &
             (1.-snwfrc)*m_qsat(MOS_STRUC(N)%MOS(M)%ct,psur,alhe))
        dqsdtc = qsattc * 5418. / ( MOS_STRUC(N)%MOS(M)%ct * MOS_STRUC(N)%MOS(M)%ct )
        dedqa = dedea * psur / epsilon
        dhsdqa = dhsdea * psur / epsilon
        pardir = par * pdir
        pardif = par * ( 1. - pdir )

        ! Enable option to use acond from forcing
        if ( mos_struc(n)%forcing_ch > 0 ) then
           ra = 1. / mos_struc(n)%mos(m)%acond    
        else 
           mos_struc(n)%mos(m)%acond = ctt*WND
        endif

        ! Calculate drag coefficient based on RA value set 
        cd = 1. / ( WND * ra )

        !=== Compute constants for longwave radiation linearization
        lwattc = 5.67e-8 * MOS_STRUC(N)%MOS(M)%ct*MOS_STRUC(N)%MOS(M)%ct*MOS_STRUC(N)%MOS(M)%ct*MOS_STRUC(N)%MOS(M)%ct
        dlwdtc =  4. * lwattc / MOS_STRUC(N)%MOS(M)%ct
        alwrad = lwattc - dlwdtc * MOS_STRUC(N)%MOS(M)%ct
        blwrad = dlwdtc

        !      if (m .eq. 1) then
        !      write(LIS_logunit,*) "hereM1 ", M,MOS_STRUC(N)%MOS(M)%CT
        !      endif 
        soilm_prev = mos_struc(n)%mos(m)%water1 + &
             mos_struc(n)%mos(m)%water2 + & 
             mos_struc(n)%mos(m)%water3

        CALL MOSTILE(RESET,LIS_rc%NTS(N),TRAINL,TRAINC,TSNOW,WND,ETURB,&
             DEDQA,DEDTC,HSTURB, DHSDQA, DHSDTC,TMP2M,QM,CD,SUNANG,PARDIR, &
             PARDIF,SWNET,LWDWN,PSUR,ZLAI,GREEN,Z2,SQSCAT,RSOIL1,RSOIL2, &
             RDC,U2FAC,QSATTC,DQSDTC,ALWRAD,BLWRAD,MOS_STRUC(N)%MOS(M)%DTCANAL, &
             MOS_STRUC(N)%MOS(M)%CT,MOS_STRUC(N)%MOS(M)%SoT,&
             MOS_STRUC(N)%MOS(M)%QA,MOS_STRUC(N)%MOS(M)%SoWET(1),&
             MOS_STRUC(N)%MOS(M)%SoWET(2), &
             MOS_STRUC(N)%MOS(M)%SoWET(3),MOS_STRUC(N)%MOS(M)%ICS,&
             MOS_STRUC(N)%MOS(M)%SNOW, &
             EVAP,SHFLUX,RUNOFF,BOMB,EINT,ESOI,EVEG,ESNO,SMELT,HLATN,HLWUP, &
             GDRAIN,RUNSRF,FWSOIL,STRDG1,STRDG2,STRDG3,STRDG4,WSMAX1,WSMAX2, &
             WSMAX3,CONDRY,GHFLUX,POR1,POR2,POR3,water1,water2,water3,VGBEEX, &
             VGPSAX,VGCSAX,VGZDEX,VGSLOX,VGPH1X,VGPH2X,VGRPLX,VGWMAX,SNWMID, &
             VGCHIL,VGZMEW,VGRST1,VGRST2,VGRST3,VGRDRS,VGZ2,VGROTD,VGDFAC, &
             VGTLL,VGTU,VGTCF1,VGTCF2,VGTCF3,tempcc,RA,&
             mos_struc(n)%mos(m)%vegt)

        !=== Restrict snow depth to maximum of 100000.0 kg/m^2 ==================
        !=== Sometimes such great values are assigned over Greenland & N. Canada
        !=== Grib output files do not handle values over this threshold well ====
        if (MOS_STRUC(N)%MOS(M)%SNOW > 100000.0) MOS_STRUC(N)%MOS(M)%SNOW = 100000.0
        call LIS_diagnoseSurfaceOutputVar(n,m,&
             LIS_MOC_SWNET,value=swnet,vlevel=1,unit="W m-2",&
             direction="DN",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_LWNET,&
             value=((-1.0)*(hlwup-lwdwn)),vlevel=1,unit="W m-2",&
             direction="DN",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_QLE,&
             value=hlatn,vlevel=1,unit="W m-2",direction="UP",&
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_QH,&
             value=shflux,vlevel=1,unit="W m-2",direction="UP",&
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_QG,&
             value=ghflux,vlevel=1,unit="W m-2",direction="UP",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SNOWF,&
             value=tsnow,vlevel=1,unit="kg m-2 s-1",direction="DN",&
             surface_type=LIS_rc%lsm_index)
        !Bowen Ratio - sensible/latent
        if(hlatn.gt.0) then 
           call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_BR,value=shflux/(hlatn),&
                vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        else
           call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_BR,&
                value=0.0,vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        endif

        !Evaporative Fraction
        if( (hlatn+shflux).ne.0) then 
           call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_EF,&
                value=hlatn/(hlatn+shflux),&
                vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        else
           !double check
           call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_EF,value=1.0, &
                vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        endif

        if (tsnow.eq.0) then
           rainf=tpcp
        else
           rainf = 0.0
        endif
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_RAINF,value=rainf,&
             vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_EVAP,value=evap,&
             vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_QS,value=runsrf,&
             vlevel=1,unit="kg m-2 s-1",direction="OUT",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_QSB,value=runoff-runsrf,&
             vlevel=1,unit="kg m-2 s-1",direction="OUT",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_QSM,value=smelt,&
             vlevel=1,unit="kg m-2 s-1",direction="S2L",surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_AVGSURFT,&
             value=mos_struc(n)%mos(m)%ct,&
             vlevel=1,unit="K",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_ALBEDO,&
             value=albave,vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SWE,&
             value=mos_struc(n)%mos(m)%snow,&
             vlevel=1,unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SNOWDEPTH,&
             value=mos_struc(n)%mos(m)%snow*10.0/1000.0,&
             vlevel=1,unit="m",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SNOWCOVER,&
             value=snwfrc,vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILTEMP,&
             value=mos_struc(n)%mos(m)%soT,vlevel=1,unit="K",&
             direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_ROOTTEMP,&
             value=mos_struc(n)%mos(m)%soT,vlevel=1,unit="K",&
             direction="-",surface_type=LIS_rc%lsm_index)
        mos_struc(n)%mos(m)%water1 = water1
        mos_struc(n)%mos(m)%water2 = water2
        mos_struc(n)%mos(m)%water3 = water3

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILMOIST,&
             value=water1,vlevel=1,unit="kg m-2",direction="-",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILMOIST,&
             value=water1/(mos_struc(n)%mos(m)%soilp(4)*LIS_CONST_RHOFW),&
             vlevel=1,unit="m^3 m-3",direction="-",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILMOIST,&
             value=water2,vlevel=2,unit="kg m-2",direction="-",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILMOIST,&
             value=water2/(mos_struc(n)%mos(m)%soilp(5)*LIS_CONST_RHOFW),&
             vlevel=2,unit="m^3 m-3",direction="-",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILMOIST,&
             value=water3,vlevel=3,unit="kg m-2",direction="-",&
             surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILMOIST,&
             value=water3/(mos_struc(n)%mos(m)%soilp(6)*LIS_CONST_RHOFW),&
             vlevel=3,unit="m^3 m-3",direction="-",&
             surface_type=LIS_rc%lsm_index)

        !Calculate MSTAV
        soilwm=0.0
        soilww=0.0
        soilwm=por1*vgzdex(1)
        soilwm=soilwm+(por2*vgzdex(2))
        soilwm=soilwm+(por3*vgzdex(3))
        soilww=(mos_struc(n)%mos(m)%sowet(1)*por1)*vgzdex(1)
        soilww=soilww+((mos_struc(n)%mos(m)%sowet(2)*por2)*vgzdex(2))
        soilww=soilww+((mos_struc(n)%mos(m)%sowet(3)*por3)*vgzdex(3))
        mos_struc(n)%mos(m)%soilwet=soilww/soilwm

        wiltpoint1=((VGPH1X/VGPSAX) ** (1.0 / (-VGBEEX)))
        wiltpoint2=((VGPH2X/VGPSAX) ** (1.0 / (-VGBEEX)))
        wiltpoint3=0.0

        !Calculate MSTAVRZ (in excess of wilting point)
        soilwm=0.0
        soilww=0.0
        soilwm=(por1-(por1*wiltpoint2))*vgzdex(1)
        soilwm=soilwm+((por2 - (por2*wiltpoint2))*vgzdex(2))
        soilww=((mos_struc(n)%mos(m)%sowet(1)*por1)-(por1*wiltpoint2))*vgzdex(1)
        soilww=soilww+((mos_struc(n)%mos(m)%sowet(2)*por2) - &
             (por2*wiltpoint2))*vgzdex(2)
        mos_struc(n)%mos(m)%soilwetrz = soilww/soilwm

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_SOILWET,value=soilww/soilwm,&
             vlevel=1,unit="-",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_ECANOP,&
             value=eint*(-1.0)/2.4548E6,&
             vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_TVEG,value=eveg*(-1.0)/2.4548E6,&
             vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_ESOIL,value=esoi*(-1.0)/2.4548E6,&
             vlevel=1,unit="kg m-2 s-1",direction="UP",surface_type=LIS_rc%lsm_index)
        ! 1m soil column
!        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_ROOTMOIST,&
!             value=(water1*20.0+&
!             water2*980.0)/1000.0,vlevel=1,unit="m^3 m-3",direction="-",&
!             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_ROOTMOIST,&
             value=(water1*20.0/(mos_struc(n)%mos(m)%soilp(4)*LIS_CONST_RHOFW)+&
             water2*980.0/(mos_struc(n)%mos(m)%soilp(5)*LIS_CONST_RHOFW))/&
             1000.0,vlevel=1,unit="m^3 m-3",direction="-",&
             surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_CANOPINT,&
             value=mos_struc(n)%mos(m)%ics,&
             vlevel=1,unit="kg m-2",direction="-",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_ACOND,value=tempac,&
             vlevel=1,unit="m s-1",direction="-",surface_type=LIS_rc%lsm_index)

        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_DELSOILMOIST,value=&
             water1+water2+water3-soilm_prev,&
             vlevel=1,unit="kg m-2",direction="INC",surface_type=LIS_rc%lsm_index)
        call LIS_diagnoseSurfaceOutputVar(n,m,LIS_MOC_DELSWE,value=&
             mos_struc(n)%mos(m)%snow-mos_struc(n)%mos(m)%swe_prev,&
             vlevel=1,unit="kg m-2",direction="INC",surface_type=LIS_rc%lsm_index)
        if(LIS_rc%tscount(n) == 0 .or. LIS_rc%tscount(n) ==1 &
             .or.LIS_rc%rstflag(n).eq.1) then
           mos_struc(n)%mos(m)%soilm_prev = mos_struc(n)%mos(m)%soilm_prev + &
                mos_struc(n)%mos(m)%soilmoist1+ & 
                mos_struc(n)%mos(m)%soilmoist2 + mos_struc(n)%mos(m)%soilmoist3
           mos_struc(n)%mos(m)%swe_prev = mos_struc(n)%mos(m)%swe
        endif

        !Calculate soil moisture (kg m-2) in top 1meter layer of soil
        !If soil is less than 1 meter deep, scale the soil moisture
        !up so that it has an effective depth of 1 meter of soil
        DEPTH=0.0
        FLAG=0
        DEPTH=DEPTH+VGZDEX(1)
        IF (DEPTH.GE.1) THEN
           FACTOR=(1.0-((DEPTH-1.0)/DEPTH))
           IF (DEPTH.EQ.1.0) THEN
              mos_struc(n)%mos(m)%soilmoist1m=mos_struc(n)%mos(m)%soilmoist1
              FLAG=1
           ELSE
              mos_struc(n)%mos(m)%soilmoist1m=mos_struc(n)%mos(m)%soilmoist1*FACTOR
              FLAG=1
           ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
           DEPTH=DEPTH+VGZDEX(2)
           IF (DEPTH.GE.1.0) THEN
              IF (DEPTH.EQ.1.0) THEN
                 mos_struc(n)%mos(m)%soilmoist1m =mos_struc(n)%mos(m)%soilmoist1 + &
                      mos_struc(n)%mos(m)%soilmoist2
                 FLAG=1
              ELSE
                 TEMP1=DEPTH-1.0
                 FACTOR=(1.0-(TEMP1/VGZDEX(2)))
                 TEMP2=mos_struc(n)%mos(m)%soilmoist2*FACTOR
                 mos_struc(n)%mos(m)%soilmoist1m=mos_struc(n)%mos(m)%soilmoist1+TEMP2
                 FLAG=1
              ENDIF
           ENDIF
        ENDIF
        IF (FLAG.EQ.0) THEN
           DEPTH=DEPTH+VGZDEX(3)
           IF (DEPTH.GE.1.0) THEN
              IF (DEPTH.EQ.1.0) THEN
                 mos_struc(n)%mos(m)%soilmoist1m = mos_struc(n)%mos(m)%soilmoist1 + &
                      mos_struc(n)%mos(m)%soilmoist2 + mos_struc(n)%mos(m)%soilmoist3
                 FLAG=1
              ELSE
                 TEMP1=DEPTH-1.0
                 FACTOR=(1.0-(TEMP1/VGZDEX(3)))
                 TEMP2=mos_struc(n)%mos(m)%soilmoist3*FACTOR
                 mos_struc(n)%mos(m)%soilmoist1m = mos_struc(n)%mos(m)%soilmoist1 + &
                      mos_struc(n)%mos(m)%soilmoist2 + TEMP2
                 FLAG=1
              ENDIF
           ELSE
              mos_struc(n)%mos(m)%soilmoist1m=(mos_struc(n)%mos(m)%soilmoist1+ &
                   mos_struc(n)%mos(m)%soilmoist2 + &
                   mos_struc(n)%mos(m)%soilmoist3)/DEPTH
              FLAG=1
           ENDIF
        ENDIF


        !reset 
        mos_struc(n)%mos(m)%tair = 0 
        mos_struc(n)%mos(m)%qair = 0 
        mos_struc(n)%mos(m)%swdown = 0 
        mos_struc(n)%mos(m)%lwdown = 0
        mos_struc(n)%mos(m)%uwind = 0
        mos_struc(n)%mos(m)%vwind = 0 
        mos_struc(n)%mos(m)%psurf = 0 
        mos_struc(n)%mos(m)%rainf = 0 
        mos_struc(n)%mos(m)%rainf_c = 0 
        mos_struc(n)%mos(m)%obsz = 0 

     enddo
     mos_struc(n)%count = mos_struc(n)%count + 1
     mos_struc(n)%forc_count = 0 
  end if
  return
end subroutine mos_main


     
     
     


