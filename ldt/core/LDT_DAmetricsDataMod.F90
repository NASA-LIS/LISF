!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_DAmetricsDataMod
!BOP
! 
! !MODULE: LDT_DAmetricsDataMod
! 
! !DESCRIPTION:
!  This module defines the data structures for storing various computed
!  statistics
! 
! !REVISION HISTORY:
!  02 Oct 2008: Sujay Kumar; Initial version
!  28 Feb 2022: Mahdi navari modified to save stratified CDFs
!
! !USES:  
  use ESMF
!EOP

  private
 
  type, public :: DAmetricsEntry

     real,    allocatable :: mask(:,:,:)
     real,    allocatable :: maxval(:,:,:)     
     real,    allocatable :: minval(:,:,:)     
     real,    allocatable :: count_drange_total(:,:,:)
     
     integer, allocatable :: cdf_bincounts(:,:,:,:)
     real,    allocatable :: delta(:,:,:)
     real,    allocatable :: xrange(:,:,:,:)
     real,    allocatable :: cdf(:,:,:,:)
     real,    allocatable :: strat_xrange(:,:,:,:)
     real,    allocatable :: strat_cdf(:,:,:,:)


     real,    allocatable :: sx_mu(:,:,:)
     real,    allocatable :: mu(:,:,:)
     integer, allocatable :: count_mu(:,:,:)

     real,    allocatable :: sx_sigma(:,:,:)
     real,    allocatable :: sxx_sigma(:,:,:)
     real,    allocatable :: sigma(:,:,:)
     integer, allocatable :: count_sigma(:,:,:)

     !-----------------------------------------------Y.Kwon
     real,    allocatable :: sx_mu_6am(:,:,:), sx_mu_6pm(:,:,:)
     real,    allocatable :: mu_6am(:,:,:), mu_6pm(:,:,:)
     integer, allocatable :: count_mu_6am(:,:,:), count_mu_6pm(:,:,:)

     real,    allocatable :: sx_sigma_6am(:,:,:), sx_sigma_6pm(:,:,:)
     real,    allocatable :: sxx_sigma_6am(:,:,:), sxx_sigma_6pm(:,:,:)
     real,    allocatable :: sigma_6am(:,:,:), sigma_6pm(:,:,:)
     integer, allocatable :: count_sigma_6am(:,:,:), count_sigma_6pm(:,:,:)
     real,    allocatable :: mask_6am(:,:,:), mask_6pm(:,:,:)
     !-----------------------------------------------Y.Kwon

     integer          :: selectOpt
     character*20     :: standard_name
  end type DAmetricsEntry
  
  type, public :: DAmetrics_struc

     type(ESMF_Alarm) :: outAlarm
     type(ESMF_Alarm) :: tavgAlarm

     logical          :: computeFlag

     type(DAmetricsEntry) :: swnet        ! Net shortwave radiation (surface) (W/m2)
     type(DAmetricsEntry) :: lwnet        ! Net longwave radiation (surface) (W/m2)
     type(DAmetricsEntry) :: qle          ! Latent Heat Flux (W/m2)
     type(DAmetricsEntry) :: qh           ! Sensible Heat Flux (W/m2)
     type(DAmetricsEntry) :: qg           ! Ground Heat Flux (W/m2)
     type(DAmetricsEntry) :: qf           ! Energy of fusion (W/m2)
     type(DAmetricsEntry) :: qv           ! Energy of sublimation (W/m2)
     type(DAmetricsEntry) :: qtau         ! Momentum flux (N/m2)
     type(DAmetricsEntry) :: qa           ! Advective flux (W/m2)
     type(DAmetricsEntry) :: delsurfheat  ! Change in surface heat storage (J/m2)
     type(DAmetricsEntry) :: delcoldcont  ! Change in snow water content (J/m2)
     type(DAmetricsEntry) :: br           ! Bowen Ratio
     type(DAmetricsEntry) :: ef           ! Evaporative Fraction

     type(DAmetricsEntry) :: snowf        ! Snowfall rate (kg/m2s)
     type(DAmetricsEntry) :: rainf        ! Rainfall rate (kg/m2s)
     type(DAmetricsEntry) :: evap         ! Evapotranspiration (kg/m2s)
     type(DAmetricsEntry) :: qs           ! Surface Runoff(kg/m2s)
     type(DAmetricsEntry) :: qrec         ! Recharge from river to the floodplain (kg/m2s)
     type(DAmetricsEntry) :: qsb          ! Subsurface Runoff (kg/m2s)
     type(DAmetricsEntry) :: qsm          ! Snowmelt (kg/m2s)
     type(DAmetricsEntry) :: qfz          ! Refreezing of water in the snowpack (kg/m2s)
     type(DAmetricsEntry) :: qst          ! Snow throughfall (kg/m2s)
     type(DAmetricsEntry) :: delsoilmoist ! DelSoilMoist
     type(DAmetricsEntry) :: delswe       ! DelSWE
     type(DAmetricsEntry) :: delsurfstor  ! Change in surface water storage (kg/m2)
     type(DAmetricsEntry) :: delintercept ! Change in interception storage (kg/m2)
     
     type(DAmetricsEntry) :: snowt        ! Snow surface temperature (K)
     type(DAmetricsEntry) :: vegt         ! Vegetation canopy temperature (K)
     type(DAmetricsEntry) :: baresoilt    ! Temperature of bare soil (K)
     type(DAmetricsEntry) :: avgsurft     ! Average Surface Temperature (K)
     type(DAmetricsEntry) :: radt         ! Surface Radiative Tempearture (K)
     type(DAmetricsEntry) :: albedo       ! Surface Albedo (-)
     type(DAmetricsEntry) :: swe          ! Snow water equivalent (kg/m2)
     type(DAmetricsEntry) :: sweveg       ! SWE intercepted by vegetation (kg/m2)
     type(DAmetricsEntry) :: snowfrac     ! Grid cell snow covered fraction
     type(DAmetricsEntry) :: snowdepth    ! Snow Depth(m)
     type(DAmetricsEntry) :: snowcover    ! Snow cover
     type(DAmetricsEntry) :: surfstor     ! Surface water storage (kg/m2)

     type(DAmetricsEntry) :: soilmoist
     type(DAmetricsEntry) :: soiltemp
     type(DAmetricsEntry) :: sliqfrac   ! fraction of SWE which is in the liquid phase
     type(DAmetricsEntry) :: smliqfrac  ! Average layer fraction of liquid
                                          ! moisture
     type(DAmetricsEntry) :: smfrozfrac ! Average layer fraction of liquid
                                          ! moisture
          
     type(DAmetricsEntry) :: soilwet      ! Total Soil Wetness (-)
     type(DAmetricsEntry) :: soilet       ! Plant transpiration from a particular root layer (W/m2)
     type(DAmetricsEntry) :: z0brd        ! Background (i.e., snow-free) roughness length (m)
     type(DAmetricsEntry) :: roughness    ! Roughness length (m)

     type(DAmetricsEntry) :: potevap      ! Potential Evapotranspiration (kg/m2s)
     type(DAmetricsEntry) :: ecanop       ! Interception evaporation (kg/m2s)
     type(DAmetricsEntry) :: tveg         ! Vegetation transpiration (kg/m2s)
     type(DAmetricsEntry) :: esoil        ! Bare soil evaporation (kg/m2s)
     type(DAmetricsEntry) :: ewater       ! Open water evaporation (kg/m2s)
     type(DAmetricsEntry) :: rootmoist    ! Root zone soil moisture (kg/m2)
     type(DAmetricsEntry) :: canopint     ! Total canopy water storage (kg/m2s)
     type(DAmetricsEntry) :: evapsnow     ! Snow evaporation (kg/m2s)
     type(DAmetricsEntry) :: subsnow      ! Snow sublimation (kg/m2s)
     type(DAmetricsEntry) :: subsurf      ! Sublimation of the snow free area (kg/m2s)
     type(DAmetricsEntry) :: acond        ! Aerodynamic conductance (m/s)

     type(DAmetricsEntry) :: etpndx       ! Ponded water evaporation (kg/m2s)
     type(DAmetricsEntry) :: infxsrt      ! Infiltration excess (kg/m2)
     type(DAmetricsEntry) :: sfcheadrt    ! Surface overland flow head (kg/m2)
     type(DAmetricsEntry) :: totalprecip  ! Total precipitation rate (kg/m2/s)
     type(DAmetricsEntry) :: rainfconv  ! Convective Rainfall rate (kg/m2/s)

     type(DAmetricsEntry) :: tws
     type(DAmetricsEntry) :: vod

     type(DAmetricsEntry) ::  windforc
     type(DAmetricsEntry) ::  rainfforc
     type(DAmetricsEntry) ::  snowfforc
     type(DAmetricsEntry) ::  tairforc
     type(DAmetricsEntry) ::  qairforc
     type(DAmetricsEntry) ::  psurfforc
     type(DAmetricsEntry) ::  swdownforc
     type(DAmetricsEntry) ::  lwdownforc
     type(DAmetricsEntry) ::  directswforc 
     type(DAmetricsEntry) ::  diffuseswforc 
     type(DAmetricsEntry) ::  nwindforc   
     type(DAmetricsEntry) ::  ewindforc   
     type(DAmetricsEntry) ::  fheightforc   
     type(DAmetricsEntry) ::  chforc   
     type(DAmetricsEntry) ::  cmforc   
     type(DAmetricsEntry) ::  emissforc   
     type(DAmetricsEntry) ::  mixratioforc   
     type(DAmetricsEntry) ::  coszenforc   
     type(DAmetricsEntry) ::  albedoforc   

     type(DAmetricsEntry) :: landmask
     type(DAmetricsEntry) :: landcover
     type(DAmetricsEntry) :: sfctype
     type(DAmetricsEntry) :: regmask
     type(DAmetricsEntry) :: lakedepth
     type(DAmetricsEntry) :: lakefrac
     type(DAmetricsEntry) :: soiltype
     type(DAmetricsEntry) :: sandfrac
     type(DAmetricsEntry) :: clayfrac
     type(DAmetricsEntry) :: siltfrac
     type(DAmetricsEntry) :: porosity
     type(DAmetricsEntry) :: soilcolor
     type(DAmetricsEntry) :: elevation
     type(DAmetricsEntry) :: slope
     type(DAmetricsEntry) :: lai
     type(DAmetricsEntry) :: sai
     type(DAmetricsEntry) :: gvf      !Y.Kwon
     type(DAmetricsEntry) :: teff     !Y.Kwon
     type(DAmetricsEntry) :: irrigtype
     type(DAmetricsEntry) :: irrigfrac
     type(DAmetricsEntry) :: snfralbedo
     type(DAmetricsEntry) :: mxsnalbedo
     type(DAmetricsEntry) :: greenness
     type(DAmetricsEntry) :: tempbot
     type(DAmetricsEntry) :: roottemp
     type(DAmetricsEntry) :: climelev
     type(DAmetricsEntry) :: climppt
     type(DAmetricsEntry) :: climtmin
     type(DAmetricsEntry) :: climtmax

     type(DAmetricsEntry) :: ccond 

     type(DAmetricsEntry) :: relsmc 
     type(DAmetricsEntry) :: rhmin 

     type(DAmetricsEntry) :: tairforc_min
     type(DAmetricsEntry) :: tairforc_max

  end type DAmetrics_struc

  type, public :: mdep
     type(DAmetricsEntry), pointer :: dataEntryPtr
  end type mdep
  
  type(DAmetrics_struc),  save :: LDT_DAmetrics
  type(mdep),     pointer    :: LDT_DAmetricsPtr(:)

  public :: LDT_DAmetrics
  public :: LDT_DAmetricsPtr

end module LDT_DAmetricsDataMod
