!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module hyssib_module
!BOP
!
! !MODULE: hyssib_module
!
!
! !DESCRIPTION:
!  The code in this file provides a description of the 
!  data structure containing the hyssib 1-d variables. The 
!  variables specified in the data structure include: 
!
!  \begin{description}
!   \item[ts]
!    time step in seconds (integer)
!   \item[maxt] 
!    maximum number of tiles per grid
!   \item[vegt] 
!    HYSSIB vegetation class index value
!   \item[count] 
!    number of timesteps between outputs
!   \item[snowtcount] 
!    number of timesteps with snow (for averaging)
!   \item[albedocount] 
!    number of timesteps for albedo averaging
!   \item[sliqfraccount]
!    number of timesteps for snow cover averaging
!   \item[vegp]
!    static monthly vegetation parameters
!   \item[vegip]
!    monthly albedo vegetation lookup table values
!   \item[tempbot]
!    bottom temperature value
!   \item[tempbot1]
!    month 1 bottom temp value
!   \item[tempbot2]
!    month 2 bottom temp value
!   \item[tempstd]
!    standard deviation of topography (m)
!   \item[sgfg]
!    flag which indicates if snow model is active
!   \item[sdens]
!    snow density
!   \item[tc]
!    canopy temperature
!   \item[tg]
!    ground temperature 
!   \item[tsn]
!    temperature of the snow pack
!   \item[www]
!    soil moisture (volumetric) for each soil layer (3)
!   \item[capac]
!    canopy/ground water content (1/2)
!   \item[snow]
!    canopy/ground snow water equivalent
!   \item[swnet]
!    net shortwave radiation flux
!   \item[lwnet]
!    net longwave radiation flux
!   \item[qle]
!    latent heat flux
!   \item[qh]
!    sensible heat flux
!   \item[qg]
!    ground heat flux
!   \item[qf]
!    energy of fusion
!   \item[qv]
!    energy of sublimation
!   \item[qtau]
!    momentum flux from/to surface
!   \item[qa]
!    advective energy
!   \item[delsurfheat]
!    change in surface heat storage
!   \item[delcoldcont]
!    change in snow cold content
!   \item[snowf]
!    snowfall
!   \item[rainf]
!    rainfall
!   \item[evap]
!    total evaporation
!   \item[qs]
!    surface runoff
!   \item[qrec]
!    recharge
!   \item[qsb]
!    subsurface runoff
!   \item[qsm]
!    snowmelt
!   \item[qfz]
!    re-freezing of water in the snow
!   \item[qst]
!    snow throughfall
!   \item[delsoilmoist]
!    change in soil moisture
!   \item[delswe]
!    change in snow water equivalent
!   \item[delintercept]
!    change in interception storage 
!   \item[snowt]
!    snow surface temperature (K)
!   \item[vegtc]
!    vegetation canopy temperature
!   \item[radteff]
!    surface radiative temperature
!   \item[albedo]
!    surface albedo
!   \item[swe]
!    snow water equivalent
!   \item[sweveg]
!    snow water equivalent intercepted by veg
!   \item[soilmoist]
!    soil moisture layers 1-3 (kg m-2)
!   \item[soiltemp]
!    soil temperature (K)
!   \item[soilwet]
!    total column soil wetness
!   \item[potevap]
!    potential evaporation
!   \item[ecanop]
!    interception evaporation
!   \item[tveg]
!    transpiration
!   \item[esoil]
!    bare soil evaporation
!   \item[rootmoist]
!    root zone soil moisture
!   \item[canopint]
!    canopy interception
!   \item[subsnow]
!    snow sublimation
!   \item[subsurf]
!    sublimation of snow free area
!   \item[acond]
!    areodynamic conductance
!   \item[ccond]
!    canopy conductance
!   \item[snowfrac]
!    percent of tile covered by snow
!   \item[sliqfrac]
!    snow liquid fraction
!   \item[ect]
!    transpiration from canopy
!   \item[eci]
!    evaporation of canopy interepted water
!   \item[egs]
!    evaporation of ground intercepted water
!   \item[hc]
!    sensible heat flux canopy
!   \item[hg]
!    sensible heat flux ground
!   \item[radnvisdir]
!    net radiation
!   \item[radnvisdif]
!    net radiation 
!   \item[radnnirdir]
!    net radiation
!   \item[radnnifdif]
!    net radiation
!   \item[salbvisdir]
!    visible component of surface albedo - direct
!   \item[salbvisdif]
!    visible component of surface albedo - diffuse
!   \item[salbnirdir]
!    near ir component of surface albedo - direct
!   \item[salbnirdif]
!    near ir component of surface albedo - diffuse
!   \item[radtc]
!    net radiation canopy
!   \item[radtg]
!    net radiation ground
!   \item[watcan]
!    canopy interception
!   \item[watgrd]
!    ground interception
!   \item[snocan]
!    swe in canopy
!   \item[snogrd]
!    swe on ground
!   \item[wet1]
!    soil wetness layer 1
!   \item[wet2]
!    soil wetness layer 2
!   \item[wet3]
!    soil wetness layer 3
!   \item[tair]
!    air temperature forcing
!   \item[qair]
!    specific humidity forcing
!   \item[swdown]
!    downward shortwave forcing
!   \item[lwdown]
!    downward longwave forcing
!   \item[uwind]
!    u-wind component forcing
!   \item[vwind]
!    v-wind component forcing
!   \item[psurf]
!    surface pressure forcing
!   \item[rainf\_in]
!    total rainfall forcing
!   \item[rainf\_cp]
!    convective rainfall forcing
!   \item[snowf\_in]
!    total snowfall forcing
!   \end{description}
!
! !REVISION HISTORY:
!
!  28 Apr 2002: K. Arsenault, added NOAH LSM 2.5 code to LDAS
!  14 Nov 2002: Sujay Kumar, Optimized version for LIS
!  21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
!  25 Aug 2007: Chuck Alonge, Updates for LIS 5.0 compliance
!  27 Oct 2010: David Mocko, changes for HY-SSiB in LIS6.1
!
!EOP
  implicit none
  
  PRIVATE
   type, public :: hyssibdec
      integer :: ts    !Timestep (seconds)
      integer :: maxt  !Maximum tiles per grid
      integer :: vegt  !UMD to HY-SSiB Vegetation Class Value
      integer :: count 
      integer :: snowtcount
      integer :: albedocount
      integer :: sliqfraccount

      real :: vegp(20)  !Static vegetation parameters,dim(HYSSiB_NVEGP)
      real :: vegip(11) !Interpolated monthly parameters,dim(HYSSiB_NVEGIP)

!-----------------------------------------------------------------------
! HY-SSiB-Constant Variables
!-----------------------------------------------------------------------
      real :: tempbot
      real :: tempbot1  
      real :: tempbot2  
      real :: tempstd
!-----------------------------------------------------------------------
! HY-SSiB-Initial Variables
!-----------------------------------------------------------------------
      real :: sgfg
      real :: sdens
!-----------------------------------------------------------------------
! HY-SSiB-State Variables
!-----------------------------------------------------------------------
      real :: tc
      real :: tg
      real :: tsn
      real :: td
      real :: www(3)
      real :: capac(2)
      real :: snow(2)
!-----------------------------------------------------------------------
! HY-SSiB-Output Variables
!-----------------------------------------------------------------------
      real :: swnet
      real :: lwnet
      real :: qle
      real :: qh
      real :: qg
      real :: qf
      real :: qv
      real :: qtau
      real :: qa
      real :: delsurfheat
      real :: delcoldcont
      real :: snowf
      real :: rainf
      real :: evap
      real :: qs
      real :: qrec
      real :: qsb
      real :: qsm
      real :: qfz
      real :: qst
      real :: delsoilmoist
      real :: delswe
      real :: delintercept
      real :: snowt
      real :: vegtc
      real :: green
      real :: baresoilt
      real :: avgsurft
      real :: radteff
      real :: albedo
      real :: swe
      real :: sweveg
      real :: soilmoist(3)
      real :: soilmoist1m
      real :: soiltemp
      real :: soilwet
      real :: soilwetrz
      real :: potevap
      real :: ecanop
      real :: tveg
      real :: esoil
      real :: rootmoist
      real :: canopint
      real :: subsnow
      real :: subsurf
      real :: acond
      real :: ccond
      real :: snowfrac
      real :: snowdepth
      real :: sliqfrac
      real :: ect
      real :: eci
      real :: egs
      real :: hc
      real :: hg
      real :: radnvisdir
      real :: radnvisdif
      real :: radnnirdir
      real :: radnnirdif
      real :: salbvisdir
      real :: salbvisdif
      real :: salbnirdir
      real :: salbnirdif
      real :: radtc
      real :: radtg
      real :: watcan
      real :: watgrd
      real :: snocan
      real :: snogrd
      real :: wet1
      real :: wet2
      real :: wet3
      real :: swdown_out
      real :: lwdown_out
!-----------------------------------------------------------------------
! HY-SSiB-Input Variables
!-----------------------------------------------------------------------
      real :: uwind
      real :: vwind
      real :: rainf_in
      real :: rainf_cp
      real :: snowf_in
      real :: tair
      real :: qair
      real :: psurf
      real :: swdown
      real :: lwdown
      real :: fheight
!=== End Variable List =================================================
   end type hyssibdec

end module hyssib_module

