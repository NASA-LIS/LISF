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
! !MODULE: clsmf25_types
! 
! !DESCRIPTION: 
!
! definition of types and associated operators for Catchment Model
!
! IMPORTANT:
! When adding a field to any of the derived types, must also update
! the associated assignment and operator definitions.
! THERE IS NO WARNING/ERROR IF OPERATOR IS NOT DEFINED FOR ALL FIELDS!
!
! !REVISION HISTORY: 
! 21 May 2003, Rolf Reichle: Initial Specification
! 25 Jan 2005, Rolf Reichle: Added cat_force_type
! 23 Nov 2012: David Mocko, Added greenness/lai from parameter file
!
!EOP
module clsmf25_types

  use clsmf25_constants  

  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: N_cat_progn
  public :: cat_progn_type, cat_diagn_type, cat_param_type, cat_force_type,&
       cat_output_type

  !ag(01Jan2021)
  public :: cat_route_type 
 
  public :: assignment (=), operator (/), operator (+)
  
  ! -------------------------------------------------------------------------
  !
  ! N_cat_progn = total # states in Catchment model
  !    (including N_snow states for snow and
  !     N_gt states for ground temperature)
  
  integer, parameter :: N_cat_progn = 25  
  
!sm  integer, parameter :: N_snow   = 3   ! # layers in snow model
!sm  integer, parameter :: N_gt     = 6   ! # layers in ground temperature model
  
  ! --------------------------------------------------------------------------
  
  ! Catchment model prognostic variables
  
  type :: cat_progn_type

     ! horizontally, the surface is divided into four fractions:
     !
     ! "1" - saturated
     ! "2" - unsaturated but not stressed
     ! "4" - stressed
     ! "S" - snow
     !
     ! ------------------------------------------------------------
     
     real :: tc1     ! surface/canopy temperature 
     real :: tc2
     real :: tc4
     
     real :: qa1     ! specific humidity in canopy air
     real :: qa2
     real :: qa4
     
     real :: capac   ! canopy interception water
     
     real :: catdef  ! catchment deficit
     real :: rzexc   ! root zone excess
     real :: srfexc  ! surface excess
     
     real, dimension(N_gt)   :: ght     ! ground heat content
     
     real, dimension(N_snow) :: wesn    ! snow water equivalent
     real, dimension(N_snow) :: htsn    ! snow heat content
     real, dimension(N_snow) :: sndz    ! snow depth
  
  end type cat_progn_type
  
  ! ---------------------------------------------------------
  
  ! Catchment model diagnostic variables
  
  type :: cat_diagn_type
     
     real :: ar1      ! area fraction of saturated zone
     real :: ar2      ! area fraction of unsaturated and unstressed zone

     real :: asnow    ! area fraction of snow
     
     real :: sfmc     ! surface moisture content
     real :: rzmc     ! root zone moisture content
     real :: prmc     ! profile moisture content
     
     real :: rzeq     ! root zone water computed from catdef only
     
     real :: tsurf    ! mean surface temperature over entire catchment 
     
     real, dimension(N_gt)   :: tp     ! temperature of soil layers
     
     real, dimension(N_snow) :: tpsn   ! temperature of snow layers
     
     real :: shflux   ! sensible heat flux
     real :: lhflux   ! total latent heat flux
     real :: ghflux   ! ground heat flux to top soil layer
     
     real :: evap     ! total evaporation 
     real :: eint     ! interception loss
     real :: esoi     ! evaporation from bare soil
     real :: eveg     ! transpiration 
     real :: esno     ! evaporation from snow
     
     real :: runoff   ! total runoff
     real :: runsrf   ! surface runoff
     real :: bflow    ! baseflow
     
     real :: snmelt   ! snow melt
     
     real :: lwup     ! outgoing/upward longwave radiation
     real :: swup     ! outgoing/upward shortwave radiation
     
     real :: qinfil   ! infiltration
     
     real :: totalb   ! total albedo

     real :: entot    ! Heat content from subroutine catchment
     real :: wtot     ! Water content from subroutine catchment

     real :: totlwc_prev   ! liquid water content
     real :: totlwc   ! liquid water content

     real :: swland   ! liquid water content
          
  end type cat_diagn_type
  
  ! ---------------------------------------------------------
  
  ! model parameters
  
  type cat_param_type
     
     real :: dpth  ! depth to bedrock from data file (dpth/=dzpr in general!)
     
     ! layer thicknesses for soil moisture model (in [mm]!!!!) 
     
     real :: dzsf  ! "surface layer"                     formerly zdep1
     real :: dzrz  ! "root zone layer"                   formerly zdep2
     real :: dzpr  ! "profile layer" (unsaturated zone)  formerly zdep3
     
     ! layer thicknesses for ground temperature model (in [m]!!!!)
     !
     ! dzgt SHOULD REPLACE data dz /.../ STATEMENT IN gndtp0 AND gndtmp
     !
     real, dimension(N_gt) :: dzgt  
     
     ! soil hydraulic parameters
     
     real :: poros   ! porosity
     real :: cond    ! saturated hydraulic conductivity
     real :: psis    ! Clapp-Hornberger parameter
     real :: bee     ! Clapp-Hornberger parameter

     real :: wpwet   ! wilting poing wetness

     real :: gnu     ! vertical decay factor for transmissivity
     
     ! constant parameters related to vegetation
     
     real :: vgwmax      ! max amount of water available to vegetation 
     
     integer :: vegcls   ! vegetation class
     integer :: cat_id   ! catchment id
     
     ! parameters specific to Catchment Model
     
     real :: bf1
     real :: bf2
     real :: bf3
     real :: cdcr1
     real :: cdcr2
     real :: ars1 
     real :: ars2
     real :: ars3
     real :: ara1
     real :: ara2
     real :: ara3
     real :: ara4
     real :: arw1
     real :: arw2
     real :: arw3
     real :: arw4
     real :: tsa1
     real :: tsa2
     real :: tsb1
     real :: tsb2
     real :: atau
     real :: btau
     
! Added climatological mean values from parameter file - dmm
     real :: grn(12)
     real :: GREENMIN
     real :: GREENMAX
     real :: LAIMIN
     real :: LAIMAX

     real :: lai(12)

  end type cat_param_type
  
  ! ---------------------------------------------------------
  !
  ! input forcings (or boundary conditions) and related variables 
  !
  ! horizontally, the surface is divided into four fractions:
  !
  ! "1" - saturated
  ! "2" - unsaturated but not stressed
  ! "4" - stressed
  ! "S" - snow
  
  type cat_force_type
     
     real :: TRAINC        ! convective rain rate
     real :: TRAINL        ! large-scale rain rate
     real :: TSNOW         ! snowfall
     real :: UM            ! wind
     real :: ETURB1 
     real :: DEDQA1 
     real :: DEDTC1 
     real :: HSTURB1
     real :: DHSDQA1 
     real :: DHSDTC1
     real :: ETURB2 
     real :: DEDQA2 
     real :: DEDTC2 
     real :: HSTURB2
     real :: DHSDQA2 
     real :: DHSDTC2
     real :: ETURB4 
     real :: DEDQA4 
     real :: DEDTC4 
     real :: HSTURB4
     real :: DHSDQA4 
     real :: DHSDTC4
     real :: ETURBS 
     real :: DEDQAS 
     real :: DEDTCS 
     real :: HSTURBS
     real :: DHSDQAS 
     real :: DHSDTCS
     real :: TM            ! 2m temperature
     real :: QM            ! 2m humidity
     real :: ra1 
     real :: ra2 
     real :: ra4 
     real :: raS 
     real :: SUNANG        ! sun angle
     real :: PARDIR        ! direct photosynthetically active radiation
     real :: PARDIF        ! diffuse photosynthetically active radiation
     real :: SWNETF        ! net shortwave radiation (?)
     real :: SWNETS        ! net shortwave radiation (?)
     real :: HLWDWN        ! downward longwave radiation
     real :: PSUR          ! surface pressure
     real :: ZLAI          ! leaf area index
     real :: GREEN         ! greenness
     real :: Z2            
     real :: SQSCAT 
     real :: RSOIL1 
     real :: RSOIL2   
     real :: RDC  
     real :: QSAT1 
     real :: DQS1 
     real :: ALW1 
     real :: BLW1
     real :: QSAT2 
     real :: DQS2 
     real :: ALW2 
     real :: BLW2
     real :: QSAT4 
     real :: DQS4 
     real :: ALW4 
     real :: BLW4
     real :: QSATS 
     real :: DQSS 
     real :: ALWS 
     real :: BLWS                  

  end type cat_force_type

  type cat_output_type
     real :: swnet
     real :: lwnet
     real :: qle
     real :: qh
     real :: qg
     real :: qf
     real :: qv
     real :: qa

     real :: rainf
     real :: snowf
     real :: evap
     real :: qs
     real :: qsb
     real :: qsm
     real :: qfz
     real :: qst
     real :: esoil

     real :: avgsurft
!<debug jvg -- catchment testing with Sarith>
     real :: albedo
!</debug jvg -- catchment testing with Sarith>
     real :: swe
     real :: snod
     real :: snocovr
     real :: tveg
     real :: ecanop
     real :: rootmoist
     real :: soilwet
     real :: acond
     real :: albsn
     real :: watertabled
     real :: tws
     
  end type cat_output_type
  
  ! ----------------------------------------------------------------
  !ag(01Jan2021)
  ! Catchment model 2-way coupling variables

  type :: cat_route_type
     !surface water storage units are in m/s (See HYMAP2_routing_run.F90 and noahmp36_getsws_hymap2.F90) 
     real :: rivsto    ! river water storage [m/s]
     real :: fldsto    ! floodplain water storage [m/s]
     real :: fldfrc    ! grid flooded fraction [-]
  end type cat_route_type

  ! ---------------------------------------------------------
  
  interface assignment (=)
     module procedure scalar2cat_diagn
     module procedure scalar2cat_progn
     module procedure scalar2cat_force
     module procedure scalar2cat_param
  end interface
  
  interface operator (/)
     module procedure cat_diagn_div_scalar
     module procedure cat_progn_div_scalar
     module procedure cat_force_div_scalar
  end interface

  interface operator (+)
     module procedure add_cat_diagn
     module procedure add_cat_progn
     module procedure add_cat_force
  end interface
  
contains
  
  subroutine scalar2cat_diagn( cat_diagn, scalar )
    
    implicit none
    
    real, intent(in) :: scalar

    integer          :: i 

    type(cat_diagn_type), intent(out) :: cat_diagn    

    cat_diagn%ar1    = scalar
    cat_diagn%ar2    = scalar

    cat_diagn%asnow  = scalar
    
    cat_diagn%sfmc   = scalar
    cat_diagn%rzmc   = scalar
    cat_diagn%prmc   = scalar    
    cat_diagn%rzeq   = scalar

    cat_diagn%tsurf  = scalar

    do i=1,N_gt
       cat_diagn%tp(i)  = scalar
    end do

    do i=1,N_snow
       cat_diagn%tpsn(i)= scalar
    end do


    cat_diagn%shflux = scalar
    cat_diagn%lhflux = scalar
    cat_diagn%ghflux = scalar

    cat_diagn%evap   = scalar
    cat_diagn%eint   = scalar
    cat_diagn%esoi   = scalar
    cat_diagn%eveg   = scalar
    cat_diagn%esno   = scalar
    
    cat_diagn%runoff = scalar
    cat_diagn%runsrf = scalar
    cat_diagn%bflow  = scalar

    cat_diagn%snmelt = scalar

    cat_diagn%lwup  = scalar
    cat_diagn%swup  = scalar
    
    cat_diagn%qinfil = scalar
    
    cat_diagn%totalb = scalar

    cat_diagn%entot = scalar
    cat_diagn%wtot = scalar

    cat_diagn%totlwc = scalar
    cat_diagn%totlwc_prev = scalar

    cat_diagn%swland = scalar
    
  end subroutine scalar2cat_diagn
  
  ! -----------------------------------------------------------

  function cat_diagn_div_scalar( cat_diagn, scalar )
    
    implicit none

    type(cat_diagn_type)             :: cat_diagn_div_scalar
    type(cat_diagn_type), intent(in) :: cat_diagn
    
    real, intent(in) :: scalar
    
    integer :: i     ! local
    
    cat_diagn_div_scalar%ar1    =     cat_diagn%ar1    / scalar
    cat_diagn_div_scalar%ar2    =     cat_diagn%ar2    / scalar

    cat_diagn_div_scalar%asnow  =     cat_diagn%asnow  / scalar

    cat_diagn_div_scalar%sfmc   =     cat_diagn%sfmc   / scalar
    cat_diagn_div_scalar%rzmc   =     cat_diagn%rzmc   / scalar
    cat_diagn_div_scalar%prmc   =     cat_diagn%prmc   / scalar
    cat_diagn_div_scalar%rzeq   =     cat_diagn%rzeq   / scalar
    
    cat_diagn_div_scalar%tsurf  =     cat_diagn%tsurf  / scalar
    
    do i=1,N_gt
       cat_diagn_div_scalar%tp(i)  =     cat_diagn%tp(i)  / scalar
    end do

    do i=1,N_snow
       cat_diagn_div_scalar%tpsn(i)=     cat_diagn%tpsn(i)/ scalar
    end do
 

    cat_diagn_div_scalar%shflux =     cat_diagn%shflux / scalar
    cat_diagn_div_scalar%lhflux =     cat_diagn%lhflux / scalar
    cat_diagn_div_scalar%ghflux =     cat_diagn%ghflux / scalar

    cat_diagn_div_scalar%evap   =     cat_diagn%evap   / scalar
    cat_diagn_div_scalar%eint   =     cat_diagn%eint   / scalar
    cat_diagn_div_scalar%esoi   =     cat_diagn%esoi   / scalar
    cat_diagn_div_scalar%eveg   =     cat_diagn%eveg   / scalar
    cat_diagn_div_scalar%esno   =     cat_diagn%esno   / scalar

    
    cat_diagn_div_scalar%runoff =     cat_diagn%runoff / scalar
    cat_diagn_div_scalar%runsrf =     cat_diagn%runsrf / scalar
    cat_diagn_div_scalar%bflow  =     cat_diagn%bflow  / scalar

    cat_diagn_div_scalar%snmelt =     cat_diagn%snmelt / scalar
    
    cat_diagn_div_scalar%lwup  =     cat_diagn%lwup  / scalar
    cat_diagn_div_scalar%swup  =     cat_diagn%swup  / scalar

    cat_diagn_div_scalar%qinfil =     cat_diagn%qinfil / scalar
    
    cat_diagn_div_scalar%totalb =     cat_diagn%totalb / scalar

    cat_diagn_div_scalar%entot =     cat_diagn%entot / scalar
    cat_diagn_div_scalar%wtot =     cat_diagn%wtot / scalar    

    cat_diagn_div_scalar%totlwc =     cat_diagn%totlwc / scalar
    cat_diagn_div_scalar%totlwc_prev =     cat_diagn%totlwc_prev / scalar

    cat_diagn_div_scalar%swland =     cat_diagn%swland / scalar

  end function cat_diagn_div_scalar

  ! -----------------------------------------------------------

  function add_cat_diagn( cat_diagn_1, cat_diagn_2 )
    
    implicit none

    type(cat_diagn_type)             :: add_cat_diagn
    type(cat_diagn_type), intent(in) :: cat_diagn_1, cat_diagn_2

    integer :: i     ! local
    
    add_cat_diagn%ar1    =     cat_diagn_1%ar1    +     cat_diagn_2%ar1    
    add_cat_diagn%ar2    =     cat_diagn_1%ar2    +     cat_diagn_2%ar2    
    
    add_cat_diagn%asnow  =     cat_diagn_1%asnow  +     cat_diagn_2%asnow  
    
    add_cat_diagn%sfmc   =     cat_diagn_1%sfmc   +     cat_diagn_2%sfmc  
    add_cat_diagn%rzmc   =     cat_diagn_1%rzmc   +     cat_diagn_2%rzmc  
    add_cat_diagn%prmc   =     cat_diagn_1%prmc   +     cat_diagn_2%prmc  
    add_cat_diagn%rzeq   =     cat_diagn_1%rzeq   +     cat_diagn_2%rzeq   
        
    add_cat_diagn%tsurf  =     cat_diagn_1%tsurf  +     cat_diagn_2%tsurf  

    do i=1,N_gt
       add_cat_diagn%tp(i)  =     cat_diagn_1%tp(i)  +     cat_diagn_2%tp(i)  
    end do
    
    do i=1,N_snow
       add_cat_diagn%tpsn(i)=     cat_diagn_1%tpsn(i)+     cat_diagn_2%tpsn(i)
    end do

    add_cat_diagn%shflux =     cat_diagn_1%shflux +     cat_diagn_2%shflux 
    add_cat_diagn%lhflux =     cat_diagn_1%lhflux +     cat_diagn_2%lhflux
    add_cat_diagn%ghflux =     cat_diagn_1%ghflux +     cat_diagn_2%ghflux 

    add_cat_diagn%evap   =     cat_diagn_1%evap   +     cat_diagn_2%evap   
    add_cat_diagn%eint   =     cat_diagn_1%eint   +     cat_diagn_2%eint   
    add_cat_diagn%esoi   =     cat_diagn_1%esoi   +     cat_diagn_2%esoi   
    add_cat_diagn%eveg   =     cat_diagn_1%eveg   +     cat_diagn_2%eveg   
    add_cat_diagn%esno   =     cat_diagn_1%esno   +     cat_diagn_2%esno   

    add_cat_diagn%runoff =     cat_diagn_1%runoff +     cat_diagn_2%runoff 
    add_cat_diagn%runsrf =     cat_diagn_1%runsrf +     cat_diagn_2%runsrf 
    add_cat_diagn%bflow  =     cat_diagn_1%bflow  +     cat_diagn_2%bflow  
    
    add_cat_diagn%snmelt =     cat_diagn_1%snmelt +     cat_diagn_2%snmelt  

    add_cat_diagn%lwup  =     cat_diagn_1%lwup  +     cat_diagn_2%lwup  
    add_cat_diagn%swup  =     cat_diagn_1%swup  +     cat_diagn_2%swup  

    add_cat_diagn%qinfil =     cat_diagn_1%qinfil +     cat_diagn_2%qinfil 

    add_cat_diagn%totalb =     cat_diagn_1%totalb +     cat_diagn_2%totalb 

    add_cat_diagn%entot =     cat_diagn_1%entot +     cat_diagn_2%entot 
    add_cat_diagn%wtot =     cat_diagn_1%wtot +     cat_diagn_2%wtot     

    add_cat_diagn%totlwc =     cat_diagn_1%totlwc +     cat_diagn_2%totlwc 
    add_cat_diagn%totlwc_prev = cat_diagn_1%totlwc_prev+cat_diagn_2%totlwc_prev

    add_cat_diagn%swland =     cat_diagn_1%swland +     cat_diagn_2%swland 
  end function add_cat_diagn

  ! *******************************************************************
  
  subroutine scalar2cat_progn( cat_progn, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_progn_type), intent(out) :: cat_progn
    
    integer :: i     ! local
    
    cat_progn%tc1    = scalar
    cat_progn%tc2    = scalar
    cat_progn%tc4    = scalar
    cat_progn%qa1    = scalar
    cat_progn%qa2    = scalar
    cat_progn%qa4    = scalar
    cat_progn%capac  = scalar
    cat_progn%catdef = scalar
    cat_progn%rzexc  = scalar
    cat_progn%srfexc = scalar
    
    do i=1,N_gt
       cat_progn%ght(i)  = scalar
    end do
    
    do i=1,N_snow
       cat_progn%wesn(i) = scalar
       cat_progn%htsn(i) = scalar
       cat_progn%sndz(i) = scalar
    end do
    
  end subroutine scalar2cat_progn
  
  ! ---------------------------------------------------
  
  function cat_progn_div_scalar( cat_progn, scalar )
    
    implicit none

    type(cat_progn_type)             :: cat_progn_div_scalar
    type(cat_progn_type), intent(in) :: cat_progn
    
    real, intent(in) :: scalar
    
    integer :: i       ! local

    cat_progn_div_scalar%tc1    =     cat_progn%tc1    / scalar
    cat_progn_div_scalar%tc2    =     cat_progn%tc2    / scalar
    cat_progn_div_scalar%tc4    =     cat_progn%tc4    / scalar
    cat_progn_div_scalar%qa1    =     cat_progn%qa1    / scalar
    cat_progn_div_scalar%qa2    =     cat_progn%qa2    / scalar
    cat_progn_div_scalar%qa4    =     cat_progn%qa4    / scalar
    cat_progn_div_scalar%capac  =     cat_progn%capac  / scalar
    cat_progn_div_scalar%catdef =     cat_progn%catdef / scalar
    cat_progn_div_scalar%rzexc  =     cat_progn%rzexc  / scalar
    cat_progn_div_scalar%srfexc =     cat_progn%srfexc / scalar
    
    do i=1,N_gt
       cat_progn_div_scalar%ght(i)  =     cat_progn%ght(i)  / scalar
    end do
    
    do i=1,N_snow
       cat_progn_div_scalar%wesn(i) =     cat_progn%wesn(i) / scalar
       cat_progn_div_scalar%htsn(i) =     cat_progn%htsn(i) / scalar
       cat_progn_div_scalar%sndz(i) =     cat_progn%sndz(i) / scalar
    end do
    
  end function cat_progn_div_scalar

  ! -----------------------------------------------------------

  function add_cat_progn( cat_progn_1, cat_progn_2 )
    
    implicit none

    type(cat_progn_type)             :: add_cat_progn
    type(cat_progn_type), intent(in) :: cat_progn_1, cat_progn_2

    integer :: i     ! local
    
    add_cat_progn%tc1    =     cat_progn_1%tc1    +     cat_progn_2%tc1
    add_cat_progn%tc2    =     cat_progn_1%tc2    +     cat_progn_2%tc2   
    add_cat_progn%tc4    =     cat_progn_1%tc4    +     cat_progn_2%tc4 
    add_cat_progn%qa1    =     cat_progn_1%qa1    +     cat_progn_2%qa1 
    add_cat_progn%qa2    =     cat_progn_1%qa2    +     cat_progn_2%qa2   
    add_cat_progn%qa4    =     cat_progn_1%qa4    +     cat_progn_2%qa4   
    add_cat_progn%capac  =     cat_progn_1%capac  +     cat_progn_2%capac   
    add_cat_progn%catdef =     cat_progn_1%catdef +     cat_progn_2%catdef 
    add_cat_progn%rzexc  =     cat_progn_1%rzexc  +     cat_progn_2%rzexc  
    add_cat_progn%srfexc =     cat_progn_1%srfexc +     cat_progn_2%srfexc 

    do i=1,N_gt
       add_cat_progn%ght(i)   = cat_progn_1%ght(i)  +   cat_progn_2%ght(i)
    end do
    
    do i=1,N_snow
       add_cat_progn%wesn(i)  = cat_progn_1%wesn(i)  +  cat_progn_2%wesn(i)  
       add_cat_progn%htsn(i)  = cat_progn_1%htsn(i)  +  cat_progn_2%htsn(i)  
       add_cat_progn%sndz(i)  = cat_progn_1%sndz(i)  +  cat_progn_2%sndz(i)  
    end do
    
    
  end function add_cat_progn
  
  ! ****************************************************
    
  subroutine scalar2cat_force( cat_force, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_force_type), intent(out) :: cat_force
    
    cat_force%TRAINC       = scalar
    cat_force%TRAINL       = scalar
    cat_force%TSNOW        = scalar
    cat_force%UM           = scalar
    cat_force%ETURB1       = scalar
    cat_force%DEDQA1       = scalar
    cat_force%DEDTC1       = scalar
    cat_force%HSTURB1      = scalar
    cat_force%DHSDQA1      = scalar
    cat_force%DHSDTC1      = scalar
    cat_force%ETURB2       = scalar
    cat_force%DEDQA2       = scalar
    cat_force%DEDTC2       = scalar
    cat_force%HSTURB2      = scalar
    cat_force%DHSDQA2      = scalar
    cat_force%DHSDTC2      = scalar
    cat_force%ETURB4       = scalar
    cat_force%DEDQA4       = scalar
    cat_force%DEDTC4       = scalar
    cat_force%HSTURB4      = scalar
    cat_force%DHSDQA4      = scalar
    cat_force%DHSDTC4      = scalar
    cat_force%ETURBS       = scalar
    cat_force%DEDQAS       = scalar
    cat_force%DEDTCS       = scalar
    cat_force%HSTURBS      = scalar
    cat_force%DHSDQAS      = scalar
    cat_force%DHSDTCS      = scalar
    cat_force%TM           = scalar
    cat_force%QM           = scalar
    cat_force%ra1          = scalar
    cat_force%ra2          = scalar
    cat_force%ra4          = scalar
    cat_force%raS          = scalar
    cat_force%SUNANG       = scalar
    cat_force%PARDIR       = scalar
    cat_force%PARDIF       = scalar
    cat_force%SWNETF       = scalar
    cat_force%SWNETS       = scalar
    cat_force%HLWDWN       = scalar
    cat_force%PSUR         = scalar
    cat_force%ZLAI         = scalar
    cat_force%GREEN        = scalar
    cat_force%Z2           = scalar
    cat_force%SQSCAT       = scalar
    cat_force%RSOIL1       = scalar
    cat_force%RSOIL2       = scalar
    cat_force%RDC          = scalar
    cat_force%QSAT1        = scalar
    cat_force%DQS1         = scalar
    cat_force%ALW1         = scalar
    cat_force%BLW1         = scalar
    cat_force%QSAT2        = scalar
    cat_force%DQS2         = scalar
    cat_force%ALW2         = scalar
    cat_force%BLW2         = scalar
    cat_force%QSAT4        = scalar
    cat_force%DQS4         = scalar
    cat_force%ALW4         = scalar
    cat_force%BLW4         = scalar
    cat_force%QSATS        = scalar
    cat_force%DQSS         = scalar
    cat_force%ALWS         = scalar
    cat_force%BLWS         = scalar
    
  end subroutine scalar2cat_force
  
  ! -------------------------------------------------------------
  
  function cat_force_div_scalar( cat_force, scalar )
    
    implicit none
    
    type(cat_force_type)             :: cat_force_div_scalar
    type(cat_force_type), intent(in) :: cat_force
    
    real, intent(in) :: scalar
    
    cat_force_div_scalar%TRAINC  = cat_force%TRAINC       / scalar
    cat_force_div_scalar%TRAINL  = cat_force%TRAINL       / scalar
    cat_force_div_scalar%TSNOW   = cat_force%TSNOW        / scalar
    cat_force_div_scalar%UM      = cat_force%UM           / scalar
    cat_force_div_scalar%ETURB1  = cat_force%ETURB1       / scalar
    cat_force_div_scalar%DEDQA1  = cat_force%DEDQA1       / scalar
    cat_force_div_scalar%DEDTC1  = cat_force%DEDTC1       / scalar
    cat_force_div_scalar%HSTURB1 = cat_force%HSTURB1      / scalar
    cat_force_div_scalar%DHSDQA1 = cat_force%DHSDQA1      / scalar
    cat_force_div_scalar%DHSDTC1 = cat_force%DHSDTC1      / scalar
    cat_force_div_scalar%ETURB2  = cat_force%ETURB2       / scalar
    cat_force_div_scalar%DEDQA2  = cat_force%DEDQA2       / scalar
    cat_force_div_scalar%DEDTC2  = cat_force%DEDTC2       / scalar
    cat_force_div_scalar%HSTURB2 = cat_force%HSTURB2      / scalar
    cat_force_div_scalar%DHSDQA2 = cat_force%DHSDQA2      / scalar
    cat_force_div_scalar%DHSDTC2 = cat_force%DHSDTC2      / scalar
    cat_force_div_scalar%ETURB4  = cat_force%ETURB4       / scalar
    cat_force_div_scalar%DEDQA4  = cat_force%DEDQA4       / scalar
    cat_force_div_scalar%DEDTC4  = cat_force%DEDTC4       / scalar
    cat_force_div_scalar%HSTURB4 = cat_force%HSTURB4      / scalar
    cat_force_div_scalar%DHSDQA4 = cat_force%DHSDQA4      / scalar
    cat_force_div_scalar%DHSDTC4 = cat_force%DHSDTC4      / scalar
    cat_force_div_scalar%ETURBS  = cat_force%ETURBS       / scalar
    cat_force_div_scalar%DEDQAS  = cat_force%DEDQAS       / scalar
    cat_force_div_scalar%DEDTCS  = cat_force%DEDTCS       / scalar
    cat_force_div_scalar%HSTURBS = cat_force%HSTURBS      / scalar
    cat_force_div_scalar%DHSDQAS = cat_force%DHSDQAS      / scalar
    cat_force_div_scalar%DHSDTCS = cat_force%DHSDTCS      / scalar
    cat_force_div_scalar%TM      = cat_force%TM           / scalar
    cat_force_div_scalar%QM      = cat_force%QM           / scalar
    cat_force_div_scalar%ra1     = cat_force%ra1          / scalar
    cat_force_div_scalar%ra2     = cat_force%ra2          / scalar
    cat_force_div_scalar%ra4     = cat_force%ra4          / scalar
    cat_force_div_scalar%raS     = cat_force%raS          / scalar
    cat_force_div_scalar%SUNANG  = cat_force%SUNANG       / scalar
    cat_force_div_scalar%PARDIR  = cat_force%PARDIR       / scalar
    cat_force_div_scalar%PARDIF  = cat_force%PARDIF       / scalar
    cat_force_div_scalar%SWNETF  = cat_force%SWNETF       / scalar
    cat_force_div_scalar%SWNETS  = cat_force%SWNETS       / scalar
    cat_force_div_scalar%HLWDWN  = cat_force%HLWDWN       / scalar
    cat_force_div_scalar%PSUR    = cat_force%PSUR         / scalar
    cat_force_div_scalar%ZLAI    = cat_force%ZLAI         / scalar
    cat_force_div_scalar%GREEN   = cat_force%GREEN        / scalar
    cat_force_div_scalar%Z2      = cat_force%Z2           / scalar
    cat_force_div_scalar%SQSCAT  = cat_force%SQSCAT       / scalar
    cat_force_div_scalar%RSOIL1  = cat_force%RSOIL1       / scalar
    cat_force_div_scalar%RSOIL2  = cat_force%RSOIL2       / scalar
    cat_force_div_scalar%RDC     = cat_force%RDC          / scalar
    cat_force_div_scalar%QSAT1   = cat_force%QSAT1        / scalar
    cat_force_div_scalar%DQS1    = cat_force%DQS1         / scalar
    cat_force_div_scalar%ALW1    = cat_force%ALW1         / scalar
    cat_force_div_scalar%BLW1    = cat_force%BLW1         / scalar
    cat_force_div_scalar%QSAT2   = cat_force%QSAT2        / scalar
    cat_force_div_scalar%DQS2    = cat_force%DQS2         / scalar
    cat_force_div_scalar%ALW2    = cat_force%ALW2         / scalar
    cat_force_div_scalar%BLW2    = cat_force%BLW2         / scalar
    cat_force_div_scalar%QSAT4   = cat_force%QSAT4        / scalar
    cat_force_div_scalar%DQS4    = cat_force%DQS4         / scalar
    cat_force_div_scalar%ALW4    = cat_force%ALW4         / scalar
    cat_force_div_scalar%BLW4    = cat_force%BLW4         / scalar
    cat_force_div_scalar%QSATS   = cat_force%QSATS        / scalar
    cat_force_div_scalar%DQSS    = cat_force%DQSS         / scalar
    cat_force_div_scalar%ALWS    = cat_force%ALWS         / scalar
    cat_force_div_scalar%BLWS    = cat_force%BLWS         / scalar

  end function cat_force_div_scalar
  
   ! -----------------------------------------------------------

  function add_cat_force( cat_force_1, cat_force_2 )
    
    implicit none
    
    type(cat_force_type)             :: add_cat_force
    type(cat_force_type), intent(in) :: cat_force_1, cat_force_2
    
    add_cat_force%TRAINC  = cat_force_1%TRAINC      + cat_force_2%TRAINC   
    add_cat_force%TRAINL  = cat_force_1%TRAINL      + cat_force_2%TRAINL   
    add_cat_force%TSNOW   = cat_force_1%TSNOW       + cat_force_2%TSNOW    
    add_cat_force%UM      = cat_force_1%UM          + cat_force_2%UM       
    add_cat_force%ETURB1  = cat_force_1%ETURB1      + cat_force_2%ETURB1   
    add_cat_force%DEDQA1  = cat_force_1%DEDQA1      + cat_force_2%DEDQA1   
    add_cat_force%DEDTC1  = cat_force_1%DEDTC1      + cat_force_2%DEDTC1   
    add_cat_force%HSTURB1 = cat_force_1%HSTURB1     + cat_force_2%HSTURB1  
    add_cat_force%DHSDQA1 = cat_force_1%DHSDQA1     + cat_force_2%DHSDQA1  
    add_cat_force%DHSDTC1 = cat_force_1%DHSDTC1     + cat_force_2%DHSDTC1  
    add_cat_force%ETURB2  = cat_force_1%ETURB2      + cat_force_2%ETURB2   
    add_cat_force%DEDQA2  = cat_force_1%DEDQA2      + cat_force_2%DEDQA2   
    add_cat_force%DEDTC2  = cat_force_1%DEDTC2      + cat_force_2%DEDTC2   
    add_cat_force%HSTURB2 = cat_force_1%HSTURB2     + cat_force_2%HSTURB2  
    add_cat_force%DHSDQA2 = cat_force_1%DHSDQA2     + cat_force_2%DHSDQA2  
    add_cat_force%DHSDTC2 = cat_force_1%DHSDTC2     + cat_force_2%DHSDTC2  
    add_cat_force%ETURB4  = cat_force_1%ETURB4      + cat_force_2%ETURB4   
    add_cat_force%DEDQA4  = cat_force_1%DEDQA4      + cat_force_2%DEDQA4   
    add_cat_force%DEDTC4  = cat_force_1%DEDTC4      + cat_force_2%DEDTC4   
    add_cat_force%HSTURB4 = cat_force_1%HSTURB4     + cat_force_2%HSTURB4  
    add_cat_force%DHSDQA4 = cat_force_1%DHSDQA4     + cat_force_2%DHSDQA4  
    add_cat_force%DHSDTC4 = cat_force_1%DHSDTC4     + cat_force_2%DHSDTC4  
    add_cat_force%ETURBS  = cat_force_1%ETURBS      + cat_force_2%ETURBS   
    add_cat_force%DEDQAS  = cat_force_1%DEDQAS      + cat_force_2%DEDQAS   
    add_cat_force%DEDTCS  = cat_force_1%DEDTCS      + cat_force_2%DEDTCS   
    add_cat_force%HSTURBS = cat_force_1%HSTURBS     + cat_force_2%HSTURBS  
    add_cat_force%DHSDQAS = cat_force_1%DHSDQAS     + cat_force_2%DHSDQAS  
    add_cat_force%DHSDTCS = cat_force_1%DHSDTCS     + cat_force_2%DHSDTCS  
    add_cat_force%TM      = cat_force_1%TM          + cat_force_2%TM       
    add_cat_force%QM      = cat_force_1%QM          + cat_force_2%QM       
    add_cat_force%ra1     = cat_force_1%ra1         + cat_force_2%ra1      
    add_cat_force%ra2     = cat_force_1%ra2         + cat_force_2%ra2      
    add_cat_force%ra4     = cat_force_1%ra4         + cat_force_2%ra4      
    add_cat_force%raS     = cat_force_1%raS         + cat_force_2%raS      
    add_cat_force%SUNANG  = cat_force_1%SUNANG      + cat_force_2%SUNANG   
    add_cat_force%PARDIR  = cat_force_1%PARDIR      + cat_force_2%PARDIR   
    add_cat_force%PARDIF  = cat_force_1%PARDIF      + cat_force_2%PARDIF   
    add_cat_force%SWNETF  = cat_force_1%SWNETF      + cat_force_2%SWNETF   
    add_cat_force%SWNETS  = cat_force_1%SWNETS      + cat_force_2%SWNETS   
    add_cat_force%HLWDWN  = cat_force_1%HLWDWN      + cat_force_2%HLWDWN   
    add_cat_force%PSUR    = cat_force_1%PSUR        + cat_force_2%PSUR     
    add_cat_force%ZLAI    = cat_force_1%ZLAI        + cat_force_2%ZLAI     
    add_cat_force%GREEN   = cat_force_1%GREEN       + cat_force_2%GREEN    
    add_cat_force%Z2      = cat_force_1%Z2          + cat_force_2%Z2       
    add_cat_force%SQSCAT  = cat_force_1%SQSCAT      + cat_force_2%SQSCAT   
    add_cat_force%RSOIL1  = cat_force_1%RSOIL1      + cat_force_2%RSOIL1   
    add_cat_force%RSOIL2  = cat_force_1%RSOIL2      + cat_force_2%RSOIL2   
    add_cat_force%RDC     = cat_force_1%RDC         + cat_force_2%RDC      
    add_cat_force%QSAT1   = cat_force_1%QSAT1       + cat_force_2%QSAT1    
    add_cat_force%DQS1    = cat_force_1%DQS1        + cat_force_2%DQS1     
    add_cat_force%ALW1    = cat_force_1%ALW1        + cat_force_2%ALW1     
    add_cat_force%BLW1    = cat_force_1%BLW1        + cat_force_2%BLW1     
    add_cat_force%QSAT2   = cat_force_1%QSAT2       + cat_force_2%QSAT2    
    add_cat_force%DQS2    = cat_force_1%DQS2        + cat_force_2%DQS2     
    add_cat_force%ALW2    = cat_force_1%ALW2        + cat_force_2%ALW2     
    add_cat_force%BLW2    = cat_force_1%BLW2        + cat_force_2%BLW2     
    add_cat_force%QSAT4   = cat_force_1%QSAT4       + cat_force_2%QSAT4    
    add_cat_force%DQS4    = cat_force_1%DQS4        + cat_force_2%DQS4     
    add_cat_force%ALW4    = cat_force_1%ALW4        + cat_force_2%ALW4     
    add_cat_force%BLW4    = cat_force_1%BLW4        + cat_force_2%BLW4     
    add_cat_force%QSATS   = cat_force_1%QSATS       + cat_force_2%QSATS    
    add_cat_force%DQSS    = cat_force_1%DQSS        + cat_force_2%DQSS     
    add_cat_force%ALWS    = cat_force_1%ALWS        + cat_force_2%ALWS     
    add_cat_force%BLWS    = cat_force_1%BLWS        + cat_force_2%BLWS     

  end function add_cat_force
  
  ! ************************************************************
  
  subroutine scalar2cat_param( cat_param, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(cat_param_type), intent(out) :: cat_param
    
    integer :: i     ! local

    ! ---------------------

    cat_param%dpth   = scalar 
    
    cat_param%dzsf   = scalar 
    cat_param%dzrz   = scalar 
    cat_param%dzpr   = scalar 
    
    do i=1,N_gt
       cat_param%dzgt(i) = scalar
    end do
    
    cat_param%poros  = scalar
    cat_param%cond   = scalar
    cat_param%psis   = scalar
    cat_param%bee    = scalar
    
    cat_param%wpwet  = scalar
    
    cat_param%gnu    = scalar 
        
    cat_param%vgwmax = scalar
    
    cat_param%vegcls = nint(scalar)
    cat_param%cat_id = nint(scalar)
        
    cat_param%bf1    = scalar
    cat_param%bf2    = scalar
    cat_param%bf3    = scalar
    cat_param%cdcr1  = scalar
    cat_param%cdcr2  = scalar
    cat_param%ars1   = scalar
    cat_param%ars2   = scalar
    cat_param%ars3   = scalar
    cat_param%ara1   = scalar
    cat_param%ara2   = scalar
    cat_param%ara3   = scalar
    cat_param%ara4   = scalar
    cat_param%arw1   = scalar
    cat_param%arw2   = scalar
    cat_param%arw3   = scalar
    cat_param%arw4   = scalar
    cat_param%tsa1   = scalar
    cat_param%tsa2   = scalar
    cat_param%tsb1   = scalar
    cat_param%tsb2   = scalar
    cat_param%atau   = scalar
    cat_param%btau   = scalar
    
  end subroutine scalar2cat_param
    
end module clsmf25_types
  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if 0

! driver routines for testing

program test_catch_types

  ! use module catch_types
  
  implicit none
  
  type(cat_diagn_type) :: cat_diagn_1, cat_diagn_2
  
  cat_diagn_1 = 1.
  cat_diagn_2 = 2.
  
  write (*,*) cat_diagn_1
  write (*,*) cat_diagn_2
  
  cat_diagn_2 = cat_diagn_1 + cat_diagn_2

  write (*,*) cat_diagn_1
  write (*,*) cat_diagn_2

end program test_catch_types

#endif

! ========================== EOF ==================================
