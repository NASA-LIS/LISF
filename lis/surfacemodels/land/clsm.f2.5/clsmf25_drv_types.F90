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
! !MODULE: clsmf25_drv_types
!
! !DESCRIPTION:
! definition of types and associated operators for Catchment Model driver
!
! IMPORTANT:
! When adding a field to any of the derived types, must also update
! the associated assignment and operator definitions.
! THERE IS NO WARNING/ERROR IF OPERATOR IS NOT DEFINED FOR ALL FIELDS
!
! !REVISION HISTORY:
! reichle, 10 May 2005
! reichle, 10 Jun 2005 - converted met_force_type to ALMA
! 23 Nov 2012: David Mocko, Added RefH to met_force_type
!
! !INTERFACE:
!EOP
module clsmf25_drv_types
  
  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: met_force_type, out_avg_type, out_avg_interval_type
  public :: assignment (=), operator (/), operator (+)
  
  ! ---------------------------------------------------------------------
  !
  ! meteorological forcing variables
  !
  ! The Catchment model driver requires forcing fields of the types
  ! and units in the structure met_force_type.
  ! Make sure to convert the native types and units of each forcing
  ! data set into the types and units of met_force_type right after
  ! the native data have been read.  See for example get_Berg_netcdf().
  
  type met_force_type
     real :: Tair                 ! 2m air temperature            [K]
     real :: Qair                 ! 2m specific humidity          [kg/kg]
     real :: Psurf                ! surface pressure              [Pa]
     real :: Rainf_C              ! convective rainfall           [kg/m2/s]
     real :: Rainf                ! total rainfall                [kg/m2/s]
     real :: Snowf                ! total snowfall                [kg/m2/s]
     real :: LWdown               ! downward longwave radiation   [W/m2]
     real :: SWdown               ! downward shortwave radiation  [W/m2]
     real :: Wind                 ! wind speed                    [m/s]
     real :: RefH                 ! reference height              [m]
     real :: SWnet
     real :: pardr
     real :: pardf
  end type met_force_type
  
  ! --------------------------------------------------------------
  
  interface assignment (=)
     module procedure scalar2met_force
  end interface
  
  interface operator (/)
     module procedure met_force_div_scalar
  end interface
  
  interface operator (+)
     module procedure add_met_force
  end interface
  
  ! ---------------------------------------------------------------
  
  type out_avg_interval_type
     logical :: xhourly
     logical :: daily
     logical :: pentad
     logical :: monthly
  end type out_avg_interval_type
  
  type out_avg_type
     type(out_avg_interval_type) :: tile
     type(out_avg_interval_type) :: grid
     type(out_avg_interval_type) :: plume
  end type out_avg_type
  
contains

  ! --------------------------------------------------

  subroutine scalar2met_force( met_force, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(met_force_type), intent(out) :: met_force
    
    met_force%Tair     = scalar
    met_force%Qair     = scalar
    met_force%Psurf    = scalar
    met_force%Rainf_C  = scalar
    met_force%Rainf    = scalar
    met_force%Snowf    = scalar
    met_force%LWdown   = scalar
    met_force%SWdown   = scalar
    met_force%Wind     = scalar
        
  end subroutine scalar2met_force
  
  ! ---------------------------------------------------
  
  function met_force_div_scalar( met_force, scalar )
    
    implicit none
    
    type(met_force_type)             :: met_force_div_scalar
    type(met_force_type), intent(in) :: met_force
    
    real, intent(in) :: scalar
    
    met_force_div_scalar%Tair     =     met_force%Tair     / scalar
    met_force_div_scalar%Qair     =     met_force%Qair     / scalar
    met_force_div_scalar%Psurf    =     met_force%Psurf    / scalar
    met_force_div_scalar%Rainf_C  =     met_force%Rainf_C  / scalar
    met_force_div_scalar%Rainf    =     met_force%Rainf    / scalar
    met_force_div_scalar%Snowf    =     met_force%Snowf    / scalar
    met_force_div_scalar%LWdown   =     met_force%LWdown   / scalar
    met_force_div_scalar%SWdown   =     met_force%SWdown   / scalar
    met_force_div_scalar%Wind     =     met_force%Wind     / scalar
        
  end function met_force_div_scalar

  ! -----------------------------------------------------------

  function add_met_force( met_force_1, met_force_2 )
    
    implicit none

    type(met_force_type)             :: add_met_force
    type(met_force_type), intent(in) :: met_force_1, met_force_2
    
    add_met_force%Tair     = met_force_1%Tair     + met_force_2%Tair    
    add_met_force%Qair     = met_force_1%Qair     + met_force_2%Qair     
    add_met_force%Psurf    = met_force_1%Psurf    + met_force_2%Psurf     
    add_met_force%Rainf_C  = met_force_1%Rainf_C  + met_force_2%Rainf_C  
    add_met_force%Rainf    = met_force_1%Rainf    + met_force_2%Rainf   
    add_met_force%Snowf    = met_force_1%Snowf    + met_force_2%Snowf    
    add_met_force%LWdown   = met_force_1%LWdown   + met_force_2%LWdown   
    add_met_force%SWdown   = met_force_1%SWdown   + met_force_2%SWdown   
    add_met_force%Wind     = met_force_1%Wind     + met_force_2%Wind    
    
  end function add_met_force

  
end module clsmf25_drv_types

