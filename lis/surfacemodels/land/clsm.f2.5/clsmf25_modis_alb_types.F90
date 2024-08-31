!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module clsmf25_modis_alb_types

  ! definition of types and associated operators for Catchment Model
  !
  ! IMPORTANT:
  ! When adding a field to any of the derived types, must also update
  ! the associated assignment and operator definitions.
  ! THERE IS NO WARNING/ERROR IF OPERATOR IS NOT DEFINED FOR ALL FIELDS!
  !
  ! reichle, 21 May 2003
  ! reichle, 25 Jan 2005 - added cat_force_type
  !
  ! --------------------------------------------------------------------------
  
  implicit none
  
  ! everything is private by default unless made public
  
  private
  
  public :: modis_alb_type
  
  public :: assignment (=), operator (/), operator (+)
  
  ! MODIS albedo scaling factors
  
  type :: modis_alb_type

     real :: albvr    !Direct visible or blacksky 0.3-0.7 
     real :: albnr    !Direct infrared or blacksky 0.7-5.0
     real :: albvf    !Diffuse visible or Whitesky 0.3-0.7 
     real :: albnf    !Diffuse infrared or Whitesky 0.7-5.0
     
  end type modis_alb_type
  
  ! ---------------------------------------------------------
  
  ! ----------------------------------------------------------------
  
  interface assignment (=)
     module procedure scalar2cat_modis
  end interface
  
  interface operator (/)
     module procedure cat_div_modis
  end interface

  interface operator (+)
     module procedure add_cat_modis
  end interface
  
contains
  
  subroutine scalar2cat_modis(alb_type, scalar )
    
    implicit none
    
    real, intent(in) :: scalar
    
    type(modis_alb_type), intent(out) :: alb_type
    
    alb_type%albvr    = scalar
    alb_type%albnr    = scalar
    alb_type%albvf    = scalar
    alb_type%albnf    = scalar

  end subroutine scalar2cat_modis
  
  ! -----------------------------------------------------------

  function  cat_div_modis ( alb_type, scalar )
    
    implicit none

    type(modis_alb_type)             :: cat_div_modis
    type(modis_alb_type), intent(in) :: alb_type    

    real, intent(in) :: scalar
    
    cat_div_modis%albvr =  alb_type%albvr /scalar
    cat_div_modis%albnr =  alb_type%albnr /scalar
    cat_div_modis%albvf =  alb_type%albvf /scalar
    cat_div_modis%albnf =  alb_type%albnf /scalar
    
  end function cat_div_modis

  ! -----------------------------------------------------------

  function add_cat_modis (cat_modis_1, cat_modis_2 )
    
    implicit none

    type(modis_alb_type)             :: add_cat_modis
    type(modis_alb_type), intent(in) :: cat_modis_1, cat_modis_2

    add_cat_modis%albvr = cat_modis_1%albvr + cat_modis_2%albvr
    add_cat_modis%albnr = cat_modis_1%albnr + cat_modis_2%albnr
    add_cat_modis%albvf = cat_modis_1%albvf + cat_modis_2%albvf
    add_cat_modis%albnf = cat_modis_1%albnf + cat_modis_2%albnf

  end function add_cat_modis

end module clsmf25_modis_alb_types
! ========================== EOF ==================================
