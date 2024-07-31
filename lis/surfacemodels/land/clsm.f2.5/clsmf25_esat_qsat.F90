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
! !MODULE: clsmf25_esat_qsat
!
! !DESCRIPTION:
!
! cleaned up qsat.f and function esat() from process\_new.f
!
! !REVISION HISTORY:
! reichle, 13 Jun 2005
!
! !INTERFACE:
module clsmf25_esat_qsat
!EOP
  use clsmf25_constants  
  
  implicit none
    
  private

  public :: esat, qsat, dxsat, e2q, q2e
  
  ! ----------------------------------------
  
  real, parameter :: C1 = 3.797915
  real, parameter :: C2 = 7.93252E-6
  real, parameter :: C3 = 2.166847E-3
  
  real, parameter :: C4 = 21.18123
  real, parameter :: C5 = 5418.
  
  real, parameter :: C1LOG = 1.33445         ! = log(C1)
  
  ! the following must be made consistent with Catchment model parameters
  ! if there is ever a module of parameters...
  
!sm  real, parameter :: ALHE = 2.4548E6
!sm  real, parameter :: ALHS = 2.8368E6
  
!sm  real, parameter :: epsilon = 18.01/28.97
  
contains
  
  ! *************************************************************
  
  real function esat(t,ALHX)
    
    implicit none
    
    real, intent(in)           :: t        ! temperature in K
    real, intent(in), optional :: ALHX
    
    if (present(ALHX)) then
       
       esat = C1*EXP((ALHX/ALHE)*(C4-C1LOG-C5/T)) / epsilon
       
    else
       
       esat = exp(C4 - C5/t) / epsilon
       
    end if
    
  end function esat
  
  ! *************************************************************
  
  real function e2q( e, P )
    
    implicit none
    
    real, intent(in) :: e, P       ! pressure in *mbar* !!!!
    
    e2q = e*epsilon/P
    
  end function e2q
  
  ! *************************************************************
  
  real function q2e( q, P )
    
    implicit none
    
    real, intent(in) :: q, P      ! pressure in *mbar* !!!!
    
    q2e = q/epsilon*P
    
  end function q2e
  
  ! *************************************************************
  
  real function qsat(T,P,ALHX)
    
    implicit none
    
    real, intent(in)           :: T        ! temperature in K
    real, intent(in)           :: P        ! pressure in *mbar* !!!!
    real, intent(in), optional :: ALHX
    
    if (present(ALHX)) then
       
       qsat = epsilon/P*esat(T,ALHX)
       
    else
       qsat = epsilon/P*esat(T)
       
    end if
    
  end function qsat
  
  ! *************************************************************
  
  real function dxsat(xsat,t,ALHX)
    
    implicit none
    
    real, intent(in)           :: xsat
    real, intent(in)           :: t       ! temperature in K
    real, intent(in), optional :: ALHX
    
    if (present(ALHX)) then
       
       dxsat = xsat * ALHX/ALHE * C5 / ( T * T )
       
    else
       
       dxsat = xsat * C5 / (T * T )
       
    end if
    
  end function dxsat
  
  ! --------------
  
end module clsmf25_esat_qsat

