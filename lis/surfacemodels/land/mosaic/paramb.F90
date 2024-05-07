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
! !ROUTINE: paramb
!  \label{paramb}
! 
! !REVISION HISTORY: 
!  11 Mar 2001 : Matt Rodell, Initial Specification
!
! !INTERFACE: 
subroutine paramb (ntiles,tex,b)

  implicit none
! !ARGUMENTS: 
  integer :: ntiles
  integer :: tex(ntiles)
  real    :: b(ntiles)

!
! !DESCRIPTION:
!   This subroutine assigns the dimensionless b parameter based on
!   values modified from the soil texture class and the table of 
!   cosby et al. [1984].
!   
!   The arguments are: 
!   \begin{description}
!   \item{ntiles}
!     number of tiles 
!    \item{tex}
!     texture class
!    \item{b}
!     b parameter
!    \end{description}
!EOP

  real, parameter :: bc = 7.62          ! cosby's b values for 
  real, parameter :: bsic = 6.12        ! each texture class, 
  real, parameter :: bsc = 9.19         ! *minus one standard
  real, parameter :: bcl = 4.43         ! deviation*
  real, parameter :: bsicl = 4.39
  real, parameter :: bscl = 3.38
  real, parameter :: bl = 3.59
  real, parameter :: bsil = 3.61
  real, parameter :: bsl = 3.34
  real, parameter :: bls = 2.31
  real, parameter :: bs = 1.41
  
  integer :: i

! assign b values.
  do i=1,ntiles
     select case (tex(i))
     case (1)
        b(i) = bc
     case (2)
        b(i) = bsic
     case (3) 
        b(i) = bsc
     case (4)
        b(i) = bcl
     case (5)
        b(i) = bsicl
     case (6)
        b(i) = bscl
     case (7)
        b(i) = bl
     case (8)
        b(i) = bsil
     case (9) 
        b(i) = bsl
     case (10) 
        b(i) = bls
     case (11)
        b(i) = bs
!    non-land points.
     case default
        b(i) = -0.99
     end select
  end do !i
  return
end subroutine paramb

