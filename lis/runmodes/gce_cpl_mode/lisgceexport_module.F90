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
! !MODULE: lisgceexport_module.F90
!
! !DESCRIPTION:
!   Defines the export state variables into LIS
!   Variable names conform to the ALMA convention
!
! !REVISION HISTORY:
!
!  09 Nov 2005: Sujay Kumar; Initial specification
!
! !INTERFACE:
module lisgceexport_module
  implicit none
  public lisgceexport
! !ARGUMENTS:  
  type lisgceexport
     real*8, allocatable :: swnet(:,:)
     real*8, allocatable :: lwnet(:,:)
     real*8, allocatable :: qle(:,:)
     real*8, allocatable :: qh(:,:)
     real*8, allocatable :: qg(:,:)
     real*8, allocatable :: avgsurft(:,:)
     real*8, allocatable :: albedo(:,:)
     real*8, allocatable :: tauu(:,:)
     real*8, allocatable :: tauv(:,:)
     real*8, allocatable :: tau(:,:)
  end type lisgceexport
!EOP
end module lisgceexport_module


