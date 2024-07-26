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
! !MODULE: lisgceimport_module.F90
!
! !DESCRIPTION:
!   Defines the import state variables into LIS
!
!
! !REVISION HISTORY:
!
!  09 Nov 2005: Sujay Kumar; Initial specification
!
! !INTERFACE:
module lisgceimport_module
  implicit none
  public lisgceimport
! !ARGUMENTS:  
  type lisgceimport
! Time information
     integer       :: yr
     integer       :: mo
     integer       :: da
     integer       :: hr
     integer       :: mn
     integer       :: ss
     real*8, allocatable :: data_tair(:,:)
     real*8, allocatable :: data_qair(:,:)
     real*8, allocatable :: data_swd(:,:)
     real*8, allocatable :: data_lwd(:,:)
     real*8, allocatable :: data_psurf(:,:)
     real*8, allocatable :: data_prcp(:,:)
     real*8, allocatable :: data_u(:,:)
     real*8, allocatable :: data_v(:,:)
  end type lisgceimport
!EOP
end module lisgceimport_module


