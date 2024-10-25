!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISWRFimport_module
! 
!BOP
! !MODULE: LISWRFimport_module
! 
! !DESCRIPTION: 
!   Defines the import state variables into LIS
! 
! !REVISION HISTORY: 
! 9 Nov 2005: Sujay Kumar, Initial Specification
!
! !USES: 
!EOP
 
  PRIVATE
  
  type, public :: liswrfimport
     logical       :: startflag
     integer       :: yr
     integer       :: mo
     integer       :: da
     integer       :: hr
     integer       :: mn
     integer       :: ss
     real, allocatable :: data_tair(:,:)
     real, allocatable :: data_qair(:,:)
     real, allocatable :: data_swd(:,:)
     real, allocatable :: data_lwd(:,:)
     real, allocatable :: data_ch(:,:)
     real, allocatable :: data_chs2(:,:)
     real, allocatable :: data_cqs2(:,:)
     real, allocatable :: data_psurf(:,:)
     real, allocatable :: data_prcp(:,:)
     real, allocatable :: data_z(:,:)
     real, allocatable :: data_q2sat(:,:)
     real, allocatable :: data_emiss(:,:)
     real, allocatable :: data_u(:,:)
     real, allocatable :: data_v(:,:)
     real, allocatable :: data_cosz(:,:)
     real, allocatable :: data_xice(:,:)
     real, allocatable :: data_tmn(:,:) ! EMK
!ccc - co2 forcing from WRF for CABLE LSM
     real, allocatable :: data_co2(:)
  end type liswrfimport

end module LISWRFimport_module
