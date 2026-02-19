!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISWRFexport_module
! 
!BOP
! !MODULE: LISWRFexport_module
! 
! !DESCRIPTION: 
!   Defines the export state variables from LIS
! 
! !REVISION HISTORY: 
! 9 Nov 2005: Sujay Kumar, Initial Specification
!
! !USES: 
!EOP
 
  PRIVATE
  
  type, public :: liswrfexport
     real, allocatable :: avgsurft(:,:)
     real, allocatable :: qh(:,:)
     real, allocatable :: eta_kinematic(:,:)
     real, allocatable :: qle(:,:)
     real, allocatable :: qg(:,:)
     real, allocatable :: rootmoist(:,:)
     real, allocatable :: soilm(:,:)
     real, allocatable :: qs(:,:)
     real, allocatable :: qsb(:,:)
     real, allocatable :: albedo(:,:)
     real, allocatable :: znt(:,:)
     real, allocatable :: snocvr(:,:)
     real, allocatable :: q1(:,:)
     real, allocatable :: snow(:,:)
     real, allocatable :: cmc(:,:)
     real, allocatable :: qsm(:,:)
     real, allocatable :: chs2(:,:)
     real, allocatable :: cqs2(:,:)
     real, allocatable :: smc1(:,:)
     real, allocatable :: smc2(:,:)
     real, allocatable :: smc3(:,:)
     real, allocatable :: smc4(:,:)
     real, allocatable :: stc1(:,:)
     real, allocatable :: stc2(:,:)
     real, allocatable :: stc3(:,:)
     real, allocatable :: stc4(:,:)     
     real, allocatable :: sh2o1(:,:)
     real, allocatable :: sh2o2(:,:)
     real, allocatable :: sh2o3(:,:)
     real, allocatable :: sh2o4(:,:)
     real, allocatable :: xland(:,:)
     real, allocatable :: emiss(:,:)
     real, allocatable :: xice(:,:)
     real, allocatable :: lispor(:,:)
     real, allocatable :: snowh(:,:) ! EMK
     real, allocatable :: relsmc1(:,:)
     real, allocatable :: relsmc2(:,:)
     real, allocatable :: relsmc3(:,:)
     real, allocatable :: relsmc4(:,:)
#ifdef WRF_HYDRO
     real, allocatable :: infxsrt(:,:)
     real, allocatable :: soldrain(:,:)
#endif
#ifdef PARFLOW
     real, allocatable :: wtrflx1(:,:)
     real, allocatable :: wtrflx2(:,:)
     real, allocatable :: wtrflx3(:,:)
     real, allocatable :: wtrflx4(:,:)
#endif

     real, allocatable :: avgsurft_t(:)
     real, allocatable :: qh_t(:)
     real, allocatable :: eta_kinematic_t(:)
     real, allocatable :: qle_t(:)
     real, allocatable :: qg_t(:)
     real, allocatable :: rootmoist_t(:)
     real, allocatable :: soilm_t(:)
     real, allocatable :: qs_t(:)
     real, allocatable :: qsb_t(:)
     real, allocatable :: albedo_t(:)
     real, allocatable :: znt_t(:)
     real, allocatable :: snocvr_t(:)
     real, allocatable :: q1_t(:)
     real, allocatable :: snow_t(:)
     real, allocatable :: cmc_t(:)
     real, allocatable :: qsm_t(:)
     real, allocatable :: chs2_t(:)
     real, allocatable :: cqs2_t(:)
     real, allocatable :: smc1_t(:)
     real, allocatable :: smc2_t(:)
     real, allocatable :: smc3_t(:)
     real, allocatable :: smc4_t(:)
     real, allocatable :: stc1_t(:)
     real, allocatable :: stc2_t(:)
     real, allocatable :: stc3_t(:)
     real, allocatable :: stc4_t(:)     
     real, allocatable :: sh2o1_t(:)
     real, allocatable :: sh2o2_t(:)
     real, allocatable :: sh2o3_t(:)
     real, allocatable :: sh2o4_t(:)
     real, allocatable :: xland_t(:)
     real, allocatable :: emiss_t(:)
     real, allocatable :: xice_t(:)
     real, allocatable :: lispor_t(:)
     real, allocatable :: snowh_t(:)
     real, allocatable :: relsmc1_t(:)
     real, allocatable :: relsmc2_t(:)
     real, allocatable :: relsmc3_t(:)
     real, allocatable :: relsmc4_t(:)
#ifdef WRF_HYDRO
     real, allocatable :: infxsrt_t(:)
     real, allocatable :: soldrain_t(:)
#endif
#ifdef PARFLOW
     real, allocatable :: wtrflx1_t(:)
     real, allocatable :: wtrflx2_t(:)
     real, allocatable :: wtrflx3_t(:)
     real, allocatable :: wtrflx4_t(:)
#endif
  end type liswrfexport

end module LISWRFexport_module
