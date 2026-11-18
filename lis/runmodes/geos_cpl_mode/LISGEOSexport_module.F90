!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISGEOSexport_module
! 
!BOP
! !MODULE: LISGEOSexport_module
! 
! !DESCRIPTION: 
!   Defines the export state variables from LIS
! 
! !REVISION HISTORY: 
! 9 Nov 2026: Ronnie Abolafia-Rosenzweig, Initial Specification
!
! !USES: 
!EOP
 
  PRIVATE
  
  type, public :: lisgeosexport
     ! --- Standard/Original LIS 2D Export Arrays ---
     real, allocatable :: avgsurft(:,:)
     real, allocatable :: qh(:,:)
     real, allocatable :: qle(:,:)
     real, allocatable :: rootmoist(:,:)
     real, allocatable :: soilm(:,:)
     real, allocatable :: qs(:,:)
     real, allocatable :: qsb(:,:)
     real, allocatable :: albedo(:,:)
     real, allocatable :: snocvr(:,:)
     real, allocatable :: snow(:,:)
     real, allocatable :: smc1(:,:)
     real, allocatable :: smc2(:,:)
     real, allocatable :: smc3(:,:)
     real, allocatable :: smc4(:,:)
     real, allocatable :: stc1(:,:)
     real, allocatable :: stc2(:,:)
     real, allocatable :: stc3(:,:)
     real, allocatable :: stc4(:,:)
     real, allocatable :: snowh(:,:)
     real, allocatable :: emiss(:,:)
     real, allocatable :: sh2o1(:,:)
     real, allocatable :: sh2o2(:,:)
     real, allocatable :: sh2o3(:,:)
     real, allocatable :: sh2o4(:,:)

     ! --- New GEOS Coupling 2D Export Arrays ---
     real, allocatable :: rhos(:,:)
     real, allocatable :: totalrunoff(:,:)
     real, allocatable :: ch_f(:,:)
     real, allocatable :: cm_f(:,:)
     real, allocatable :: avgsurfq(:,:)
     real, allocatable :: snowdepth_sca(:,:)
     real, allocatable :: albedovisdir(:,:)
     real, allocatable :: albedonirdir(:,:)
     real, allocatable :: albedovisdif(:,:)
     real, allocatable :: albedonirdif(:,:)  
     real, allocatable :: z0(:,:)
     real, allocatable :: zpd(:,:)
     real, allocatable :: ustarveg(:,:)
     real, allocatable :: ustarbare(:,:)
     real, allocatable :: ustar(:,:)           
     real, allocatable :: taux(:,:)
     real, allocatable :: tauy(:,:)
     real, allocatable :: rtmass(:,:)
     real, allocatable :: stmass(:,:)
     real, allocatable :: wood(:,:)
     real, allocatable :: stblcp(:,:)
     real, allocatable :: fastcp(:,:)
     real, allocatable :: lai(:,:)
     real, allocatable :: sai(:,:)
     real, allocatable :: vai(:,:)
     real, allocatable :: evap(:,:)
     real, allocatable :: greenness(:,:)
     real, allocatable :: landcover(:,:)
     real, allocatable :: nee(:,:)
     real, allocatable :: canopint(:,:)
     real, allocatable :: tair_t2m(:,:)

     ! --- Standard/Original LIS 1D Tile Arrays ---
     real, allocatable :: avgsurft_t(:)
     real, allocatable :: qh_t(:)
     real, allocatable :: qle_t(:)
     real, allocatable :: rootmoist_t(:)
     real, allocatable :: soilm_t(:)
     real, allocatable :: qs_t(:)
     real, allocatable :: qsb_t(:)
     real, allocatable :: albedo_t(:)
     real, allocatable :: snocvr_t(:)
     real, allocatable :: snow_t(:)
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
     real, allocatable :: emiss_t(:)
     real, allocatable :: snowh_t(:)
     
    ! --- New GEOS Coupling 1D Tile Arrays ---
     real, allocatable :: rhos_t(:)
     real, allocatable :: totalrunoff_t(:)
     real, allocatable :: ch_f_t(:)
     real, allocatable :: cm_f_t(:)
     real, allocatable :: avgsurfq_t(:)
     real, allocatable :: snowdepth_sca_t(:)
     real, allocatable :: albedovisdir_t(:)
     real, allocatable :: albedonirdir_t(:)
     real, allocatable :: albedovisdif_t(:)
     real, allocatable :: albedonirdif_t(:)
     real, allocatable :: z0_t(:)
     real, allocatable :: zpd_t(:)
     real, allocatable :: ustarveg_t(:)
     real, allocatable :: ustarbare_t(:)
     real, allocatable :: ustar_t(:)          
     real, allocatable :: taux_t(:)
     real, allocatable :: tauy_t(:)
     real, allocatable :: rtmass_t(:)
     real, allocatable :: stmass_t(:)
     real, allocatable :: wood_t(:)
     real, allocatable :: stblcp_t(:)
     real, allocatable :: fastcp_t(:)
     real, allocatable :: lai_t(:)
     real, allocatable :: sai_t(:)
     real, allocatable :: vai_t(:)
     real, allocatable :: evap_t(:)
     real, allocatable :: greenness_t(:)
     real, allocatable :: landcover_t(:)
     real, allocatable :: nee_t(:)
     real, allocatable :: canopint_t(:)
     real, allocatable :: tair_t2m_t(:)
#ifdef WRF_HYDRO
     real, allocatable :: infxsrt(:,:)
     real, allocatable :: soldrain(:,:)
     real, allocatable :: infxsrt_t(:)
     real, allocatable :: soldrain_t(:)
#endif
#ifdef PARFLOW
     real, allocatable :: wtrflx1(:,:)
     real, allocatable :: wtrflx2(:,:)
     real, allocatable :: wtrflx3(:,:)
     real, allocatable :: wtrflx4(:,:)
     real, allocatable :: wtrflx1_t(:)
     real, allocatable :: wtrflx2_t(:)
     real, allocatable :: wtrflx3_t(:)
     real, allocatable :: wtrflx4_t(:)
#endif
  end type lisgeosexport

end module LISGEOSexport_module