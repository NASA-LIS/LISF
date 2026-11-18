!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISGEOSGridCompMod
! 
!BOP
! !MODULE: LISGEOSGridCompMod
! 
! !DESCRIPTION: 
!   This module defines the handling of the fields being exchanged between
!   LIS and GEOS. 
! 
! !REVISION HISTORY: 
! 9 Nov 2026: Ronnie Abolafia-Rosenzweig, Initial Specification
!
! !USES: 
  use LISGEOSimport_module
  use LISGEOSexport_module
!EOP

  PRIVATE

  public :: LISGEOS_alloc_states
  public :: LISGEOS_reset_states
  public :: LISGEOS_dealloc_states

  public :: LISGEOS_import
  public :: LISGEOS_export

  
  type(LISGEOSimport), allocatable, target :: LISGEOS_import(:)
  type(LISGEOSexport), allocatable, target :: LISGEOS_export(:)

  contains 
    subroutine LISGEOS_alloc_states

      use LIS_coreMod, only : LIS_rc
      
      implicit none
      
      integer :: n

      allocate(LISGEOS_import(LIS_rc%nnest))
      allocate(LISGEOS_export(LIS_rc%nnest))
      
      do n=1, LIS_rc%nnest
         ! Import states      
         allocate(LISGEOS_import(n)%data_tair(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_qair(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_swd(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_lwd(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_ch(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_psurf(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_prcp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_z(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_q2sat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_emiss(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_u(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_v(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_cosz(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_import(n)%data_tmn(LIS_rc%lnc(n),LIS_rc%lnr(n))) 
         allocate(LISGEOS_import(n)%data_co2(1))
         LISGEOS_import(n)%startflag = .true.

         ! Export states - 1D Tile Arrays
         allocate(LISGEOS_export(n)%avgsurft_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%qh_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%qle_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%albedo_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%smc1_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%smc2_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%smc3_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%smc4_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%stc1_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%stc2_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%stc3_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%stc4_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%snocvr_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%snow_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%snowh_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%rootmoist_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%soilm_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%qs_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%qsb_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%emiss_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%sh2o1_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%sh2o2_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%sh2o3_t(LIS_rc%ntiles(n)))

         ! Export states - New GEOS Coupling 1D Tile Arrays
         allocate(LISGEOS_export(n)%rhos_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%totalrunoff_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%ch_f_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%qair_t2m_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%snowdepth_sca_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%albedovisdir_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%albedonirdir_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%albedovisdif_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%albedonirdif_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%z0_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%zpd_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%ustarveg_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%ustarbare_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%taux_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%tauy_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%rtmass_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%stmass_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%wood_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%stblcp_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%fastcp_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%lai_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%sai_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%evap_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%greenness_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%landcover_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%nee_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%canopint_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%tair_t2m_t(LIS_rc%ntiles(n)))

         ! Export states - 2D Grid Arrays
         allocate(LISGEOS_export(n)%avgsurft(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%qh(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%qle(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%albedo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%smc1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%smc2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%smc3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%smc4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%stc1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%stc2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%stc3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%stc4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%snocvr(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%snow(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%snowh(LIS_rc%lnc(n),LIS_rc%lnr(n))) 
         allocate(LISGEOS_export(n)%rootmoist(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%soilm(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%qs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%qsb(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%emiss(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%sh2o1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%sh2o2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%sh2o3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%sh2o4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         
         ! Export states - New GEOS Coupling 2D Grid Arrays
         allocate(LISGEOS_export(n)%rhos(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%totalrunoff(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%ch_f(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%qair_t2m(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%snowdepth_sca(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%albedovisdir(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%albedonirdir(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%albedovisdif(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%albedonirdif(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%z0(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%zpd(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%ustarveg(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%ustarbare(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%taux(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%tauy(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%rtmass(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%stmass(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%wood(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%stblcp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%fastcp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%lai(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%sai(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%evap(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%greenness(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%landcover(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%nee(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%canopint(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%tair_t2m(LIS_rc%lnc(n),LIS_rc%lnr(n)))

#ifdef GEOS_HYDRO
         allocate(LISGEOS_export(n)%infxsrt_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%soldrain_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%infxsrt(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%soldrain(LIS_rc%lnc(n),LIS_rc%lnr(n)))
#endif
#ifdef PARFLOW
         allocate(LISGEOS_export(n)%wtrflx1_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%wtrflx2_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%wtrflx3_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%wtrflx4_t(LIS_rc%ntiles(n)))
         allocate(LISGEOS_export(n)%wtrflx1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%wtrflx2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%wtrflx3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISGEOS_export(n)%wtrflx4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
#endif
      enddo
    end subroutine LISGEOS_alloc_states

    subroutine LISGEOS_reset_states

      use LIS_coreMod, only : LIS_rc
      
      implicit none
      
      integer :: n

      do n=1, LIS_rc%nnest
         ! Import states
         LISGEOS_import(n)%data_tair = LIS_rc%udef
         LISGEOS_import(n)%data_qair = LIS_rc%udef
         LISGEOS_import(n)%data_swd = LIS_rc%udef
         LISGEOS_import(n)%data_lwd = LIS_rc%udef
         LISGEOS_import(n)%data_ch = LIS_rc%udef
         LISGEOS_import(n)%data_psurf = LIS_rc%udef
         LISGEOS_import(n)%data_prcp = LIS_rc%udef
         LISGEOS_import(n)%data_z = LIS_rc%udef
         LISGEOS_import(n)%data_q2sat = LIS_rc%udef
         LISGEOS_import(n)%data_emiss = LIS_rc%udef
         LISGEOS_import(n)%data_u = LIS_rc%udef
         LISGEOS_import(n)%data_v = LIS_rc%udef
         LISGEOS_import(n)%data_cosz = LIS_rc%udef
         LISGEOS_import(n)%data_tmn = LIS_rc%udef
         
         ! Export states - 1D Tile Arrays
         LISGEOS_export(n)%avgsurft_t = LIS_rc%udef
         LISGEOS_export(n)%qh_t = LIS_rc%udef
         LISGEOS_export(n)%qle_t = LIS_rc%udef
         LISGEOS_export(n)%albedo_t = LIS_rc%udef
         LISGEOS_export(n)%snocvr_t = LIS_rc%udef
         LISGEOS_export(n)%snow_t = LIS_rc%udef
         LISGEOS_export(n)%smc1_t = LIS_rc%udef
         LISGEOS_export(n)%smc2_t = LIS_rc%udef
         LISGEOS_export(n)%smc3_t = LIS_rc%udef
         LISGEOS_export(n)%smc4_t = LIS_rc%udef
         LISGEOS_export(n)%stc1_t = LIS_rc%udef
         LISGEOS_export(n)%stc2_t = LIS_rc%udef
         LISGEOS_export(n)%stc3_t = LIS_rc%udef
         LISGEOS_export(n)%stc4_t = LIS_rc%udef
         LISGEOS_export(n)%rootmoist_t = LIS_rc%udef
         LISGEOS_export(n)%soilm_t = LIS_rc%udef
         LISGEOS_export(n)%qs_t = LIS_rc%udef
         LISGEOS_export(n)%qsb_t = LIS_rc%udef
         LISGEOS_export(n)%snowh_t = LIS_rc%udef
         LISGEOS_export(n)%emiss_t = LIS_rc%udef
         LISGEOS_export(n)%sh2o1_t = LIS_rc%udef
         LISGEOS_export(n)%sh2o2_t = LIS_rc%udef
         LISGEOS_export(n)%sh2o3_t = LIS_rc%udef
         LISGEOS_export(n)%sh2o4_t = LIS_rc%udef

         ! Export states - New GEOS Coupling 1D Tile Arrays
         LISGEOS_export(n)%rhos_t = LIS_rc%udef
         LISGEOS_export(n)%totalrunoff_t = LIS_rc%udef
         LISGEOS_export(n)%ch_f_t = LIS_rc%udef
         LISGEOS_export(n)%qair_t2m_t = LIS_rc%udef
         LISGEOS_export(n)%snowdepth_sca_t = LIS_rc%udef
         LISGEOS_export(n)%albedovisdir_t = LIS_rc%udef
         LISGEOS_export(n)%albedonirdir_t = LIS_rc%udef
         LISGEOS_export(n)%albedovisdif_t = LIS_rc%udef
         LISGEOS_export(n)%albedonirdif_t = LIS_rc%udef
         LISGEOS_export(n)%z0_t = LIS_rc%udef
         LISGEOS_export(n)%zpd_t = LIS_rc%udef
         LISGEOS_export(n)%ustarveg_t = LIS_rc%udef
         LISGEOS_export(n)%ustarbare_t = LIS_rc%udef
         LISGEOS_export(n)%taux_t = LIS_rc%udef
         LISGEOS_export(n)%tauy_t = LIS_rc%udef
         LISGEOS_export(n)%rtmass_t = LIS_rc%udef
         LISGEOS_export(n)%stmass_t = LIS_rc%udef
         LISGEOS_export(n)%wood_t = LIS_rc%udef
         LISGEOS_export(n)%stblcp_t = LIS_rc%udef
         LISGEOS_export(n)%fastcp_t = LIS_rc%udef
         LISGEOS_export(n)%lai_t = LIS_rc%udef
         LISGEOS_export(n)%sai_t = LIS_rc%udef
         LISGEOS_export(n)%evap_t = LIS_rc%udef
         LISGEOS_export(n)%greenness_t = LIS_rc%udef
         LISGEOS_export(n)%landcover_t = LIS_rc%udef
         LISGEOS_export(n)%nee_t = LIS_rc%udef
         LISGEOS_export(n)%canopint_t = LIS_rc%udef
         LISGEOS_export(n)%tair_t2m_t = LIS_rc%udef

         ! Export states - 2D Grid Arrays
         LISGEOS_export(n)%avgsurft = LIS_rc%udef
         LISGEOS_export(n)%qh = LIS_rc%udef
         LISGEOS_export(n)%qle = LIS_rc%udef
         LISGEOS_export(n)%albedo = LIS_rc%udef
         LISGEOS_export(n)%smc1 = LIS_rc%udef
         LISGEOS_export(n)%smc2 = LIS_rc%udef
         LISGEOS_export(n)%smc3 = LIS_rc%udef
         LISGEOS_export(n)%smc4 = LIS_rc%udef
         LISGEOS_export(n)%stc1 = LIS_rc%udef
         LISGEOS_export(n)%stc2 = LIS_rc%udef
         LISGEOS_export(n)%stc3 = LIS_rc%udef
         LISGEOS_export(n)%stc4 = LIS_rc%udef
         LISGEOS_export(n)%snocvr = LIS_rc%udef
         LISGEOS_export(n)%snow = LIS_rc%udef
         LISGEOS_export(n)%snowh = LIS_rc%udef
         LISGEOS_export(n)%rootmoist = LIS_rc%udef
         LISGEOS_export(n)%soilm = LIS_rc%udef
         LISGEOS_export(n)%qs = LIS_rc%udef
         LISGEOS_export(n)%qsb = LIS_rc%udef
         LISGEOS_export(n)%emiss = LIS_rc%udef
         LISGEOS_export(n)%sh2o1 = LIS_rc%udef
         LISGEOS_export(n)%sh2o2 = LIS_rc%udef
         LISGEOS_export(n)%sh2o3 = LIS_rc%udef
         LISGEOS_export(n)%sh2o4 = LIS_rc%udef
        
         ! Export states - New GEOS Coupling 2D Grid Arrays
         LISGEOS_export(n)%rhos = LIS_rc%udef
         LISGEOS_export(n)%totalrunoff = LIS_rc%udef
         LISGEOS_export(n)%ch_f = LIS_rc%udef
         LISGEOS_export(n)%qair_t2m = LIS_rc%udef
         LISGEOS_export(n)%snowdepth_sca = LIS_rc%udef
         LISGEOS_export(n)%albedovisdir = LIS_rc%udef
         LISGEOS_export(n)%albedonirdir = LIS_rc%udef
         LISGEOS_export(n)%albedovisdif = LIS_rc%udef
         LISGEOS_export(n)%albedonirdif = LIS_rc%udef
         LISGEOS_export(n)%z0 = LIS_rc%udef
         LISGEOS_export(n)%zpd = LIS_rc%udef
         LISGEOS_export(n)%ustarveg = LIS_rc%udef
         LISGEOS_export(n)%ustarbare = LIS_rc%udef
         LISGEOS_export(n)%taux = LIS_rc%udef
         LISGEOS_export(n)%tauy = LIS_rc%udef
         LISGEOS_export(n)%rtmass = LIS_rc%udef
         LISGEOS_export(n)%stmass = LIS_rc%udef
         LISGEOS_export(n)%wood = LIS_rc%udef
         LISGEOS_export(n)%stblcp = LIS_rc%udef
         LISGEOS_export(n)%fastcp = LIS_rc%udef
         LISGEOS_export(n)%lai = LIS_rc%udef
         LISGEOS_export(n)%sai = LIS_rc%udef
         LISGEOS_export(n)%evap = LIS_rc%udef
         LISGEOS_export(n)%greenness = LIS_rc%udef
         LISGEOS_export(n)%landcover = LIS_rc%udef
         LISGEOS_export(n)%nee = LIS_rc%udef
         LISGEOS_export(n)%canopint = LIS_rc%udef
         LISGEOS_export(n)%tair_t2m = LIS_rc%udef

#ifdef GEOS_HYDRO
         LISGEOS_export(n)%infxsrt_t = LIS_rc%udef
         LISGEOS_export(n)%soldrain_t = LIS_rc%udef
         LISGEOS_export(n)%infxsrt = LIS_rc%udef
         LISGEOS_export(n)%soldrain = LIS_rc%udef
#endif
#ifdef PARFLOW
         LISGEOS_export(n)%wtrflx1_t = LIS_rc%udef
         LISGEOS_export(n)%wtrflx2_t = LIS_rc%udef
         LISGEOS_export(n)%wtrflx3_t = LIS_rc%udef
         LISGEOS_export(n)%wtrflx4_t = LIS_rc%udef
         LISGEOS_export(n)%wtrflx1 = LIS_rc%udef
         LISGEOS_export(n)%wtrflx2 = LIS_rc%udef
         LISGEOS_export(n)%wtrflx3 = LIS_rc%udef
         LISGEOS_export(n)%wtrflx4 = LIS_rc%udef
#endif
      enddo
    end subroutine LISGEOS_reset_states

    subroutine LISGEOS_dealloc_states
      use LIS_coreMod, only : LIS_rc
      implicit none
      integer :: n

      if (allocated(LISGEOS_import)) then
         deallocate(LISGEOS_import)
      endif

      if (allocated(LISGEOS_export)) then
         deallocate(LISGEOS_export)
      endif
    end subroutine LISGEOS_dealloc_states

end module LISGEOSGridCompMod