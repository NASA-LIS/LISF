!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISWRFGridCompMod
! 
!BOP
! !MODULE: LISWRFGridCompMod
! 
! !DESCRIPTION: 
!   This module defines the handling of the fields being exchanged between
!   LIS and WRF. 
! 
! !REVISION HISTORY: 
! 9 Nov 2008: Sujay Kumar, Initial Specification
!
! !USES: 
  use LISWRFimport_module
  use LISWRFexport_module
!EOP

  PRIVATE

  public :: LISWRF_alloc_states
  public :: LISWRF_reset_states

  public :: LISWRF_import
  public :: LISWRF_export

  
  type(LISWRFimport), allocatable, target :: LISWRF_import(:)
  type(LISWRFexport), allocatable, target :: LISWRF_export(:)

  contains 
    subroutine LISWRF_alloc_states

      use LIS_coreMod, only : LIS_rc
      
      implicit none
      
      integer :: n 

      allocate(LISWRF_import(LIS_rc%nnest))
      allocate(LISWRF_export(LIS_rc%nnest))
      
      do n=1, LIS_rc%nnest
         ! Import states      
         allocate(LISWRF_import(n)%data_tair(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_qair(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_swd(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_lwd(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_ch(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_chs2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_cqs2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_psurf(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_prcp(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_z(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_q2sat(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_emiss(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_u(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_v(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_cosz(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_xice(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_import(n)%data_tmn(LIS_rc%lnc(n),LIS_rc%lnr(n))) ! EMK
!ccc - co2 forcing from WRF for CABLE LSM
         allocate(LISWRF_import(n)%data_co2(1))
         LISWRF_import(n)%startflag = .true. 

!Export states
         allocate(LISWRF_export(n)%avgsurft(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%qh(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%eta_kinematic(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%qle(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%qg(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%albedo(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%znt(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%q1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%smc1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%smc2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%smc3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%smc4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%stc1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%stc2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%stc3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%stc4(LIS_rc%lnc(n),LIS_rc%lnr(n)))     
         allocate(LISWRF_export(n)%chs2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%cqs2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%snocvr(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%snow(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%snowh(LIS_rc%lnc(n),LIS_rc%lnr(n))) ! EMK
         allocate(LISWRF_export(n)%lispor(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%rootmoist(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%soilm(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%qs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%qsb(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%cmc(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%qsm(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%emiss(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%xice(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%sh2o1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%sh2o2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%sh2o3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%sh2o4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%relsmc1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%relsmc2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%relsmc3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%relsmc4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%xland(LIS_rc%lnc(n),LIS_rc%lnr(n)))
#ifdef WRF_HYDRO
         allocate(LISWRF_export(n)%infxsrt(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%soldrain(LIS_rc%lnc(n),LIS_rc%lnr(n)))
#endif
#ifdef PARFLOW
         allocate(LISWRF_export(n)%wtrflx1(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%wtrflx2(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%wtrflx3(LIS_rc%lnc(n),LIS_rc%lnr(n)))
         allocate(LISWRF_export(n)%wtrflx4(LIS_rc%lnc(n),LIS_rc%lnr(n)))
#endif

!Export states
         allocate(LISWRF_export(n)%avgsurft_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%qh_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%eta_kinematic_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%qle_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%qg_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%albedo_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%znt_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%q1_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%smc1_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%smc2_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%smc3_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%smc4_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%stc1_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%stc2_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%stc3_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%stc4_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%chs2_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%cqs2_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%snocvr_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%snow_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%snowh_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%lispor_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%rootmoist_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%soilm_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%qs_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%qsb_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%cmc_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%qsm_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%emiss_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%xice_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%sh2o1_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%sh2o2_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%sh2o3_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%sh2o4_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%relsmc1_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%relsmc2_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%relsmc3_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%relsmc4_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%xland_t(LIS_rc%ntiles(n)))
#ifdef WRF_HYDRO
         allocate(LISWRF_export(n)%infxsrt_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%soldrain_t(LIS_rc%ntiles(n)))
#endif
#ifdef PARFLOW
         allocate(LISWRF_export(n)%wtrflx1_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%wtrflx2_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%wtrflx3_t(LIS_rc%ntiles(n)))
         allocate(LISWRF_export(n)%wtrflx4_t(LIS_rc%ntiles(n)))
#endif
      enddo
    end subroutine LISWRF_alloc_states

    subroutine LISWRF_reset_states()

      use LIS_coreMod, only : LIS_rc

      implicit none

      integer :: n

      do n=1, LIS_rc%nnest
         ! Import states
         LISWRF_import(n)%data_tair = LIS_rc%udef
         LISWRF_import(n)%data_qair = LIS_rc%udef
         LISWRF_import(n)%data_swd = LIS_rc%udef
         LISWRF_import(n)%data_lwd = LIS_rc%udef
         LISWRF_import(n)%data_ch = LIS_rc%udef
         LISWRF_import(n)%data_chs2 = LIS_rc%udef
         LISWRF_import(n)%data_cqs2 = LIS_rc%udef
         LISWRF_import(n)%data_psurf = LIS_rc%udef
         LISWRF_import(n)%data_prcp = LIS_rc%udef
         LISWRF_import(n)%data_z = LIS_rc%udef
         LISWRF_import(n)%data_q2sat = LIS_rc%udef
         LISWRF_import(n)%data_emiss = LIS_rc%udef
         LISWRF_import(n)%data_u = LIS_rc%udef
         LISWRF_import(n)%data_v = LIS_rc%udef
         LISWRF_import(n)%data_cosz = LIS_rc%udef
         LISWRF_import(n)%data_xice = LIS_rc%udef
         LISWRF_import(n)%data_tmn = LIS_rc%udef
         LISWRF_import(n)%data_co2(1) = LIS_rc%udef

!Export states
         LISWRF_export(n)%avgsurft = LIS_rc%udef
         LISWRF_export(n)%qh = LIS_rc%udef
         LISWRF_export(n)%eta_kinematic = LIS_rc%udef
         LISWRF_export(n)%qle = LIS_rc%udef
         LISWRF_export(n)%qg = LIS_rc%udef
         LISWRF_export(n)%albedo = LIS_rc%udef
         LISWRF_export(n)%znt = LIS_rc%udef
         LISWRF_export(n)%q1 = LIS_rc%udef
         LISWRF_export(n)%smc1 = LIS_rc%udef
         LISWRF_export(n)%smc2 = LIS_rc%udef
         LISWRF_export(n)%smc3 = LIS_rc%udef
         LISWRF_export(n)%smc4 = LIS_rc%udef
         LISWRF_export(n)%stc1 = LIS_rc%udef
         LISWRF_export(n)%stc2 = LIS_rc%udef
         LISWRF_export(n)%stc3 = LIS_rc%udef
         LISWRF_export(n)%stc4 = LIS_rc%udef
         LISWRF_export(n)%chs2 = LIS_rc%udef
         LISWRF_export(n)%cqs2 = LIS_rc%udef
         LISWRF_export(n)%snocvr = LIS_rc%udef
         LISWRF_export(n)%snow = LIS_rc%udef
         LISWRF_export(n)%snowh = LIS_rc%udef
         LISWRF_export(n)%lispor = LIS_rc%udef
         LISWRF_export(n)%rootmoist = LIS_rc%udef
         LISWRF_export(n)%soilm = LIS_rc%udef
         LISWRF_export(n)%qs = LIS_rc%udef
         LISWRF_export(n)%qsb = LIS_rc%udef
         LISWRF_export(n)%cmc = LIS_rc%udef
         LISWRF_export(n)%qsm = LIS_rc%udef
         LISWRF_export(n)%emiss = LIS_rc%udef
         LISWRF_export(n)%xice = LIS_rc%udef
         LISWRF_export(n)%sh2o1 = LIS_rc%udef
         LISWRF_export(n)%sh2o2 = LIS_rc%udef
         LISWRF_export(n)%sh2o3 = LIS_rc%udef
         LISWRF_export(n)%sh2o4 = LIS_rc%udef
         LISWRF_export(n)%relsmc1 = LIS_rc%udef
         LISWRF_export(n)%relsmc2 = LIS_rc%udef
         LISWRF_export(n)%relsmc3 = LIS_rc%udef
         LISWRF_export(n)%relsmc4 = LIS_rc%udef
         LISWRF_export(n)%xland = LIS_rc%udef
#ifdef WRF_HYDRO
         LISWRF_export(n)%infxsrt = LIS_rc%udef
         LISWRF_export(n)%soldrain = LIS_rc%udef
#endif
#ifdef PARFLOW
         LISWRF_export(n)%wtrflx1 = LIS_rc%udef
         LISWRF_export(n)%wtrflx2 = LIS_rc%udef
         LISWRF_export(n)%wtrflx3 = LIS_rc%udef
         LISWRF_export(n)%wtrflx4 = LIS_rc%udef
#endif

!Export states
         LISWRF_export(n)%avgsurft_t = LIS_rc%udef
         LISWRF_export(n)%qh_t = LIS_rc%udef
         LISWRF_export(n)%eta_kinematic_t = LIS_rc%udef
         LISWRF_export(n)%qle_t = LIS_rc%udef
         LISWRF_export(n)%qg_t = LIS_rc%udef
         LISWRF_export(n)%albedo_t = LIS_rc%udef
         LISWRF_export(n)%znt_t = LIS_rc%udef
         LISWRF_export(n)%q1_t = LIS_rc%udef
         LISWRF_export(n)%smc1_t = LIS_rc%udef
         LISWRF_export(n)%smc2_t = LIS_rc%udef
         LISWRF_export(n)%smc3_t = LIS_rc%udef
         LISWRF_export(n)%smc4_t = LIS_rc%udef
         LISWRF_export(n)%stc1_t = LIS_rc%udef
         LISWRF_export(n)%stc2_t = LIS_rc%udef
         LISWRF_export(n)%stc3_t = LIS_rc%udef
         LISWRF_export(n)%stc4_t = LIS_rc%udef
         LISWRF_export(n)%chs2_t = LIS_rc%udef
         LISWRF_export(n)%cqs2_t = LIS_rc%udef
         LISWRF_export(n)%snocvr_t = LIS_rc%udef
         LISWRF_export(n)%snow_t = LIS_rc%udef
         LISWRF_export(n)%snowh_t = LIS_rc%udef
         LISWRF_export(n)%lispor_t = LIS_rc%udef
         LISWRF_export(n)%rootmoist_t = LIS_rc%udef
         LISWRF_export(n)%soilm_t = LIS_rc%udef
         LISWRF_export(n)%qs_t = LIS_rc%udef
         LISWRF_export(n)%qsb_t = LIS_rc%udef
         LISWRF_export(n)%cmc_t = LIS_rc%udef
         LISWRF_export(n)%qsm_t = LIS_rc%udef
         LISWRF_export(n)%emiss_t = LIS_rc%udef
         LISWRF_export(n)%xice_t = LIS_rc%udef
         LISWRF_export(n)%sh2o1_t = LIS_rc%udef
         LISWRF_export(n)%sh2o2_t = LIS_rc%udef
         LISWRF_export(n)%sh2o3_t = LIS_rc%udef
         LISWRF_export(n)%sh2o4_t = LIS_rc%udef
         LISWRF_export(n)%relsmc1_t = LIS_rc%udef
         LISWRF_export(n)%relsmc2_t = LIS_rc%udef
         LISWRF_export(n)%relsmc3_t = LIS_rc%udef
         LISWRF_export(n)%relsmc4_t = LIS_rc%udef
         LISWRF_export(n)%xland_t = LIS_rc%udef
#ifdef WRF_HYDRO
         LISWRF_export(n)%infxsrt_t = LIS_rc%udef
         LISWRF_export(n)%soldrain_t = LIS_rc%udef
#endif
#ifdef PARFLOW
         LISWRF_export(n)%wtrflx1_t = LIS_rc%udef
         LISWRF_export(n)%wtrflx2_t = LIS_rc%udef
         LISWRF_export(n)%wtrflx3_t = LIS_rc%udef
         LISWRF_export(n)%wtrflx4_t = LIS_rc%udef
#endif
      enddo
    end subroutine LISWRF_reset_states

end module LISWRFGridCompMod
