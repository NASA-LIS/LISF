!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

subroutine tile_to_trifctl(n,t,pft)
  use trifctl 
  use jules50_lsmMod
  use ancil_info, only : npft_trif 
  use jules_surface_types_mod,  only: npft, nnvg, ntype
  implicit none 
  integer :: n, t, pft
 
  asteps_since_triffid        = jules50_struc(n)%jules50(t)%asteps_since_triffid     ! Number of atmospheric timesteps since last call to TRIFFID
  resp_s_acc_soilt(1,1,:, :)  = jules50_struc(n)%jules50(t)%resp_s_acc_soilt(:,:)            ! Accumulated RESP_S
  gpp_gb(1)                   = jules50_struc(n)%jules50(t)%gpp                      ! Gross primary productivity (kg C/m2/s)
  npp_gb(1)                   = jules50_struc(n)%jules50(t)%npp                      ! Net primary productivity (kg C/m2/s)
  resp_p_gb(1)                = jules50_struc(n)%jules50(t)%resp_p                   ! Plant respiration (kg C/m2/s)
  resp_s_soilt(1,1,:, :)      = jules50_struc(n)%jules50(t)%resp_s_soilt(:,:)                ! Soil respiration (kg C/m2/s)
  cv_gb(1)                    = jules50_struc(n)%jules50(t)%cv                       ! Gridbox mean vegetation carbon (kg C/m2)
  lit_c_mn_gb(1)              = jules50_struc(n)%jules50(t)%lit_c_mn                 ! Gridbox mean carbon litter (kg C/m2/360days)
  resp_s_dr_out_gb(1,:,:)     = jules50_struc(n)%jules50(t)%resp_s_dr_out_gb(:,:)         ! Mean soil respiration for driving TRIFFID (kg C/m2/360days)
  frac_agr_gb(1)              = jules50_struc(n)%jules50(t)%frac_agr                 ! Fraction of agriculture
  
  if(pft .le. npft) then
    g_leaf_acc_pft(1, pft)      = jules50_struc(n)%jules50(t)%g_leaf_acc(pft)           ! Accumulated leaf turnover rate
    g_leaf_phen_acc_pft(1, pft) = jules50_struc(n)%jules50(t)%g_leaf_phen_acc(pft)      ! Accumulated leaf turnover rate including phenology
    g_leaf_pft(1, pft)          = jules50_struc(n)%jules50(t)%g_leaf(pft)               ! Leaf turnover rate (/360days)
    g_leaf_phen_pft(1, pft)     = jules50_struc(n)%jules50(t)%g_leaf_phen(pft)          ! Mean leaf turnover rate over phenology period(/360days)
    resp_p_pft(1, pft)       = jules50_struc(n)%jules50(t)%resp_p_ft(pft)            ! Plant respiration on PFTs (kg C/m2/s)
    gpp_pft(1, pft)          = jules50_struc(n)%jules50(t)%gpp_ft(pft)               ! Gross primary productivity on PFTs (kg C/m2/s)
    npp_pft(1, pft)          = jules50_struc(n)%jules50(t)%npp_ft(pft)               ! Net primary productivity on PFTs (kg C/m2/s)
    resp_w_pft(1, pft)       = jules50_struc(n)%jules50(t)%resp_w_ft(pft)            ! Wood maintenance respiration (kg C/m2/s)
    lai_phen_pft(1, pft)        = jules50_struc(n)%jules50(t)%lai_phen(pft)             ! LAI of PFTs after phenology. Required as separate variable for top-level argument list matching with VEG_IC2A
    c_veg_pft(1, pft)           = jules50_struc(n)%jules50(t)%c_veg(pft)                ! Total carbon content of the vegetation (kg C/m2)
    g_leaf_day_pft(1, pft)      = jules50_struc(n)%jules50(t)%g_leaf_day(pft)           ! Mean leaf turnover rate for input to PHENOL (/360days)
    g_leaf_dr_out_pft(1, pft)   = jules50_struc(n)%jules50(t)%g_leaf_dr_out(pft)        ! Mean leaf turnover rate for driving TRIFFID (/360days)
    lit_c_pft(1, pft)           = jules50_struc(n)%jules50(t)%lit_c(pft)                ! Carbon Litter (kg C/m2/360days)
    npp_dr_out_pft(1, pft)      = jules50_struc(n)%jules50(t)%npp_dr_out(pft)           ! Mean NPP for driving TRIFFID (kg C/m2/360days)
    resp_w_dr_out_pft(1, pft)   = jules50_struc(n)%jules50(t)%resp_w_dr_out(pft)        ! Mean wood respiration for driving TRIFFID (kg C/m2/360days)
  endif

  if(npft_trif>1 .and. pft .le. npft_trif) then
     npp_acc_pft(1, pft)     = jules50_struc(n)%jules50(t)%npp_ft_acc(pft)            ! Accumulated NPP_FT
     resp_w_acc_pft(1, pft)  = jules50_struc(n)%jules50(t)%resp_w_ft_acc(pft)         ! Accum RESP_W_FT
  else                      
     npp_acc_pft(1, :)       = jules50_struc(n)%jules50(t)%npp_ft_acc(:)         ! Accumulated NPP_FT
     resp_w_acc_pft(1, :)    = jules50_struc(n)%jules50(t)%resp_w_ft_acc(:)      ! Accum RESP_W_FT
  endif                     

end subroutine tile_to_trifctl
