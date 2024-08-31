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
! !ROUTINE: jules50_finalize
! \label{jules50_finalize}
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES 5.0 
!
! !INTERFACE:
subroutine jules50_finalize()
! !USES:
   use LIS_coreMod,    only : LIS_rc
   use jules50_lsmMod, only : jules50_struc

! !ARGUMENTS:

!
! !DESCRIPTION:
!
! This routine cleans up the allocated memory structures in
! the jules50 (forcing-only option)
!EOP
   implicit none

   integer :: n, t

#if 0
   do n = 1, lis_rc%nnest
      do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         deallocate(jules50_struc(n)%jules50(t)%nsnow)
         deallocate(jules50_struc(n)%jules50(t)%tsoil_deep)
         deallocate(jules50_struc(n)%jules50(t)%sice)
         deallocate(jules50_struc(n)%jules50(t)%sliq)
         deallocate(jules50_struc(n)%jules50(t)%snowdepth)
         deallocate(jules50_struc(n)%jules50(t)%tsnow)
         deallocate(jules50_struc(n)%jules50(t)%rgrainl)
         deallocate(jules50_struc(n)%jules50(t)%rho_snow_grnd)
         deallocate(jules50_struc(n)%jules50(t)%rho_snow)
         deallocate(jules50_struc(n)%jules50(t)%ds)
         deallocate(jules50_struc(n)%jules50(t)%ns)
         deallocate(jules50_struc(n)%jules50(t)%canht_ft)
         deallocate(jules50_struc(n)%jules50(t)%canopy)
         deallocate(jules50_struc(n)%jules50(t)%cs)
         deallocate(jules50_struc(n)%jules50(t)%di_ncat)
         deallocate(jules50_struc(n)%jules50(t)%k_sice)
         deallocate(jules50_struc(n)%jules50(t)%gc)
         deallocate(jules50_struc(n)%jules50(t)%lai)
         deallocate(jules50_struc(n)%jules50(t)%rgrain)
         deallocate(jules50_struc(n)%jules50(t)%smcl)
         deallocate(jules50_struc(n)%jules50(t)%snow_tile)
         deallocate(jules50_struc(n)%jules50(t)%snow_grnd)
         deallocate(jules50_struc(n)%jules50(t)%snow_mass_sea_ncat)
         deallocate(jules50_struc(n)%jules50(t)%t_soil)
         deallocate(jules50_struc(n)%jules50(t)%tstar_tile)
         deallocate(jules50_struc(n)%jules50(t)%p_s_b)
         deallocate(jules50_struc(n)%jules50(t)%p_s_sathh)
         deallocate(jules50_struc(n)%jules50(t)%p_s_hcap)
         deallocate(jules50_struc(n)%jules50(t)%p_s_hcon)
         deallocate(jules50_struc(n)%jules50(t)%p_s_satcon)
         deallocate(jules50_struc(n)%jules50(t)%p_s_smvccl)
         deallocate(jules50_struc(n)%jules50(t)%p_s_smvcst)
         deallocate(jules50_struc(n)%jules50(t)%p_s_smvcwt)
         deallocate(jules50_struc(n)%jules50(t)%p_s_catch)
         deallocate(jules50_struc(n)%jules50(t)%p_s_catch_snow)
         deallocate(jules50_struc(n)%jules50(t)%p_s_infil_tile)
         deallocate(jules50_struc(n)%jules50(t)%p_s_z0_tile)
         deallocate(jules50_struc(n)%jules50(t)%p_s_z0h_tile_bare)
         deallocate(jules50_struc(n)%jules50(t)%p_s_sthu)
         deallocate(jules50_struc(n)%jules50(t)%p_s_sthf)
         deallocate(jules50_struc(n)%jules50(t)%p_s_sthu_min)
         deallocate(jules50_struc(n)%jules50(t)%sice_pts_ncat)
         deallocate(jules50_struc(n)%jules50(t)%sice_index_ncat)
         deallocate(jules50_struc(n)%jules50(t)%sice_frac_ncat)
         deallocate(jules50_struc(n)%jules50(t)%tile_index)
         deallocate(jules50_struc(n)%jules50(t)%tile_pts)
         deallocate(jules50_struc(n)%jules50(t)%frac)
         deallocate(jules50_struc(n)%jules50(t)%ice_fract_ncat)
         deallocate(jules50_struc(n)%jules50(t)%ti_cat)
         deallocate(jules50_struc(n)%jules50(t)%pond_frac_cat)
         deallocate(jules50_struc(n)%jules50(t)%pond_depth_cat)
         deallocate(jules50_struc(n)%jules50(t)%isoprene_ft)
         deallocate(jules50_struc(n)%jules50(t)%terpene_ft)
         deallocate(jules50_struc(n)%jules50(t)%methanol_ft)
         deallocate(jules50_struc(n)%jules50(t)%acetone_ft)
         deallocate(jules50_struc(n)%jules50(t)%g_leaf_acc)
         deallocate(jules50_struc(n)%jules50(t)%npp_ft_acc)
         deallocate(jules50_struc(n)%jules50(t)%g_leaf_phen_acc)
         deallocate(jules50_struc(n)%jules50(t)%resp_w_ft_acc)
         deallocate(jules50_struc(n)%jules50(t)%resp_s_acc)
         deallocate(jules50_struc(n)%jules50(t)%g_leaf)
         deallocate(jules50_struc(n)%jules50(t)%g_leaf_phen)
         deallocate(jules50_struc(n)%jules50(t)%gpp_ft)
         deallocate(jules50_struc(n)%jules50(t)%npp_ft)
         deallocate(jules50_struc(n)%jules50(t)%resp_p_ft)
         deallocate(jules50_struc(n)%jules50(t)%resp_s)
         deallocate(jules50_struc(n)%jules50(t)%resp_w_ft)
         deallocate(jules50_struc(n)%jules50(t)%lai_phen)
         deallocate(jules50_struc(n)%jules50(t)%c_veg)
         deallocate(jules50_struc(n)%jules50(t)%g_leaf_day)
         deallocate(jules50_struc(n)%jules50(t)%g_leaf_dr_out)
         deallocate(jules50_struc(n)%jules50(t)%lit_c)
         deallocate(jules50_struc(n)%jules50(t)%npp_dr_out)
         deallocate(jules50_struc(n)%jules50(t)%resp_w_dr_out)
         deallocate(jules50_struc(n)%jules50(t)%resp_s_dr_out)
         deallocate(jules50_struc(n)%jules50(t)%unload_backgrnd)
         deallocate(jules50_struc(n)%jules50(t)%alb_tile)
         deallocate(jules50_struc(n)%jules50(t)%fsmc)
         deallocate(jules50_struc(n)%jules50(t)%ftl_tile)
         deallocate(jules50_struc(n)%jules50(t)%le_tile)
         deallocate(jules50_struc(n)%jules50(t)%fqw_tile)
         deallocate(jules50_struc(n)%jules50(t)%fqw_ice)
         deallocate(jules50_struc(n)%jules50(t)%ftl_ice)
         deallocate(jules50_struc(n)%jules50(t)%esoil_tile)
         deallocate(jules50_struc(n)%jules50(t)%sea_ice_htf)
         deallocate(jules50_struc(n)%jules50(t)%surf_htf_tile)
         deallocate(jules50_struc(n)%jules50(t)%land_albedo)
         deallocate(jules50_struc(n)%jules50(t)%ei_tile)
         deallocate(jules50_struc(n)%jules50(t)%ecan_tile)
         deallocate(jules50_struc(n)%jules50(t)%ext)
         deallocate(jules50_struc(n)%jules50(t)%radnet_tile)
         deallocate(jules50_struc(n)%jules50(t)%sw_tile)
         deallocate(jules50_struc(n)%jules50(t)%emis_tile)
         deallocate(jules50_struc(n)%jules50(t)%melt_tile)
      enddo
      deallocate(jules50_struc(n)%jules50)
   enddo
   deallocate(jules50_struc)
#endif

end subroutine jules50_finalize
