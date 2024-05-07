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
! !ROUTINE: jules51_finalize
! \label{jules51_finalize}
! !REVISION HISTORY:
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES.5.1 
!
! !INTERFACE:
subroutine jules51_finalize()
! !USES:
   use LIS_coreMod,    only : LIS_rc
   use jules51_lsmMod, only : jules51_struc

! !ARGUMENTS:

!
! !DESCRIPTION:
!
! This routine cleans up the allocated memory structures in
! the jules51 (forcing-only option)
!EOP
   implicit none

   integer :: n, t

#if 0
   do n = 1, lis_rc%nnest
      do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         deallocate(jules51_struc(n)%jules51(t)%nsnow)
         deallocate(jules51_struc(n)%jules51(t)%tsoil_deep)
         deallocate(jules51_struc(n)%jules51(t)%sice)
         deallocate(jules51_struc(n)%jules51(t)%sliq)
         deallocate(jules51_struc(n)%jules51(t)%snowdepth)
         deallocate(jules51_struc(n)%jules51(t)%tsnow)
         deallocate(jules51_struc(n)%jules51(t)%rgrainl)
         deallocate(jules51_struc(n)%jules51(t)%rho_snow_grnd)
         deallocate(jules51_struc(n)%jules51(t)%rho_snow)
         deallocate(jules51_struc(n)%jules51(t)%ds)
         deallocate(jules51_struc(n)%jules51(t)%ns)
         deallocate(jules51_struc(n)%jules51(t)%canht_ft)
         deallocate(jules51_struc(n)%jules51(t)%canopy)
         deallocate(jules51_struc(n)%jules51(t)%cs)
         deallocate(jules51_struc(n)%jules51(t)%di_ncat)
         deallocate(jules51_struc(n)%jules51(t)%k_sice)
         deallocate(jules51_struc(n)%jules51(t)%gc)
         deallocate(jules51_struc(n)%jules51(t)%lai)
         deallocate(jules51_struc(n)%jules51(t)%rgrain)
         deallocate(jules51_struc(n)%jules51(t)%smcl)
         deallocate(jules51_struc(n)%jules51(t)%snow_tile)
         deallocate(jules51_struc(n)%jules51(t)%snow_grnd)
         deallocate(jules51_struc(n)%jules51(t)%snow_mass_sea_ncat)
         deallocate(jules51_struc(n)%jules51(t)%t_soil)
         deallocate(jules51_struc(n)%jules51(t)%tstar_tile)
         deallocate(jules51_struc(n)%jules51(t)%p_s_b)
         deallocate(jules51_struc(n)%jules51(t)%p_s_sathh)
         deallocate(jules51_struc(n)%jules51(t)%p_s_hcap)
         deallocate(jules51_struc(n)%jules51(t)%p_s_hcon)
         deallocate(jules51_struc(n)%jules51(t)%p_s_satcon)
         deallocate(jules51_struc(n)%jules51(t)%p_s_smvccl)
         deallocate(jules51_struc(n)%jules51(t)%p_s_smvcst)
         deallocate(jules51_struc(n)%jules51(t)%p_s_smvcwt)
         deallocate(jules51_struc(n)%jules51(t)%p_s_catch)
         deallocate(jules51_struc(n)%jules51(t)%p_s_catch_snow)
         deallocate(jules51_struc(n)%jules51(t)%p_s_infil_tile)
         deallocate(jules51_struc(n)%jules51(t)%p_s_z0_tile)
         deallocate(jules51_struc(n)%jules51(t)%p_s_z0h_tile_bare)
         deallocate(jules51_struc(n)%jules51(t)%p_s_sthu)
         deallocate(jules51_struc(n)%jules51(t)%p_s_sthf)
         deallocate(jules51_struc(n)%jules51(t)%p_s_sthu_min)
         deallocate(jules51_struc(n)%jules51(t)%sice_pts_ncat)
         deallocate(jules51_struc(n)%jules51(t)%sice_index_ncat)
         deallocate(jules51_struc(n)%jules51(t)%sice_frac_ncat)
         deallocate(jules51_struc(n)%jules51(t)%tile_index)
         deallocate(jules51_struc(n)%jules51(t)%tile_pts)
         deallocate(jules51_struc(n)%jules51(t)%frac)
         deallocate(jules51_struc(n)%jules51(t)%ice_fract_ncat)
         deallocate(jules51_struc(n)%jules51(t)%ti_cat)
         deallocate(jules51_struc(n)%jules51(t)%pond_frac_cat)
         deallocate(jules51_struc(n)%jules51(t)%pond_depth_cat)
         deallocate(jules51_struc(n)%jules51(t)%isoprene_ft)
         deallocate(jules51_struc(n)%jules51(t)%terpene_ft)
         deallocate(jules51_struc(n)%jules51(t)%methanol_ft)
         deallocate(jules51_struc(n)%jules51(t)%acetone_ft)
         deallocate(jules51_struc(n)%jules51(t)%g_leaf_acc)
         deallocate(jules51_struc(n)%jules51(t)%npp_ft_acc)
         deallocate(jules51_struc(n)%jules51(t)%g_leaf_phen_acc)
         deallocate(jules51_struc(n)%jules51(t)%resp_w_ft_acc)
         deallocate(jules51_struc(n)%jules51(t)%resp_s_acc)
         deallocate(jules51_struc(n)%jules51(t)%g_leaf)
         deallocate(jules51_struc(n)%jules51(t)%g_leaf_phen)
         deallocate(jules51_struc(n)%jules51(t)%gpp_ft)
         deallocate(jules51_struc(n)%jules51(t)%npp_ft)
         deallocate(jules51_struc(n)%jules51(t)%resp_p_ft)
         deallocate(jules51_struc(n)%jules51(t)%resp_s)
         deallocate(jules51_struc(n)%jules51(t)%resp_w_ft)
         deallocate(jules51_struc(n)%jules51(t)%lai_phen)
         deallocate(jules51_struc(n)%jules51(t)%c_veg)
         deallocate(jules51_struc(n)%jules51(t)%g_leaf_day)
         deallocate(jules51_struc(n)%jules51(t)%g_leaf_dr_out)
         deallocate(jules51_struc(n)%jules51(t)%lit_c)
         deallocate(jules51_struc(n)%jules51(t)%npp_dr_out)
         deallocate(jules51_struc(n)%jules51(t)%resp_w_dr_out)
         deallocate(jules51_struc(n)%jules51(t)%resp_s_dr_out)
         deallocate(jules51_struc(n)%jules51(t)%unload_backgrnd)
         deallocate(jules51_struc(n)%jules51(t)%alb_tile)
         deallocate(jules51_struc(n)%jules51(t)%fsmc)
         deallocate(jules51_struc(n)%jules51(t)%ftl_tile)
         deallocate(jules51_struc(n)%jules51(t)%le_tile)
         deallocate(jules51_struc(n)%jules51(t)%fqw_tile)
         deallocate(jules51_struc(n)%jules51(t)%fqw_ice)
         deallocate(jules51_struc(n)%jules51(t)%ftl_ice)
         deallocate(jules51_struc(n)%jules51(t)%esoil_tile)
         deallocate(jules51_struc(n)%jules51(t)%sea_ice_htf)
         deallocate(jules51_struc(n)%jules51(t)%surf_htf_tile)
         deallocate(jules51_struc(n)%jules51(t)%land_albedo)
         deallocate(jules51_struc(n)%jules51(t)%ei_tile)
         deallocate(jules51_struc(n)%jules51(t)%ecan_tile)
         deallocate(jules51_struc(n)%jules51(t)%ext)
         deallocate(jules51_struc(n)%jules51(t)%radnet_tile)
         deallocate(jules51_struc(n)%jules51(t)%sw_tile)
         deallocate(jules51_struc(n)%jules51(t)%emis_tile)
         deallocate(jules51_struc(n)%jules51(t)%melt_tile)
      enddo
      deallocate(jules51_struc(n)%jules51)
   enddo
   deallocate(jules51_struc)
#endif

end subroutine jules51_finalize
