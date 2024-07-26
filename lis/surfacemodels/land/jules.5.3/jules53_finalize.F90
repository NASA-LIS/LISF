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
! !ROUTINE: jules53_finalize
! \label{jules53_finalize}
!
! !INTERFACE:
subroutine jules53_finalize()
! !USES:
   use LIS_coreMod,    only : LIS_rc
   use jules53_lsmMod, only : jules53_struc

! !ARGUMENTS:

!
! !DESCRIPTION:
!
! This routine cleans up the allocated memory structures in
! the jules53 (forcing-only option)
!EOP
   implicit none

   integer :: n, t

#if 0
   do n = 1, lis_rc%nnest
      do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         deallocate(jules53_struc(n)%jules53(t)%nsnow)
         deallocate(jules53_struc(n)%jules53(t)%tsoil_deep)
         deallocate(jules53_struc(n)%jules53(t)%sice)
         deallocate(jules53_struc(n)%jules53(t)%sliq)
         deallocate(jules53_struc(n)%jules53(t)%snowdepth)
         deallocate(jules53_struc(n)%jules53(t)%tsnow)
         deallocate(jules53_struc(n)%jules53(t)%rgrainl)
         deallocate(jules53_struc(n)%jules53(t)%rho_snow_grnd)
         deallocate(jules53_struc(n)%jules53(t)%rho_snow)
         deallocate(jules53_struc(n)%jules53(t)%ds)
         deallocate(jules53_struc(n)%jules53(t)%ns)
         deallocate(jules53_struc(n)%jules53(t)%canht_ft)
         deallocate(jules53_struc(n)%jules53(t)%canopy)
         deallocate(jules53_struc(n)%jules53(t)%cs)
         deallocate(jules53_struc(n)%jules53(t)%di_ncat)
         deallocate(jules53_struc(n)%jules53(t)%k_sice)
         deallocate(jules53_struc(n)%jules53(t)%gc)
         deallocate(jules53_struc(n)%jules53(t)%lai)
         deallocate(jules53_struc(n)%jules53(t)%rgrain)
         deallocate(jules53_struc(n)%jules53(t)%smcl)
         deallocate(jules53_struc(n)%jules53(t)%snow_tile)
         deallocate(jules53_struc(n)%jules53(t)%snow_grnd)
         deallocate(jules53_struc(n)%jules53(t)%snow_mass_sea_ncat)
         deallocate(jules53_struc(n)%jules53(t)%t_soil)
         deallocate(jules53_struc(n)%jules53(t)%tstar_tile)
         deallocate(jules53_struc(n)%jules53(t)%p_s_b)
         deallocate(jules53_struc(n)%jules53(t)%p_s_sathh)
         deallocate(jules53_struc(n)%jules53(t)%p_s_hcap)
         deallocate(jules53_struc(n)%jules53(t)%p_s_hcon)
         deallocate(jules53_struc(n)%jules53(t)%p_s_satcon)
         deallocate(jules53_struc(n)%jules53(t)%p_s_smvccl)
         deallocate(jules53_struc(n)%jules53(t)%p_s_smvcst)
         deallocate(jules53_struc(n)%jules53(t)%p_s_smvcwt)
         deallocate(jules53_struc(n)%jules53(t)%p_s_catch)
         deallocate(jules53_struc(n)%jules53(t)%p_s_catch_snow)
         deallocate(jules53_struc(n)%jules53(t)%p_s_infil_tile)
         deallocate(jules53_struc(n)%jules53(t)%p_s_z0_tile)
         deallocate(jules53_struc(n)%jules53(t)%p_s_z0h_tile_bare)
         deallocate(jules53_struc(n)%jules53(t)%p_s_sthu)
         deallocate(jules53_struc(n)%jules53(t)%p_s_sthf)
         deallocate(jules53_struc(n)%jules53(t)%p_s_sthu_min)
         deallocate(jules53_struc(n)%jules53(t)%sice_pts_ncat)
         deallocate(jules53_struc(n)%jules53(t)%sice_index_ncat)
         deallocate(jules53_struc(n)%jules53(t)%sice_frac_ncat)
         deallocate(jules53_struc(n)%jules53(t)%tile_index)
         deallocate(jules53_struc(n)%jules53(t)%tile_pts)
         deallocate(jules53_struc(n)%jules53(t)%frac)
         deallocate(jules53_struc(n)%jules53(t)%ice_fract_ncat)
         deallocate(jules53_struc(n)%jules53(t)%ti_cat)
         deallocate(jules53_struc(n)%jules53(t)%pond_frac_cat)
         deallocate(jules53_struc(n)%jules53(t)%pond_depth_cat)
         deallocate(jules53_struc(n)%jules53(t)%isoprene_ft)
         deallocate(jules53_struc(n)%jules53(t)%terpene_ft)
         deallocate(jules53_struc(n)%jules53(t)%methanol_ft)
         deallocate(jules53_struc(n)%jules53(t)%acetone_ft)
         deallocate(jules53_struc(n)%jules53(t)%g_leaf_acc)
         deallocate(jules53_struc(n)%jules53(t)%npp_ft_acc)
         deallocate(jules53_struc(n)%jules53(t)%g_leaf_phen_acc)
         deallocate(jules53_struc(n)%jules53(t)%resp_w_ft_acc)
         deallocate(jules53_struc(n)%jules53(t)%resp_s_acc)
         deallocate(jules53_struc(n)%jules53(t)%g_leaf)
         deallocate(jules53_struc(n)%jules53(t)%g_leaf_phen)
         deallocate(jules53_struc(n)%jules53(t)%gpp_ft)
         deallocate(jules53_struc(n)%jules53(t)%npp_ft)
         deallocate(jules53_struc(n)%jules53(t)%resp_p_ft)
         deallocate(jules53_struc(n)%jules53(t)%resp_s)
         deallocate(jules53_struc(n)%jules53(t)%resp_w_ft)
         deallocate(jules53_struc(n)%jules53(t)%lai_phen)
         deallocate(jules53_struc(n)%jules53(t)%c_veg)
         deallocate(jules53_struc(n)%jules53(t)%g_leaf_day)
         deallocate(jules53_struc(n)%jules53(t)%g_leaf_dr_out)
         deallocate(jules53_struc(n)%jules53(t)%lit_c)
         deallocate(jules53_struc(n)%jules53(t)%npp_dr_out)
         deallocate(jules53_struc(n)%jules53(t)%resp_w_dr_out)
         deallocate(jules53_struc(n)%jules53(t)%resp_s_dr_out)
         deallocate(jules53_struc(n)%jules53(t)%unload_backgrnd)
         deallocate(jules53_struc(n)%jules53(t)%alb_tile)
         deallocate(jules53_struc(n)%jules53(t)%fsmc)
         deallocate(jules53_struc(n)%jules53(t)%ftl_tile)
         deallocate(jules53_struc(n)%jules53(t)%le_tile)
         deallocate(jules53_struc(n)%jules53(t)%fqw_tile)
         deallocate(jules53_struc(n)%jules53(t)%fqw_ice)
         deallocate(jules53_struc(n)%jules53(t)%ftl_ice)
         deallocate(jules53_struc(n)%jules53(t)%esoil_tile)
         deallocate(jules53_struc(n)%jules53(t)%sea_ice_htf)
         deallocate(jules53_struc(n)%jules53(t)%surf_htf_tile)
         deallocate(jules53_struc(n)%jules53(t)%land_albedo)
         deallocate(jules53_struc(n)%jules53(t)%ei_tile)
         deallocate(jules53_struc(n)%jules53(t)%ecan_tile)
         deallocate(jules53_struc(n)%jules53(t)%ext)
         deallocate(jules53_struc(n)%jules53(t)%radnet_tile)
         deallocate(jules53_struc(n)%jules53(t)%sw_tile)
         deallocate(jules53_struc(n)%jules53(t)%emis_tile)
         deallocate(jules53_struc(n)%jules53(t)%melt_tile)
      enddo
      deallocate(jules53_struc(n)%jules53)
   enddo
   deallocate(jules53_struc)
#endif

end subroutine jules53_finalize
