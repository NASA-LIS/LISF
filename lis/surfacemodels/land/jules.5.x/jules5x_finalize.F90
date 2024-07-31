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
! !ROUTINE: jules5x_finalize
! \label{jules5x_finalize}
!
! !INTERFACE:
subroutine jules5x_finalize()
! !USES:
   use LIS_coreMod,    only : LIS_rc
   use jules5x_lsmMod, only : jules5x_struc

! !ARGUMENTS:

!
! !DESCRIPTION:
!
! This routine cleans up the allocated memory structures in
! the jules5x (forcing-only option)
!EOP
   implicit none

   integer :: n, t

#if 0
   do n = 1, lis_rc%nnest
      do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         deallocate(jules5x_struc(n)%jules5x(t)%nsnow)
         deallocate(jules5x_struc(n)%jules5x(t)%tsoil_deep)
         deallocate(jules5x_struc(n)%jules5x(t)%sice)
         deallocate(jules5x_struc(n)%jules5x(t)%sliq)
         deallocate(jules5x_struc(n)%jules5x(t)%snowdepth)
         deallocate(jules5x_struc(n)%jules5x(t)%tsnow)
         deallocate(jules5x_struc(n)%jules5x(t)%rgrainl)
         deallocate(jules5x_struc(n)%jules5x(t)%rho_snow_grnd)
         deallocate(jules5x_struc(n)%jules5x(t)%rho_snow)
         deallocate(jules5x_struc(n)%jules5x(t)%ds)
         deallocate(jules5x_struc(n)%jules5x(t)%ns)
         deallocate(jules5x_struc(n)%jules5x(t)%canht_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%canopy)
         deallocate(jules5x_struc(n)%jules5x(t)%cs)
         deallocate(jules5x_struc(n)%jules5x(t)%di_ncat)
         deallocate(jules5x_struc(n)%jules5x(t)%k_sice)
         deallocate(jules5x_struc(n)%jules5x(t)%gc)
         deallocate(jules5x_struc(n)%jules5x(t)%lai)
         deallocate(jules5x_struc(n)%jules5x(t)%rgrain)
         deallocate(jules5x_struc(n)%jules5x(t)%smcl)
         deallocate(jules5x_struc(n)%jules5x(t)%snow_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%snow_grnd)
         deallocate(jules5x_struc(n)%jules5x(t)%snow_mass_sea_ncat)
         deallocate(jules5x_struc(n)%jules5x(t)%t_soil)
         deallocate(jules5x_struc(n)%jules5x(t)%tstar_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_b)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_sathh)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_hcap)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_hcon)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_satcon)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_smvccl)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_smvcst)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_smvcwt)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_catch)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_catch_snow)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_infil_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_z0_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_z0h_tile_bare)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_sthu)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_sthf)
         deallocate(jules5x_struc(n)%jules5x(t)%p_s_sthu_min)
         deallocate(jules5x_struc(n)%jules5x(t)%sice_pts_ncat)
         deallocate(jules5x_struc(n)%jules5x(t)%sice_index_ncat)
         deallocate(jules5x_struc(n)%jules5x(t)%sice_frac_ncat)
         deallocate(jules5x_struc(n)%jules5x(t)%tile_index)
         deallocate(jules5x_struc(n)%jules5x(t)%tile_pts)
         deallocate(jules5x_struc(n)%jules5x(t)%frac)
         deallocate(jules5x_struc(n)%jules5x(t)%ice_fract_ncat)
         deallocate(jules5x_struc(n)%jules5x(t)%ti_cat)
         deallocate(jules5x_struc(n)%jules5x(t)%pond_frac_cat)
         deallocate(jules5x_struc(n)%jules5x(t)%pond_depth_cat)
         deallocate(jules5x_struc(n)%jules5x(t)%isoprene_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%terpene_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%methanol_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%acetone_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%g_leaf_acc)
         deallocate(jules5x_struc(n)%jules5x(t)%npp_ft_acc)
         deallocate(jules5x_struc(n)%jules5x(t)%g_leaf_phen_acc)
         deallocate(jules5x_struc(n)%jules5x(t)%resp_w_ft_acc)
         deallocate(jules5x_struc(n)%jules5x(t)%resp_s_acc)
         deallocate(jules5x_struc(n)%jules5x(t)%g_leaf)
         deallocate(jules5x_struc(n)%jules5x(t)%g_leaf_phen)
         deallocate(jules5x_struc(n)%jules5x(t)%gpp_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%npp_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%resp_p_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%resp_s)
         deallocate(jules5x_struc(n)%jules5x(t)%resp_w_ft)
         deallocate(jules5x_struc(n)%jules5x(t)%lai_phen)
         deallocate(jules5x_struc(n)%jules5x(t)%c_veg)
         deallocate(jules5x_struc(n)%jules5x(t)%g_leaf_day)
         deallocate(jules5x_struc(n)%jules5x(t)%g_leaf_dr_out)
         deallocate(jules5x_struc(n)%jules5x(t)%lit_c)
         deallocate(jules5x_struc(n)%jules5x(t)%npp_dr_out)
         deallocate(jules5x_struc(n)%jules5x(t)%resp_w_dr_out)
         deallocate(jules5x_struc(n)%jules5x(t)%resp_s_dr_out)
         deallocate(jules5x_struc(n)%jules5x(t)%unload_backgrnd)
         deallocate(jules5x_struc(n)%jules5x(t)%alb_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%fsmc)
         deallocate(jules5x_struc(n)%jules5x(t)%ftl_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%le_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%fqw_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%fqw_ice)
         deallocate(jules5x_struc(n)%jules5x(t)%ftl_ice)
         deallocate(jules5x_struc(n)%jules5x(t)%esoil_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%sea_ice_htf)
         deallocate(jules5x_struc(n)%jules5x(t)%surf_htf_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%land_albedo)
         deallocate(jules5x_struc(n)%jules5x(t)%ei_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%ecan_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%ext)
         deallocate(jules5x_struc(n)%jules5x(t)%radnet_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%sw_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%emis_tile)
         deallocate(jules5x_struc(n)%jules5x(t)%melt_tile)
      enddo
      deallocate(jules5x_struc(n)%jules5x)
   enddo
   deallocate(jules5x_struc)
#endif

end subroutine jules5x_finalize
