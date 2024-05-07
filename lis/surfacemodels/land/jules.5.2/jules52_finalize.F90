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
! !ROUTINE: jules52_finalize
! \label{jules52_finalize}
!
! !INTERFACE:
subroutine jules52_finalize()
! !USES:
   use LIS_coreMod,    only : LIS_rc
   use jules52_lsmMod, only : jules52_struc

! !ARGUMENTS:

!
! !DESCRIPTION:
!
! This routine cleans up the allocated memory structures in
! the jules52 (forcing-only option)
!EOP
   implicit none

   integer :: n, t

#if 0
   do n = 1, lis_rc%nnest
      do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         deallocate(jules52_struc(n)%jules52(t)%nsnow)
         deallocate(jules52_struc(n)%jules52(t)%tsoil_deep)
         deallocate(jules52_struc(n)%jules52(t)%sice)
         deallocate(jules52_struc(n)%jules52(t)%sliq)
         deallocate(jules52_struc(n)%jules52(t)%snowdepth)
         deallocate(jules52_struc(n)%jules52(t)%tsnow)
         deallocate(jules52_struc(n)%jules52(t)%rgrainl)
         deallocate(jules52_struc(n)%jules52(t)%rho_snow_grnd)
         deallocate(jules52_struc(n)%jules52(t)%rho_snow)
         deallocate(jules52_struc(n)%jules52(t)%ds)
         deallocate(jules52_struc(n)%jules52(t)%ns)
         deallocate(jules52_struc(n)%jules52(t)%canht_ft)
         deallocate(jules52_struc(n)%jules52(t)%canopy)
         deallocate(jules52_struc(n)%jules52(t)%cs)
         deallocate(jules52_struc(n)%jules52(t)%di_ncat)
         deallocate(jules52_struc(n)%jules52(t)%k_sice)
         deallocate(jules52_struc(n)%jules52(t)%gc)
         deallocate(jules52_struc(n)%jules52(t)%lai)
         deallocate(jules52_struc(n)%jules52(t)%rgrain)
         deallocate(jules52_struc(n)%jules52(t)%smcl)
         deallocate(jules52_struc(n)%jules52(t)%snow_tile)
         deallocate(jules52_struc(n)%jules52(t)%snow_grnd)
         deallocate(jules52_struc(n)%jules52(t)%snow_mass_sea_ncat)
         deallocate(jules52_struc(n)%jules52(t)%t_soil)
         deallocate(jules52_struc(n)%jules52(t)%tstar_tile)
         deallocate(jules52_struc(n)%jules52(t)%p_s_b)
         deallocate(jules52_struc(n)%jules52(t)%p_s_sathh)
         deallocate(jules52_struc(n)%jules52(t)%p_s_hcap)
         deallocate(jules52_struc(n)%jules52(t)%p_s_hcon)
         deallocate(jules52_struc(n)%jules52(t)%p_s_satcon)
         deallocate(jules52_struc(n)%jules52(t)%p_s_smvccl)
         deallocate(jules52_struc(n)%jules52(t)%p_s_smvcst)
         deallocate(jules52_struc(n)%jules52(t)%p_s_smvcwt)
         deallocate(jules52_struc(n)%jules52(t)%p_s_catch)
         deallocate(jules52_struc(n)%jules52(t)%p_s_catch_snow)
         deallocate(jules52_struc(n)%jules52(t)%p_s_infil_tile)
         deallocate(jules52_struc(n)%jules52(t)%p_s_z0_tile)
         deallocate(jules52_struc(n)%jules52(t)%p_s_z0h_tile_bare)
         deallocate(jules52_struc(n)%jules52(t)%p_s_sthu)
         deallocate(jules52_struc(n)%jules52(t)%p_s_sthf)
         deallocate(jules52_struc(n)%jules52(t)%p_s_sthu_min)
         deallocate(jules52_struc(n)%jules52(t)%sice_pts_ncat)
         deallocate(jules52_struc(n)%jules52(t)%sice_index_ncat)
         deallocate(jules52_struc(n)%jules52(t)%sice_frac_ncat)
         deallocate(jules52_struc(n)%jules52(t)%tile_index)
         deallocate(jules52_struc(n)%jules52(t)%tile_pts)
         deallocate(jules52_struc(n)%jules52(t)%frac)
         deallocate(jules52_struc(n)%jules52(t)%ice_fract_ncat)
         deallocate(jules52_struc(n)%jules52(t)%ti_cat)
         deallocate(jules52_struc(n)%jules52(t)%pond_frac_cat)
         deallocate(jules52_struc(n)%jules52(t)%pond_depth_cat)
         deallocate(jules52_struc(n)%jules52(t)%isoprene_ft)
         deallocate(jules52_struc(n)%jules52(t)%terpene_ft)
         deallocate(jules52_struc(n)%jules52(t)%methanol_ft)
         deallocate(jules52_struc(n)%jules52(t)%acetone_ft)
         deallocate(jules52_struc(n)%jules52(t)%g_leaf_acc)
         deallocate(jules52_struc(n)%jules52(t)%npp_ft_acc)
         deallocate(jules52_struc(n)%jules52(t)%g_leaf_phen_acc)
         deallocate(jules52_struc(n)%jules52(t)%resp_w_ft_acc)
         deallocate(jules52_struc(n)%jules52(t)%resp_s_acc)
         deallocate(jules52_struc(n)%jules52(t)%g_leaf)
         deallocate(jules52_struc(n)%jules52(t)%g_leaf_phen)
         deallocate(jules52_struc(n)%jules52(t)%gpp_ft)
         deallocate(jules52_struc(n)%jules52(t)%npp_ft)
         deallocate(jules52_struc(n)%jules52(t)%resp_p_ft)
         deallocate(jules52_struc(n)%jules52(t)%resp_s)
         deallocate(jules52_struc(n)%jules52(t)%resp_w_ft)
         deallocate(jules52_struc(n)%jules52(t)%lai_phen)
         deallocate(jules52_struc(n)%jules52(t)%c_veg)
         deallocate(jules52_struc(n)%jules52(t)%g_leaf_day)
         deallocate(jules52_struc(n)%jules52(t)%g_leaf_dr_out)
         deallocate(jules52_struc(n)%jules52(t)%lit_c)
         deallocate(jules52_struc(n)%jules52(t)%npp_dr_out)
         deallocate(jules52_struc(n)%jules52(t)%resp_w_dr_out)
         deallocate(jules52_struc(n)%jules52(t)%resp_s_dr_out)
         deallocate(jules52_struc(n)%jules52(t)%unload_backgrnd)
         deallocate(jules52_struc(n)%jules52(t)%alb_tile)
         deallocate(jules52_struc(n)%jules52(t)%fsmc)
         deallocate(jules52_struc(n)%jules52(t)%ftl_tile)
         deallocate(jules52_struc(n)%jules52(t)%le_tile)
         deallocate(jules52_struc(n)%jules52(t)%fqw_tile)
         deallocate(jules52_struc(n)%jules52(t)%fqw_ice)
         deallocate(jules52_struc(n)%jules52(t)%ftl_ice)
         deallocate(jules52_struc(n)%jules52(t)%esoil_tile)
         deallocate(jules52_struc(n)%jules52(t)%sea_ice_htf)
         deallocate(jules52_struc(n)%jules52(t)%surf_htf_tile)
         deallocate(jules52_struc(n)%jules52(t)%land_albedo)
         deallocate(jules52_struc(n)%jules52(t)%ei_tile)
         deallocate(jules52_struc(n)%jules52(t)%ecan_tile)
         deallocate(jules52_struc(n)%jules52(t)%ext)
         deallocate(jules52_struc(n)%jules52(t)%radnet_tile)
         deallocate(jules52_struc(n)%jules52(t)%sw_tile)
         deallocate(jules52_struc(n)%jules52(t)%emis_tile)
         deallocate(jules52_struc(n)%jules52(t)%melt_tile)
      enddo
      deallocate(jules52_struc(n)%jules52)
   enddo
   deallocate(jules52_struc)
#endif

end subroutine jules52_finalize
