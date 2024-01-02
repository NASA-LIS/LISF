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
! !ROUTINE: jules54_finalize
! \label{jules54_finalize}
!
! !INTERFACE:
subroutine jules54_finalize()
! !USES:
   use LIS_coreMod,    only : LIS_rc
   use jules54_lsmMod, only : jules54_struc

! !ARGUMENTS:

!
! !DESCRIPTION:
!
! This routine cleans up the allocated memory structures in
! the jules54 (forcing-only option)
!EOP
   implicit none

   integer :: n, t

#if 0
   do n = 1, lis_rc%nnest
      do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         deallocate(jules54_struc(n)%jules54(t)%nsnow)
         deallocate(jules54_struc(n)%jules54(t)%tsoil_deep)
         deallocate(jules54_struc(n)%jules54(t)%sice)
         deallocate(jules54_struc(n)%jules54(t)%sliq)
         deallocate(jules54_struc(n)%jules54(t)%snowdepth)
         deallocate(jules54_struc(n)%jules54(t)%tsnow)
         deallocate(jules54_struc(n)%jules54(t)%rgrainl)
         deallocate(jules54_struc(n)%jules54(t)%rho_snow_grnd)
         deallocate(jules54_struc(n)%jules54(t)%rho_snow)
         deallocate(jules54_struc(n)%jules54(t)%ds)
         deallocate(jules54_struc(n)%jules54(t)%ns)
         deallocate(jules54_struc(n)%jules54(t)%canht_ft)
         deallocate(jules54_struc(n)%jules54(t)%canopy)
         deallocate(jules54_struc(n)%jules54(t)%cs)
         deallocate(jules54_struc(n)%jules54(t)%di_ncat)
         deallocate(jules54_struc(n)%jules54(t)%k_sice)
         deallocate(jules54_struc(n)%jules54(t)%gc)
         deallocate(jules54_struc(n)%jules54(t)%lai)
         deallocate(jules54_struc(n)%jules54(t)%rgrain)
         deallocate(jules54_struc(n)%jules54(t)%smcl)
         deallocate(jules54_struc(n)%jules54(t)%snow_tile)
         deallocate(jules54_struc(n)%jules54(t)%snow_grnd)
         deallocate(jules54_struc(n)%jules54(t)%snow_mass_sea_ncat)
         deallocate(jules54_struc(n)%jules54(t)%t_soil)
         deallocate(jules54_struc(n)%jules54(t)%tstar_tile)
         deallocate(jules54_struc(n)%jules54(t)%p_s_b)
         deallocate(jules54_struc(n)%jules54(t)%p_s_sathh)
         deallocate(jules54_struc(n)%jules54(t)%p_s_hcap)
         deallocate(jules54_struc(n)%jules54(t)%p_s_hcon)
         deallocate(jules54_struc(n)%jules54(t)%p_s_satcon)
         deallocate(jules54_struc(n)%jules54(t)%p_s_smvccl)
         deallocate(jules54_struc(n)%jules54(t)%p_s_smvcst)
         deallocate(jules54_struc(n)%jules54(t)%p_s_smvcwt)
         deallocate(jules54_struc(n)%jules54(t)%p_s_catch)
         deallocate(jules54_struc(n)%jules54(t)%p_s_catch_snow)
         deallocate(jules54_struc(n)%jules54(t)%p_s_infil_tile)
         deallocate(jules54_struc(n)%jules54(t)%p_s_z0_tile)
         deallocate(jules54_struc(n)%jules54(t)%p_s_z0h_tile_bare)
         deallocate(jules54_struc(n)%jules54(t)%p_s_sthu)
         deallocate(jules54_struc(n)%jules54(t)%p_s_sthf)
         deallocate(jules54_struc(n)%jules54(t)%p_s_sthu_min)
         deallocate(jules54_struc(n)%jules54(t)%sice_pts_ncat)
         deallocate(jules54_struc(n)%jules54(t)%sice_index_ncat)
         deallocate(jules54_struc(n)%jules54(t)%sice_frac_ncat)
         deallocate(jules54_struc(n)%jules54(t)%tile_index)
         deallocate(jules54_struc(n)%jules54(t)%tile_pts)
         deallocate(jules54_struc(n)%jules54(t)%frac)
         deallocate(jules54_struc(n)%jules54(t)%ice_fract_ncat)
         deallocate(jules54_struc(n)%jules54(t)%ti_cat)
         deallocate(jules54_struc(n)%jules54(t)%pond_frac_cat)
         deallocate(jules54_struc(n)%jules54(t)%pond_depth_cat)
         deallocate(jules54_struc(n)%jules54(t)%isoprene_ft)
         deallocate(jules54_struc(n)%jules54(t)%terpene_ft)
         deallocate(jules54_struc(n)%jules54(t)%methanol_ft)
         deallocate(jules54_struc(n)%jules54(t)%acetone_ft)
         deallocate(jules54_struc(n)%jules54(t)%g_leaf_acc)
         deallocate(jules54_struc(n)%jules54(t)%npp_ft_acc)
         deallocate(jules54_struc(n)%jules54(t)%g_leaf_phen_acc)
         deallocate(jules54_struc(n)%jules54(t)%resp_w_ft_acc)
         deallocate(jules54_struc(n)%jules54(t)%resp_s_acc)
         deallocate(jules54_struc(n)%jules54(t)%g_leaf)
         deallocate(jules54_struc(n)%jules54(t)%g_leaf_phen)
         deallocate(jules54_struc(n)%jules54(t)%gpp_ft)
         deallocate(jules54_struc(n)%jules54(t)%npp_ft)
         deallocate(jules54_struc(n)%jules54(t)%resp_p_ft)
         deallocate(jules54_struc(n)%jules54(t)%resp_s)
         deallocate(jules54_struc(n)%jules54(t)%resp_w_ft)
         deallocate(jules54_struc(n)%jules54(t)%lai_phen)
         deallocate(jules54_struc(n)%jules54(t)%c_veg)
         deallocate(jules54_struc(n)%jules54(t)%g_leaf_day)
         deallocate(jules54_struc(n)%jules54(t)%g_leaf_dr_out)
         deallocate(jules54_struc(n)%jules54(t)%lit_c)
         deallocate(jules54_struc(n)%jules54(t)%npp_dr_out)
         deallocate(jules54_struc(n)%jules54(t)%resp_w_dr_out)
         deallocate(jules54_struc(n)%jules54(t)%resp_s_dr_out)
         deallocate(jules54_struc(n)%jules54(t)%unload_backgrnd)
         deallocate(jules54_struc(n)%jules54(t)%alb_tile)
         deallocate(jules54_struc(n)%jules54(t)%fsmc)
         deallocate(jules54_struc(n)%jules54(t)%ftl_tile)
         deallocate(jules54_struc(n)%jules54(t)%le_tile)
         deallocate(jules54_struc(n)%jules54(t)%fqw_tile)
         deallocate(jules54_struc(n)%jules54(t)%fqw_ice)
         deallocate(jules54_struc(n)%jules54(t)%ftl_ice)
         deallocate(jules54_struc(n)%jules54(t)%esoil_tile)
         deallocate(jules54_struc(n)%jules54(t)%sea_ice_htf)
         deallocate(jules54_struc(n)%jules54(t)%surf_htf_tile)
         deallocate(jules54_struc(n)%jules54(t)%land_albedo)
         deallocate(jules54_struc(n)%jules54(t)%ei_tile)
         deallocate(jules54_struc(n)%jules54(t)%ecan_tile)
         deallocate(jules54_struc(n)%jules54(t)%ext)
         deallocate(jules54_struc(n)%jules54(t)%radnet_tile)
         deallocate(jules54_struc(n)%jules54(t)%sw_tile)
         deallocate(jules54_struc(n)%jules54(t)%emis_tile)
         deallocate(jules54_struc(n)%jules54(t)%melt_tile)
      enddo
      deallocate(jules54_struc(n)%jules54)
   enddo
   deallocate(jules54_struc)
#endif

end subroutine jules54_finalize
