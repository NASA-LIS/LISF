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
! !ROUTINE: clm2_setwrfexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from CLM to WRF in coupled mode
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 17 Nov 2008; Sujay Kumar, Modified for the ESMF coupled version
! 
! !INTERFACE:
subroutine clm2_setwrfexport(n, expState)
! !USES:
  use ESMF
  use LIS_precisionMod
  use LIS_coreMod,    only : LIS_rc, LIS_surface
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_logMod,     only : LIS_verify
  use LISWRFGridCompMod, only : LISWRF_export
  use clm2_lsmMod
  use clm2_varcon, only : denh2o, denice, hvap, hsub, hfus, istwet 
  use clm2_varpar, only : nlevsoi
!EOP

  implicit none
  integer, intent(in) :: n 
  type(ESMF_State)      :: expState
  integer               :: i,j,k,t,m
  real                  :: temp(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real :: asurft(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real :: snowt(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real :: soilmr(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real :: soilmtc(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real :: soilmoist(LIS_rc%npatch(n,LIS_rc%lsm_index),1:nlevsoi)
  real :: q1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real :: RCH

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snowt(t)=0.
     if (clm2_struc(n)%clm(t)%itypwat/=istwet)then 
        if(clm2_struc(n)%clm(t)%snl < 0)then
           snowt(t)=clm2_struc(n)%clm(t)%t_soisno(clm2_struc(n)%clm(t)%snl+1)
        endif
     endif
     if(snowt(t)==0.)snowt(t)=LIS_rc%udef  
  enddo
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(snowt(t).ne.LIS_rc%udef)then
        asurft(t)=clm2_struc(n)%clm(t)%frac_sno*snowt(t)+ & 
             clm2_struc(n)%clm(t)%frac_veg_nosno*clm2_struc(n)%clm(t)%t_veg+  & 
             (1-(clm2_struc(n)%clm(t)%frac_sno+clm2_struc(n)%clm(t)%frac_veg_nosno))* & 
             clm2_struc(n)%clm(t)%t_grnd
     else
        asurft(t)=clm2_struc(n)%clm(t)%frac_veg_nosno*clm2_struc(n)%clm(t)%t_veg+ & 
             (1-clm2_struc(n)%clm(t)%frac_veg_nosno)*clm2_struc(n)%clm(t)%t_grnd
     endif
  enddo

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%avgsurft_t,asurft)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qh_t,clm2_struc(n)%clm%eflx_sh_tot)

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%eta_kinematic_t,clm2_struc(n)%clm%eflx_lh_tot/2.501E6)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qle_t,clm2_struc(n)%clm%eflx_lh_tot)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qg_t,clm2_struc(n)%clm%eflx_soil_grnd)

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%lispor,clm2_struc(n)%clm%watsat(1))

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%chs2_t,clm2_struc(n)%clm%chs2)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%cqs2_t,clm2_struc(n)%clm%cqs2)
!----------------------------------------------------------------------
! Root Zone Soil Moisture (kg m-2)
! Calculation of root zone soil moisture 
!----------------------------------------------------------------------
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     soilmr(t)=0.
     do m=1,nlevsoi
        soilmr(t)=soilmr(t)+clm2_struc(n)%clm(t)%rootfr(m)*clm2_struc(n)%clm(t)%h2osoi_liq(m)
     enddo
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%rootmoist_t,soilmr) 
  do m=1,nlevsoi 
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        soilmoist(t,m)=clm2_struc(n)%clm(t)%h2osoi_liq(m)+clm2_struc(n)%clm(t)%h2osoi_ice(m)
     enddo
  enddo
  
  do m=1,nlevsoi
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        soilmtc(t)=soilmoist(t,m)
     enddo
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%soilm_t,soilmtc)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qs_t,clm2_struc(n)%clm%qflx_surf+clm2_struc(n)%clm%qflx_qrgwl)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qsb_t,clm2_struc(n)%clm%qflx_drain)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%albedo_t,clm2_struc(n)%clm%surfalb)

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(clm2_struc(n)%clm(i)%itypveg.ne.LIS_rc%waterclass) then 
        temp(i) = clm2_struc(n)%clm(i)%z0mr*clm2_struc(n)%clm(i)%htop
     endif
     if(temp(i).eq.0) then 
        temp(i) = 0.001
     endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%znt_t,temp)

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%snocvr_t,real(1-clm2_struc(n)%clm%frac_veg_nosno))
!mixing  ratio
  do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(clm2_struc(n)%clm(k)%forc_q.gt.1E-19) then  
        RCH = clm2_struc(n)%clm(k)%forc_pbot/(287.04*clm2_struc(n)%clm(k)%forc_t*(1+0.61*clm2_struc(n)%clm(k)%forc_q)) 
        q1(k) = clm2_struc(n)%clm(k)%forc_q+(clm2_struc(n)%clm(k)%eflx_lh_tot/2.501E6)/(RCH*clm2_struc(n)%clm(k)%forc_ch)
     else
        q1(k) = 0.0
     endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%q1_t,q1)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%snow_t,clm2_struc(n)%clm%h2osno)
!canopy interception
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%cmc_t,clm2_struc(n)%clm%h2ocan)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qsm_t,clm2_struc(n)%clm%qflx_snomelt)

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc4_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc4_t,temp)

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o4_t,temp)

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = LIS_surface(n,LIS_rc%lsm_index)%tile(i)%vegt
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%xland_t,temp)

end subroutine clm2_setwrfexport
 
