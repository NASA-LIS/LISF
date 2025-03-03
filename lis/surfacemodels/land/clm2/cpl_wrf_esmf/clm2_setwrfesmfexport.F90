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
! !ROUTINE: clm2_setwrfesmfexport
! \label{clm2_setwrfesmfexport}
!
! !DESCRIPTION:  
!  LIS CLM data writer: Writes clm output in binary format
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine clm2_setwrfesmfexport(n)
! !USES:
#if 0 
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use LIS_precisionMod
  use lisWRFGridCompMod, only : liswrf_export
  use LIS_historyMod, only : LIS_tile2grid
  use clm2_lsmMod, only : clm2_struc
  use clm2_varcon, only : denh2o, denice, hvap, hsub, hfus, istwet 
  use clm2_varpar, only : nlevsoi
#endif
  implicit none
  integer :: n
#if 0 
  integer :: i,j,k,t,m
  real :: temp(LIS_rc%ntiles(n))
  real :: asurft(LIS_rc%ntiles(n))
  real :: snowt(LIS_rc%ntiles(n))
  real :: soilmr(LIS_rc%ntiles(n))
  real :: soilmtc(LIS_rc%ntiles(n))
  real :: soilmoist(LIS_rc%ntiles(n),1:nlevsoi)
  real :: q1(LIS_rc%ntiles(n))
  real :: RCH
#endif
!EOP

#if 0 
  do t=1,LIS_rc%ntiles(n)
     snowt(t)=0.
     if (clm2_struc(n)%clm(t)%itypwat/=istwet)then 
        if(clm2_struc(n)%clm(t)%snl < 0)then
           snowt(t)=clm2_struc(n)%clm(t)%t_soisno(clm2_struc(n)%clm(t)%snl+1)
        endif
     endif
     if(snowt(t)==0.)snowt(t)=LIS_rc%udef  
  enddo
  do t=1,LIS_rc%ntiles(n)
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
  call LIS_tile2grid(n,liswrf_export(n)%avgsurft,asurft)
  call LIS_tile2grid(n,liswrf_export(n)%qh,clm2_struc(n)%clm%eflx_sh_tot)

  call LIS_tile2grid(n,liswrf_export(n)%eta_kinematic,clm2_struc(n)%clm%eflx_lh_tot/2.501E6)
  call LIS_tile2grid(n,liswrf_export(n)%qle,clm2_struc(n)%clm%eflx_lh_tot)
  call LIS_tile2grid(n,liswrf_export(n)%qg,clm2_struc(n)%clm%eflx_soil_grnd)
!----------------------------------------------------------------------
! Root Zone Soil Moisture (kg m-2)
! Calculation of root zone soil moisture 
!----------------------------------------------------------------------
  do t=1,LIS_rc%ntiles(n)
     soilmr(t)=0.
     do m=1,nlevsoi
        soilmr(t)=soilmr(t)+clm2_struc(n)%clm(t)%rootfr(m)*clm2_struc(n)%clm(t)%h2osoi_liq(m)
     enddo
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%rootmoist,soilmr) 
  do m=1,nlevsoi 
     do t=1,LIS_rc%ntiles(n)
        soilmoist(t,m)=clm2_struc(n)%clm(t)%h2osoi_liq(m)+clm2_struc(n)%clm(t)%h2osoi_ice(m)
     enddo
  enddo
  
  do m=1,nlevsoi
     do t=1,LIS_rc%ntiles(n)
        soilmtc(t)=soilmoist(t,m)
     enddo
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%soilm,soilmtc)
  call LIS_tile2grid(n,liswrf_export(n)%qs,clm2_struc(n)%clm%qflx_surf+clm2_struc(n)%clm%qflx_qrgwl)
  call LIS_tile2grid(n,liswrf_export(n)%qsb,clm2_struc(n)%clm%qflx_drain)
  call LIS_tile2grid(n,liswrf_export(n)%albedo,clm2_struc(n)%clm%surfalb)

  do i=1,LIS_rc%ntiles(n)
     if(clm2_struc(n)%clm(i)%itypveg.ne.LIS_rc%waterclass) then 
        temp(i) = clm2_struc(n)%clm(i)%z0mr*clm2_struc(n)%clm(i)%htop
     endif
     if(temp(i).eq.0) then 
        temp(i) = 0.001
     endif
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%znt,temp)

  call LIS_tile2grid(n,liswrf_export(n)%snocvr,real(1-clm2_struc(n)%clm%frac_veg_nosno))
!mixing  ratio
  do k=1,LIS_rc%ntiles(n)
     if(clm2_struc(n)%clm(k)%forc_q.gt.1E-19) then  
        RCH = clm2_struc(n)%clm(k)%forc_pbot/(287.04*clm2_struc(n)%clm(k)%forc_t*(1+0.61*clm2_struc(n)%clm(k)%forc_q)) 
        q1(k) = clm2_struc(n)%clm(k)%forc_q+(clm2_struc(n)%clm(k)%eflx_lh_tot/2.501E6)/(RCH*clm2_struc(n)%clm(k)%forc_ch)
     else
        q1(k) = 0.0
     endif
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%q1,q1)
  call LIS_tile2grid(n,liswrf_export(n)%snow,clm2_struc(n)%clm%h2osno)
!canopy interception
  call LIS_tile2grid(n,liswrf_export(n)%cmc,clm2_struc(n)%clm%h2ocan)
  call LIS_tile2grid(n,liswrf_export(n)%qsm,clm2_struc(n)%clm%qflx_snomelt)

  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,1)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%smc1,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,2)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%smc2,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,3)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%smc3,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,4)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%smc4,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(1)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%stc1,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(2)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%stc2,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(3)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%stc3,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = clm2_struc(n)%clm(i)%t_soisno(4)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%stc4,temp)

  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,1)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%sh2o1,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,2)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%sh2o2,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,3)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%sh2o3,temp)
  do i=1,LIS_rc%ntiles(n)
     temp(i) = soilmoist(i,4)
  enddo
  call LIS_tile2grid(n,liswrf_export(n)%sh2o4,temp)

  do j=1,LIS_rc%lnr(n)
     do i=1,LIS_rc%lnc(n)
        k = LIS_domain(n)%gindex(i,j) 
        if(LIS_domain(n)%gindex(i,j).ne.-1) then 
           liswrf_export(n)%xland(i,j) = LIS_domain(n)%tile(k)%vegt
        endif
     enddo
  enddo
#endif
end subroutine clm2_setwrfexport
 









