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
! !ROUTINE: cable_setwrfexport.F90
!
! !REVISION HISTORY:
!  02 Dec 2003; Sujay Kumar, Initial Version
!  17 Nov 2008; Sujay Kumar, Modified for the ESMF coupled version
!  17 Oct 2011: Claire Carouge (ccc), CABLE LSM improvements
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
! 
! !INTERFACE:
subroutine cable_setwrfexport(n, expState)
! !USES:
  use ESMF
  use LIS_precisionMod
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_logMod,     only : LIS_verify
  use LISWRFGridCompMod, only : LISWRF_export
  use cable_lsmMod
  use cable_dimensions, only : ms

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
  type(ESMF_State)      :: expState
! !DESCRIPTION:  
!  Defines the export states from CABLE to WRF in coupled mode
!EOP
  integer               :: i,j,k,t,m
  real, allocatable :: temp(:)
  real, allocatable :: soilmoist(:,:)
  real, allocatable :: q1(:)
  real :: RCH

  allocate(temp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(soilmoist(LIS_rc%npatch(n,LIS_rc%lsm_index),1:4))
  allocate(q1(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%avgsurft_t,cable_struc(n)%cable%radt)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qh_t,cable_struc(n)%cable%qh)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%eta_kinematic_t,cable_struc(n)%cable%qle/2.501E6)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qle_t,cable_struc(n)%cable%qle)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qg_t,-cable_struc(n)%cable%qg)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%chs2_t,cable_struc(n)%cable%forc_chs2)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%cqs2_t,cable_struc(n)%cable%forc_cqs2)
! Not used by WRF, only for output. Not currently calculated by CABLE.
! Left for reference for future use. ccc
!----------------------------------------------------------------------
! Root Zone Soil Moisture (kg m-2)
! Calculation of root zone soil moisture 
!----------------------------------------------------------------------
!  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!
!     soilmr(t)=0.
!     do m=1,nlevsoi
!        soilmr(t)=soilmr(t)+clm2_struc(n)%clm(t)%rootfr(m)*clm_struc(n)%clm(t)%h2osoi_liq(m)
!     enddo
!  enddo
!  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%rootmoist,soilmr) 

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%soilm_t,cable_struc(n)%cable%wbtot)  !mm
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qs_t,cable_struc(n)%cable%qs)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qsb_t,cable_struc(n)%cable%qsb)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%albedo_t,cable_struc(n)%cable%albedo)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%znt_t,cable_struc(n)%cable%z0m)

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%snocvr_t,real(cable_struc(n)%cable%isflag))
!mixing  ratio
  do k=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(cable_struc(n)%cable(k)%qair.gt.1E-19) then  
        RCH = cable_struc(n)%cable(k)%psurf/(287.04*cable_struc(n)%cable(k)%tair*(1+0.61*cable_struc(n)%cable(k)%qair)) 
        q1(k) = cable_struc(n)%cable(k)%qair+(cable_struc(n)%cable(k)%qle/2.501E6)/(RCH*cable_struc(n)%cable(k)%forc_ch)
     else
        q1(k) = 0.0
     endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%q1_t,q1)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%snow_t,cable_struc(n)%cable%snowd)
!canopy interception
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%cmc_t,cable_struc(n)%cable%cansto)
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%qsm_T,cable_struc(n)%cable%smelt)

  ! Need soil moisture in kg m-2
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilmoist(i,1) = cable_struc(n)%cable(i)%wb(1)*1000.*cable_struc(n)%cable(i)%zse(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc1_t,soilmoist(:,1))
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilmoist(i,2) = cable_struc(n)%cable(i)%wb(2)*1000.*cable_struc(n)%cable(i)%zse(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc2_t,soilmoist(:,2))
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilmoist(i,3) = cable_struc(n)%cable(i)%wb(3)*1000.*cable_struc(n)%cable(i)%zse(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc3_t,soilmoist(:,3))
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilmoist(i,4) = cable_struc(n)%cable(i)%wb(4)*1000.*cable_struc(n)%cable(i)%zse(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%smc4_t,soilmoist(:,4))

  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc1_t,cable_struc(n)%cable%tgg(1))
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc2_t,cable_struc(n)%cable%tgg(2))
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc3_t,cable_struc(n)%cable%tgg(3))
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%stc4_t,cable_struc(n)%cable%tgg(4))

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,1)-cable_struc(n)%cable(i)%wbice(1)*1000.*cable_struc(n)%cable(i)%zse(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,2)-cable_struc(n)%cable(i)%wbice(2)*1000.*cable_struc(n)%cable(i)%zse(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,3)-cable_struc(n)%cable(i)%wbice(3)*1000.*cable_struc(n)%cable(i)%zse(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = soilmoist(i,4)-cable_struc(n)%cable(i)%wbice(4)*1000.*cable_struc(n)%cable(i)%zse(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%sh2o4_t,temp)

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index) 
     temp(i) = real(cable_struc(n)%cable(i)%vegtype)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,liswrf_export(n)%xland_t,temp)

  deallocate(temp)
  deallocate(soilmoist)
  deallocate(q1)

end subroutine cable_setwrfexport
 
