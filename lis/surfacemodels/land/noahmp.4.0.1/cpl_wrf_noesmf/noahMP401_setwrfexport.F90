!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: noahMP401_setwrfexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from Noah to WRF in coupled mode
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 17 Nov 2008; Sujay Kumar, Modified for the ESMF coupled version
! 28 Jun 2018; Chandana Gangodagamage Modified for the NoahMP.3.6
! 
! !INTERFACE:
subroutine noahMP401_setwrfexport(n)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_logMod
  use LISWRFGridCompMod, only : LISWRF_export
!  use noah33_lsmMod
  use NoahMP401_lsmMod
  
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
!EOP
  integer               :: i,j,k,t
  real, allocatable     :: temp(:)

  allocate(temp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%avgsurft_t,&
!       NOAHMP401_struc(n)%noahmp401%t1)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%qh_t,&
!       NOAHMP401_struc(n)%noahmp401%qh)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%eta_kinematic_t,&
!       noah33_struc(n)%noah%eta_kinematic)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%qle_t,&
!       noah33_struc(n)%noah%qle)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%qg_t,&
!       noah33_struc(n)%noah%qg)
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%albedo_t,&
       NOAHMP401_struc(n)%noahmp401%albedo)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%znt_t,&
!       noah33_struc(n)%noah%z0)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%q1_t,&
!       noah33_struc(n)%noah%q1)
#ifdef WRF_HYDRO
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%infxsrt_t,&
       NOAHMP401_struc(n)%noahmp401%infxs1rt)
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%soldrain_t,&
       NOAHMP401_struc(n)%noahmp401%soldrain1rt)
#endif

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%smc(1)
     !  print*, 'LISsideout:', noah33_struc(n)%noah(i)%chs2, noah33_struc(n)%noah(i)%cqs2
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc1_t,&
       temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%smc(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc2_t,&
       temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%smc(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc3_t,&
       temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%smc(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc4_t,&
       temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%tslb(NOAHMP401_struc(n)%nsnow+1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%tslb(NOAHMP401_struc(n)%nsnow+2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%tslb(NOAHMP401_struc(n)%nsnow+3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%tslb(NOAHMP401_struc(n)%nsnow+4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc4_t,temp)

!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%chs2_t,&
!       noah33_struc(n)%noah%chs2)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%cqs2_t,&
!       noah33_struc(n)%noah%cqs2)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%snocvr_t,&
!       noah33_struc(n)%noah%sca)
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%snow_t,&
       NOAHMP401_struc(n)%noahmp401%sneqv*1000.0)
  ! NUWRF EMK
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%snowh_t,&
       NOAHMP401_struc(n)%noahmp401%snowh)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%lispor_t,&
!       noah33_struc(n)%noah%smcmax)


!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%rootmoist_t,&
!       noah33_struc(n)%noah%rootmoist) !units? 
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%soilm_t,&
!       noah33_struc(n)%noah%soilm)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%qs_t,&
!       noah33_struc(n)%noah%qs)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%qsb_t,&
!       noah33_struc(n)%noah%qsb)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%cmc_t,&
!       noah33_struc(n)%noah%cmc)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%qsm_t,&
!       noah33_struc(n)%noah%qsm)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%emiss_t,&
!       noah33_struc(n)%noah%emiss)
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%xice_t,&
!       noah33_struc(n)%noah%xice)

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%sh2o(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%sh2o(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%sh2o(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index) 
     temp(i) = NOAHMP401_struc(n)%noahmp401(i)%sh2o(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o4_t,temp)

!  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     temp(i) = noah33_struc(n)%noah(i)%relsmc(1)
!  enddo
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%relsmc1_t,temp)
!  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     temp(i) = noah33_struc(n)%noah(i)%relsmc(2)
!  enddo
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%relsmc2_t,temp)
!  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     temp(i) = noah33_struc(n)%noah(i)%relsmc(3)
!  enddo
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%relsmc3_t,temp)
!  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index) 
!     temp(i) = noah33_struc(n)%noah(i)%relsmc(4)
!  enddo
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%relsmc4_t,temp)

!  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index) 
!     temp(i) = real(noah33_struc(n)%noah(i)%vegt)
!  enddo
!  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%xland_t,temp)

  deallocate(temp)
  !  do j=1,LIS_rc%lnr(n)
  !     do i=1,LIS_rc%lnc(n)
  !        k = LIS_domain(n)%gindex(i,j) 
  !        if(LIS_domain(n)%gindex(i,j).ne.-1) then 
  !           LISWRF_export(n)%xland(i,j) = noah33_struc(n)%noah(k)%vegt
  !        endif
  !     enddo
  !  enddo

end subroutine noahMP401_setwrfexport
 
