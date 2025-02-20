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
! !ROUTINE: noahMP36_setwrfexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from Noah to WRF in coupled mode.
!   surface albedo (albedo)
!   soil moisture content (smc 1:4)
!   soil temperature (stc 1:4)
!   snow water equivalent (snow)
!   snow height (snowh)
!   volumetric liquid soil moisture (sh2o 1:4)
!   if WRF_HYDRO
!    infiltration excess (infxsrt)
!    soil drainage (soldrain)
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 17 Nov 2008; Sujay Kumar, Modified for the ESMF coupled version
! 28 Jun 2018; Chandana Gangodagamage Modified for the NoahMP.3.6
! 19 Jan 2021; Dan Rosen, Remove commented code blocks
! 
! !INTERFACE:
subroutine noahMP36_setwrfexport(n)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_logMod
  use LISWRFGridCompMod, only : LISWRF_export
!  use noah33_lsmMod
  use noahMP36_lsmMod
  
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
!EOP
  integer               :: i,j,k,t
  real, allocatable     :: temp(:)

  ! surface albedo
  allocate(temp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%albedo_t,&
       NOAHMP36_struc(n)%noahmp36%albedo)
#ifdef WRF_HYDRO
  ! infiltration excess
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%infxsrt_t,&
       NOAHMP36_struc(n)%noahmp36%infxs1rt)
  ! soil drainage
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%soldrain_t,&
       NOAHMP36_struc(n)%noahmp36%soldrain1rt)
#endif
  ! soil moisture content layers 1:4
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%smc(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc1_t,&
       temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%smc(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc2_t,&
       temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%smc(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc3_t,&
       temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%smc(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%smc4_t,&
       temp)
  ! soil temperature layers 1:4
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sstc(NOAHMP36_struc(n)%nsnow+1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sstc(NOAHMP36_struc(n)%nsnow+2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sstc(NOAHMP36_struc(n)%nsnow+3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sstc(NOAHMP36_struc(n)%nsnow+4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%stc4_t,temp)
  ! snow water equivalent
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%snow_t,&
       NOAHMP36_struc(n)%noahmp36%sneqv*1000.0)
  ! snow height NUWRF EMK 
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%snowh_t,&
       NOAHMP36_struc(n)%noahmp36%snowh)
  ! volumetric liquid soil moisture layers 1:4
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sh2o(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sh2o(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sh2o(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index) 
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%sh2o(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%sh2o4_t,temp)

#ifdef PARFLOW
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%wtrflx(1)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%wtrflx1_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%wtrflx(2)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%wtrflx2_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%wtrflx(3)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%wtrflx3_t,temp)
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     temp(i) = NOAHMP36_struc(n)%noahmp36(i)%wtrflx(4)
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,LISWRF_export(n)%wtrflx4_t,temp)
#endif

  deallocate(temp)

end subroutine noahMP36_setwrfexport

