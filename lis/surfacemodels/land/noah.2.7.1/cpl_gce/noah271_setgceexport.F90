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
! !ROUTINE: noah271_setgceexport.F90
!
! !DESCRIPTION:  
!  Defines the export states from Noah to GCE in coupled mode
!
! !REVISION HISTORY:
! 02 Dec 2003; Sujay Kumar, Initial Version
! 
! !INTERFACE:
subroutine noah271_setgceexport(n, expState)
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc,LIS_domain
  use LIS_historyMod, only : LIS_tile2grid
  use lisGCEGridCompMod, only : lisgce_export
  use noah271_lsmMod
  
  implicit none
  integer, intent(in) :: n
  type(ESMF_State)    :: expState
 
  integer :: i
  real :: temp(LIS_rc%ntiles(n))
  real :: t2v, rho
  real :: temp1(LIS_rc%lnc(n), LIS_rc%lnr(n))
!EOP
! to avoid the real*8 -> real*4 problem in the tile2grid call
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%swnet)
  lisgce_export%swnet = temp1
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%lwnet)
  lisgce_export%lwnet = temp1

  do i=1,LIS_rc%ntiles(n)
     t2v = noah271_struc(n)%noah(i)%tair*(1+0.61*&
          noah271_struc(n)%noah(i)%qair)
     rho = noah271_struc(n)%noah(i)%psurf/(287.04*t2v)
     temp(i) = noah271_struc(n)%noah(i)%qle/(2.5E6*rho)
  enddo 
  call LIS_tile2grid(n,temp1,temp)
  lisgce_export%qle = temp1

  do i=1,LIS_rc%ntiles(n)
     t2v = noah271_struc(n)%noah(i)%tair*(1+0.61*&
          noah271_struc(n)%noah(i)%qair)
     rho = noah271_struc(n)%noah(i)%psurf/(287.04*t2v)
     temp(i) = noah271_struc(n)%noah(i)%qh/(1005.0*rho)
  enddo 
  call LIS_tile2grid(n,temp1,temp)
  lisgce_export%qh = temp1
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%qg)
  lisgce_export%qg = temp1
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%t1)
  lisgce_export%avgsurft = temp1
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%albedo)
  lisgce_export%albedo = temp1
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%tauu)
  lisgce_export%tauu = temp1
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%tauv)
  lisgce_export%tauv = temp1
  call LIS_tile2grid(n,temp1,noah271_struc(n)%noah%tau)
  lisgce_export%tau = temp1
end subroutine noah271_setgceexport
 









