!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: interp_pildas
! \label{interp_pildas}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 14 Jul 2016: Mahdi Navari - Modified for PILDAS
!
! !INTERFACE:
subroutine interp_pildas(n,findex,input_data,output_2d)
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_spatialDownscalingMod
  use pildas_forcingMod, only :pildas_struc
  
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: findex
  real, intent(in)      :: input_data(pildas_struc(n)%nc,pildas_struc(n)%nr)
  real, intent(inout)   :: output_2d(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:
!   This subroutine interpolates a given GEOS5 data field 
!   to the LIS grid. 
!EOP
  integer :: iret
  integer :: mo
  integer :: count1,i,j,c,r,t
  real    :: output_data(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical*1 :: output_bitmap(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical*1 :: input_bitmap(pildas_struc(n)%nc*pildas_struc(n)%nr)
  real      :: f(pildas_struc(n)%nc*pildas_struc(n)%nr)
  
  input_bitmap = .false.
  do r=1,pildas_struc(n)%nr
     do c=1,pildas_struc(n)%nc
        t = c+(r-1)*pildas_struc(n)%nc
        f(t) = input_data(c,r)
        if(f(t).gt.-1E10) then 
           input_bitmap(t) = .true. 
        endif
     enddo
  end do
  
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

!-----------------------------------------------------------------------
! Initialize output bitmap. 
!-----------------------------------------------------------------------
  output_bitmap = .true.

!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------
  if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then 
     call bilinear_interp(LIS_rc%gridDesc(n,:),input_bitmap,f,&
          output_bitmap,&
          output_data,pildas_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          pildas_struc(n)%w111, pildas_struc(n)%w121,&
          pildas_struc(n)%w211,pildas_struc(n)%w221,&
          pildas_struc(n)%n111,pildas_struc(n)%n121,&
          pildas_struc(n)%n211,pildas_struc(n)%n221,LIS_rc%udef,iret)
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then 
     call neighbor_interp(LIS_rc%gridDesc(n,:),input_bitmap,&
          f, output_bitmap, output_data, pildas_struc(n)%mi, mo, & 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          pildas_struc(n)%n113, LIS_rc%udef, iret)
  endif

!-----------------------------------------------------------------------    
! Transform to a 2D array.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, LIS_rc%lnr(n)
     do i = 1, LIS_rc%lnc(n)
        output_2d(i,j) = output_data(i+count1)
     enddo
     count1 = count1 + LIS_rc%lnc(n)
  enddo

end subroutine interp_pildas
