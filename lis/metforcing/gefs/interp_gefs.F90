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
! !ROUTINE: interp_gefs
! \label{interp_gefs}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
!
! !INTERFACE:
subroutine interp_gefs( n, findex, input_data, pcp_flag, output_2d )
! !USES:
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_spatialDownscalingMod
  use gefs_forcingMod, only : gefs_struc
  
  implicit none

! !ARGUMENTS:   
  integer, intent(in)   :: n 
  integer, intent(in)   :: findex
  real, intent(in)      :: input_data(gefs_struc(n)%nc*gefs_struc(n)%nr)
  logical, intent(in)   :: pcp_flag
  real, intent(inout)   :: output_2d(LIS_rc%lnc(n),LIS_rc%lnr(n))
!
! !DESCRIPTION:
!   This subroutine interpolates a given GEFS data field 
!   to the LIS grid. 
!EOP
  integer   :: count1,i,j,c,r,t
  integer   :: iret
  integer   :: ndata
  integer   :: mo
  logical*1 :: input_bitmap(gefs_struc(n)%nc*gefs_struc(n)%nr)
  logical*1 :: output_bitmap(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real      :: f(gefs_struc(n)%nc*gefs_struc(n)%nr)
  real      :: output_data(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

  input_bitmap = .false.
  do t = 1, gefs_struc(n)%mi
     if(input_data(t).lt.1E15) then 
        input_bitmap(t) = .true. 
     endif
  end do
  
!-----------------------------------------------------------------------
! Initialize output bitmap. 
!-----------------------------------------------------------------------
  output_bitmap = .true.

!-----------------------------------------------------------------------  
! Interpolate to LIS grid
!-----------------------------------------------------------------------
  if( LIS_rc%met_interp(findex).eq."bilinear") then 
      call bilinear_interp(LIS_rc%gridDesc(n,:), &
                    input_bitmap, input_data,    &
                    output_bitmap, output_data,  &
                    gefs_struc(n)%mi, mo,        &
                    LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                    gefs_struc(n)%w111, gefs_struc(n)%w121,&
                    gefs_struc(n)%w211, gefs_struc(n)%w221,&
                    gefs_struc(n)%n111, gefs_struc(n)%n121,&
                    gefs_struc(n)%n211, gefs_struc(n)%n221,&
                    LIS_rc%udef, iret)


  elseif( LIS_rc%met_interp(findex).eq."budget-bilinear") then 

    if( pcp_flag ) then  ! Budget-bilinear for precip fields
      call conserv_interp(LIS_rc%gridDesc(n,:), &
                    input_bitmap, input_data,   &
                    output_bitmap, output_data, &
                    gefs_struc(n)%mi, mo, &
                    LIS_domain(n)%lat, LIS_domain(n)%lon, &
                    gefs_struc(n)%w112,gefs_struc(n)%w122,&
                    gefs_struc(n)%w212,gefs_struc(n)%w222,&
                    gefs_struc(n)%n112,gefs_struc(n)%n122,&
                    gefs_struc(n)%n212,gefs_struc(n)%n222,&
                    LIS_rc%udef, iret)
    else
      call bilinear_interp(LIS_rc%gridDesc(n,:), &
                    input_bitmap, input_data,    &
                    output_bitmap, output_data,  &
                    gefs_struc(n)%mi, mo,        &
                    LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                    gefs_struc(n)%w111, gefs_struc(n)%w121,&
                    gefs_struc(n)%w211, gefs_struc(n)%w221,&
                    gefs_struc(n)%n111, gefs_struc(n)%n121,&
                    gefs_struc(n)%n211, gefs_struc(n)%n221,&
                    LIS_rc%udef, iret)
    endif
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

end subroutine interp_gefs
