!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: interp_rdhm356_temper
! \label{interp_rdhm356_temper}
!
! !REVISION HISTORY:
! 03 May 2010: Soni Yatheendradas; Precip and Temper. input grids now can have different
!                                  extents/directories and different from the run-domain 
!                                  extent, as per the new input grids posted onto the 
!                                  DMIP2 website for Sierra Nevada
! 18 Dec 2013: Shugong Wang; implementation for RDHM356
!
! !INTERFACE: 
subroutine interp_rdhm356_temper (n,findex, ksec1, nrdhm356, &
                                   ppt_field, lb, lis_gds, &
                                    nc, nr, varfield )

! !USES:  
  use rdhm356_forcingMod, only : rdhm356_struc_temper
  use LIS_coreMod, only: LIS_rc, LIS_domain

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer                :: nc      ! Number of columns (in the E-W dimension) in the LIS grid
  integer                :: nr      ! Number of rows (in the N-S dimension) in the LIS grid
  integer                :: nrdhm356   ! Number of points in original dmip II grid
  integer                :: ksec1(100)
  real                   :: lis_gds(50) 
  real                   :: ppt_field(nrdhm356) 
  logical*1              :: lb(nrdhm356)
  real, dimension(nc,nr) :: varfield 

! !DESCRIPTION:
!   This subroutine interpolates a given STG2 field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[ksec1]
!  grib decoding array
! \item[nrdhm356]
!  number of elements in the input grid
! \item[ppt\_field]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    Spatially interpolate the forcing data using bilinear interpolation, or
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative (budget bil.) interpolation
! \end{description}
!EOP

  integer :: iret
  integer :: mo
  integer :: count1, i, j

  real    :: lis1d(nc*nr)   ! LIS Grid  
  logical*1 :: lo(nc*nr)    ! LIS Grid


!--------------------------------------------------------------------
! NOTE:: Recommended to use budget bilinear for precip forcing fields
!--------------------------------------------------------------------
  rdhm356_struc_temper(n)%mi = nrdhm356
  mo = nc * nr   ! LIS Grid

  lo = .true.    !-Initialize output bitmap
  lb = .true.   

  do i=1, nc*nr
    if(ppt_field(i) == LIS_rc%udef) then
      lb(i) = .false. 
    endif
  enddo

  if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then 
     call bilinear_interp( lis_gds, lb, ppt_field, lo, &
                          lis1d, rdhm356_struc_temper(n)%mi, mo, &
                          LIS_domain(n)%lat, LIS_domain(n)%lon,&
                          rdhm356_struc_temper(n)%w111, rdhm356_struc_temper(n)%w121, &
                          rdhm356_struc_temper(n)%w211, rdhm356_struc_temper(n)%w221, &
                          rdhm356_struc_temper(n)%n111, rdhm356_struc_temper(n)%n121, &
                          rdhm356_struc_temper(n)%n211, rdhm356_struc_temper(n)%n221, &
                          LIS_rc%udef, iret)
  elseif (trim(LIS_rc%met_interp(findex)) .eq. "budget-bilinear") then 
     call conserv_interp( lis_gds, lb, ppt_field, lo, &
                          lis1d, rdhm356_struc_temper(n)%mi, mo, & 
                          LIS_domain(n)%lat, LIS_domain(n)%lon,&
                          rdhm356_struc_temper(n)%w112, rdhm356_struc_temper(n)%w122,&
                          rdhm356_struc_temper(n)%w212, rdhm356_struc_temper(n)%w222,&
                          rdhm356_struc_temper(n)%n112, rdhm356_struc_temper(n)%n122,&
                          rdhm356_struc_temper(n)%n212, rdhm356_struc_temper(n)%n222,&
                          LIS_rc%udef,iret)
  elseif( trim(LIS_rc%met_interp(findex)) == "neighbor" ) then 
     call neighbor_interp(lis_gds,lb,ppt_field,lo,lis1d,&
                          rdhm356_struc_temper(n)%mi,mo,&
                          LIS_domain(n)%lat, LIS_domain(n)%lon,&
                          rdhm356_struc_temper(n)%n113, LIS_rc%udef, iret)
  endif

!--------------------------------------------------------------------    
! Create 2D array (on LIS LIS_domain) for main program. 
!--------------------------------------------------------------------    
   count1 = 0
   do j = 1, nr
      do i = 1, nc 
         varfield(i,j) = lis1d(i+count1)
      enddo
      count1 = count1 + nc
   enddo


end subroutine interp_rdhm356_temper

