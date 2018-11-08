!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: interp_rdhm356_precip
! \label{interp_rdhm356_precip}
!
! !REVISION HISTORY:
! 03 May 2010: Soni Yatheendradas; Precip and Temper. input grids now can have different
!                                  extents/directories and different from the run-domain 
!                                  extent, as per the new input grids posted onto the 
!                                  DMIP2 website for Sierra Nevada
! 18 Dec 2013: Shugong Wang; implementation for RDHM356
!
! !INTERFACE: 
subroutine interp_rdhm356_precip (n,findex, ksec1, nrdhm356, &
                                   ppt_field, lb, lis_gds, &
                                   nc, nr, varfield )

! !USES:  
  use rdhm356_forcingMod, only : rdhm356_struc_precip
  use LDT_coreMod, only: LDT_rc, LDT_domain

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer                :: nc      ! Number of columns (in the E-W dimension) in the LDT grid
  integer                :: nr      ! Number of rows (in the N-S dimension) in the LDT grid
  integer                :: nrdhm356   ! Number of points in original dmip II grid
  integer                :: ksec1(100)
  real                   :: lis_gds(20) 
  real                   :: ppt_field(nrdhm356) 
  logical*1              :: lb(nrdhm356)
  real, dimension(nc,nr) :: varfield 

! !DESCRIPTION:
!   This subroutine interpolates a given STG2 field 
!   to the LDT grid. 
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
!  array description of the LDT grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LDT grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LDT grid
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

  real    :: lis1d(nc*nr)   ! LDT Grid  
  logical*1 :: lo(nc*nr)    ! LDT Grid


!--------------------------------------------------------------------
! NOTE:: Recommended to use budget bilinear for precip forcing fields
!--------------------------------------------------------------------
  rdhm356_struc_precip(n)%mi = nrdhm356
  mo = nc * nr   ! LDT Grid

  lo = .true.    !-Initialize output bitmap
  lb = .true.   
  
  do i=1, nc*nr
    if(ppt_field(i) == LDT_rc%udef) then
      lb(i) = .false. 
    endif
  enddo

  if (trim(LDT_rc%met_gridtransform(findex)) .eq. "bilinear") then 
     call bilinear_interp( lis_gds, lb, ppt_field, lo, &
                          lis1d, rdhm356_struc_precip(n)%mi, mo, &
                          LDT_domain(n)%lat, LDT_domain(n)%lon,&
                          rdhm356_struc_precip(n)%w111, rdhm356_struc_precip(n)%w121, &
                          rdhm356_struc_precip(n)%w211, rdhm356_struc_precip(n)%w221, &
                          rdhm356_struc_precip(n)%n111, rdhm356_struc_precip(n)%n121, &
                          rdhm356_struc_precip(n)%n211, rdhm356_struc_precip(n)%n221, &
                          LDT_rc%udef, iret)
  elseif (trim(LDT_rc%met_gridtransform(findex)) .eq. "budget-bilinear") then 
     call conserv_interp( lis_gds, lb, ppt_field, lo, &
                          lis1d, rdhm356_struc_precip(n)%mi, mo, & 
                          LDT_domain(n)%lat, LDT_domain(n)%lon,&
                          rdhm356_struc_precip(n)%w112, rdhm356_struc_precip(n)%w122,&
                          rdhm356_struc_precip(n)%w212, rdhm356_struc_precip(n)%w222,&
                          rdhm356_struc_precip(n)%n112, rdhm356_struc_precip(n)%n122,&
                          rdhm356_struc_precip(n)%n212, rdhm356_struc_precip(n)%n222,&
                          LDT_rc%udef,iret)
  elseif( trim(LDT_rc%met_gridtransform(findex)) == "neighbor" ) then
     call neighbor_interp(lis_gds,lb,ppt_field,lo,lis1d,&
                          rdhm356_struc_precip(n)%mi,mo,&
                          LDT_domain(n)%lat, LDT_domain(n)%lon,&
                          rdhm356_struc_precip(n)%n113, LDT_rc%udef, iret)
  endif

!--------------------------------------------------------------------    
! Create 2D array (on LDT LDT_domain) for main program. 
!--------------------------------------------------------------------    
   count1 = 0
   do j = 1, nr
      do i = 1, nc 
         varfield(i,j) = lis1d(i+count1)
      enddo
      count1 = count1 + nc
   enddo


end subroutine interp_rdhm356_precip

