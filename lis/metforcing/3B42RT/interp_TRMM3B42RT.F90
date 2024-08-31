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
! !ROUTINE: interp_TRMM3B42RT
! \label{interp_TRMM3B42RT}
!
! !INTERFACE: 
subroutine interp_TRMM3B42RT(n, nx, ny, finput, lis_gds, nc, nr, varfield, findex)
! !USES:
  use TRMM3B42RT_forcingMod, only : TRMM3B42RT_struc
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod,    only : LIS_logunit, LIS_endrun 

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n 
  integer             :: nx
  integer             :: ny
  integer             :: nc
  integer             :: nr
  real                :: lis_gds(50)
  real, dimension(nx,ny) :: finput
  integer, intent(in) :: findex ! SY

!
! !DESCRIPTION:
!   This subroutine interpolates a given TRMM 3B42RT field
!   to the LIS grid.
!  The arguments are:
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nx]
!  number of columns (in the east-west dimension) in the TRMM 3B42RT grid
! \item[ny]
!  number of rows (in the north-south dimension) in the TRMM 3B42RT grid
! \item[finput]
!  input data array to be interpolated
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
!
!  The routines invoked are:
!  \begin{description}
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP
  logical*1, dimension(nx*ny)  :: lb
  logical*1, dimension(nc*nr)  :: lo
  real, dimension(nc,nr) :: varfield

  real, dimension(nx*ny) :: f
  integer :: ngdas
  integer :: iret
  integer :: mo
  integer :: count1,i,j,v

  real, dimension(nc*nr) :: lis1d


!=== End variable declarations
!--------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!--------------------------------------------------------------------
  ngdas = nx * ny
  mo = nc*nr
   v = 0
   Do j=1, ny
     Do i=1, nx
      v = v+ 1
      f(v) = finput(i, j) 
     End Do
   End Do
  
!--------------------------------------------------------------------  
! Initialize output bitmap. Important for soil moisture and temp.
!--------------------------------------------------------------------  
  lb = .true.
  lo = .true.

  where ( f < 0 )
     lb = .false.
  endwhere
  

  TRMM3B42RT_struc(n)%mi = ngdas

  if(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then ! SY
   call conserv_interp(lis_gds,lb,f,lo,lis1d,TRMM3B42RT_struc(n)%mi,mo,&
        LIS_domain(n)%lat, LIS_domain(n)%lon, &
        TRMM3B42RT_struc(n)%w112,&
        TRMM3B42RT_struc(n)%w122,TRMM3B42RT_struc(n)%w212,&
        TRMM3B42RT_struc(n)%w222,&
        TRMM3B42RT_struc(n)%n112,TRMM3B42RT_struc(n)%n122,&
        TRMM3B42RT_struc(n)%n212,&
        TRMM3B42RT_struc(n)%n222,-9999.0,iret)
  elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then ! SY
  ! SY: begin for nearest neighbor
   call neighbor_interp(lis_gds,lb,f,lo,lis1d,TRMM3B42RT_struc(n)%mi,mo,&
        LIS_domain(n)%lat, LIS_domain(n)%lon, &
        TRMM3B42RT_struc(n)%n112, LIS_rc%udef,iret)
  ! SY: end for nearest neighbor
  else
   write(LIS_logunit,*) 'This interpolation not defined for TRMM data'
   write(LIS_logunit,*) 'Program stopping ... '
   call LIS_endrun()
  endif

!--------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between GDAS & LDAS. For LDAS land 
! points not included in GDAS geography dataset only.
!--------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo
end subroutine interp_TRMM3B42RT

