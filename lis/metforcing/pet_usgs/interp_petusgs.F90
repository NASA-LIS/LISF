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
! !ROUTINE: interp_petusgs
! \label{interp_petusgs}
!
!
! !INTERFACE: 
subroutine interp_petusgs(n, findex, nx, ny, finput, lis_gds, nc, nr, varfield)

! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_logMod,         only : LIS_logunit
  use petusgs_forcingMod, only : petusgs_struc

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
  integer             :: nx
  integer             :: ny
  integer             :: nc
  integer             :: nr
  real, dimension(nx,ny) :: finput
  real, dimension(nc,nr) :: varfield
!
! !DESCRIPTION:
!   This subroutine interpolates a given USGS PET field 
!   to the LIS grid. 
!
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nx]
!  number of columns (in the east-west dimension) in the USGS PET grid
! \item[ny]
!  number of rows (in the north-south dimension) in the USGS PET grid
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
!  The routines invoked are: 
!  \begin{description}
!   \item[conserv\_interp](\ref{conserv_interp}) \newline
!     spatially interpolate the forcing data using conservative interpolation
!
! \end{description}
!EOP
  real                         :: lis_gds(50)
  logical*1, dimension(nx*ny)  :: lb
  logical*1, dimension(nc*nr)  :: lo
  real,      dimension(nx*ny)  :: f
  integer                      :: iret
  integer                      :: mo
  integer                      :: count1,i,j,v
  real,       dimension(nc*nr) :: lis1d

!=================  End variable declarations  ======================

   mo = nc*nr
   v = 0
   do j=1, ny
     do i=1, nx
        v = v+ 1
        f(v) = finput(i, j)
     end do
   end do

!--------------------------------------------------------------------
! Initialize output bitmap.
!--------------------------------------------------------------------
  lb = .true.
  lo = .true.

  petusgs_struc(n)%mi = nx * ny

  if( LIS_rc%met_interp(findex) == "bilinear" ) then       ! Bilinear interpolation

     call bilinear_interp( lis_gds,lb,f,lo,lis1d,petusgs_struc(n)%mi,mo, &
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          petusgs_struc(n)%w111, &
          petusgs_struc(n)%w121,petusgs_struc(n)%w211,petusgs_struc(n)%w221,    &
          petusgs_struc(n)%n111,petusgs_struc(n)%n121,petusgs_struc(n)%n211,    &
          petusgs_struc(n)%n221,LIS_rc%udef, iret )

  else if( LIS_rc%met_interp(findex) == "budget-bilinear" ) then  ! Conservative-budget

    call conserv_interp( lis_gds,lb,f,lo,lis1d,petusgs_struc(n)%mi,mo, &
         LIS_domain(n)%lat, LIS_domain(n)%lon,&
         petusgs_struc(n)%w111, &
         petusgs_struc(n)%w121,petusgs_struc(n)%w211,petusgs_struc(n)%w221,    &
         petusgs_struc(n)%n111,petusgs_struc(n)%n121,petusgs_struc(n)%n211,    &
         petusgs_struc(n)%n221,LIS_rc%udef, iret )

  else if( LIS_rc%met_interp(findex) == "neighbor" ) then  ! Nearest neighbor

    call neighbor_interp( lis_gds,lb,f,lo,lis1d,petusgs_struc(n)%mi,mo, &
         LIS_domain(n)%lat, LIS_domain(n)%lon,&
         petusgs_struc(n)%n111,    &
         LIS_rc%udef, iret )

  end if

!--------------------------------------------------------------------
! Create 2D array for main program.  
!--------------------------------------------------------------------
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_petusgs

