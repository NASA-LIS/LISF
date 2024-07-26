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
! !ROUTINE: AGRMET_interp_cmorph
! \label{AGRMET_interp_cmorph}
! !REVISION HISTORY:
! 05 May 2013; Moved to AGRMET from suppforcing...Ryan Ruhge/16WS/WXE/SEMS
!
!
! !INTERFACE: 
subroutine AGRMET_interp_cmorph(n, nx, ny, finput, lis_gds, nc, nr, varfield)
! !USES:
  
  use LIS_coreMod,       only : LIS_domain
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer             :: nx
  integer             :: ny
  integer             :: nc
  integer             :: nr
  real, dimension(nx,ny) :: finput
  real, dimension(nc,nr) :: varfield
!
! !DESCRIPTION:
!   This subroutine interpolates a given CMORPH field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[nx]
!  number of columns (in the east-west dimension) in the CMORPH grid
! \item[ny]
!  number of rows (in the north-south dimension) in the CMORPH grid
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
  real :: lis_gds(50)
  logical*1, dimension(nx*ny)  :: lb
  logical*1, dimension(nc*nr)  :: lo
  real, dimension(nx*ny) :: f
  integer             :: ngdas

  integer :: ip, ipopt(20),km,iret
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
     ip=3
     ipopt(1)=-1
     ipopt(2)=-1
     km=1
   f = -9999.0
   lb = .false.
   v = 0
   Do j=1, ny
     Do i=1, nx
      v = v+ 1
      if (finput(i,j) .ne. -9999.0) then
        f(v) = finput(i, j)
        lb(v) = .true.
      endif
     End Do
   End Do

!--------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!--------------------------------------------------------------------
  lo = .true.



  agrmet_struc(n)%micmor = ngdas

  call conserv_interp(lis_gds,lb,f,lo,lis1d,&
       agrmet_struc(n)%micmor,mo,&
       LIS_domain(n)%lat, LIS_domain(n)%lon, &
       agrmet_struc(n)%w112cmor, agrmet_struc(n)%w122cmor,&
       agrmet_struc(n)%w212cmor,agrmet_struc(n)%w222cmor,&
       agrmet_struc(n)%n112cmor,agrmet_struc(n)%n122cmor,&
       agrmet_struc(n)%n212cmor,&
       agrmet_struc(n)%n222cmor,-9999.0, iret)


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

end subroutine AGRMET_interp_cmorph

