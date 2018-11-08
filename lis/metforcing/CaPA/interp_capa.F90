!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: interp_capa
! \label{interp_capa}
!
! !INTERFACE: 
subroutine interp_capa(n, ncapa, f,lb, lis_gds, nc, nr, varfield)
! !USES:  
  use LIS_coreMod,     only : LIS_domain
  use capa_forcingMod, only : capa_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer                :: nc
  integer                :: nr
  integer                :: ncapa
  real                   :: lis_gds(50)
  real                   :: f(ncapa)
  logical*1              :: lb(ncapa)
  real, dimension(nc,nr) :: varfield
!
! !DESCRIPTION:
!   This subroutine interpolates a given CAPA field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[ncapa]
!  number of elements in the input grid
! \item[f]
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
!
!  The routines invoked are: 
!  \begin{description}
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: mo
  integer :: count1,i,j

  real, dimension(nc*nr) :: lis1d

  logical*1 :: lo(nc*nr)

!=== End variable declarations

  mo = nc*nr

!--------------------------------------------------------------------  
! Initialize output bitmap.
!--------------------------------------------------------------------  
  lo = .true.
  

  capa_struc(n)%mi = ncapa
  call conserv_interp(lis_gds,lb,f,lo,lis1d,capa_struc(n)%mi,mo,&
       LIS_domain(n)%lat, LIS_domain(n)%lon,&
       capa_struc(n)%w112, & ! EMK bug fix
       capa_struc(n)%w122,capa_struc(n)%w212,capa_struc(n)%w222,&
       capa_struc(n)%n112,capa_struc(n)%n122,capa_struc(n)%n212,&
       capa_struc(n)%n222,-9999.0,iret)

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

end subroutine interp_capa

