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
! !ROUTINE: interp_cmap
! \label{interp_cmap}
! !REVISION HISTORY: 
!  20 Jan 2006:  Yudong Tian: Initial Implementation
!  10 Oct 2006:  Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
!  10 Oct 2014:  KR Arsenault: Added to LDT
!
! !INTERFACE: 
subroutine interp_cmap(n, ncmap, f, lb, ldt_gds, nc, nr, &
                       varfield)
! !USES:  
  use cmap_forcingMod, only : cmap_struc
  use LDT_coreMod,     only : LDT_rc, LDT_domain

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer                :: nc
  integer                :: nr
  integer                :: ncmap
  real                   :: ldt_gds(20)
  real                   :: f(ncmap)
  logical*1              :: lb(ncmap)
  real, dimension(nc,nr) :: varfield
!
! !DESCRIPTION:
!   This subroutine interpolates a given CMAP field 
!   to the LDT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!   index of the nest
! \item[ncmap]
!   number of elements in the input grid
! \item[f]
!   input data array to be interpolated
! \item[lb]
!   input bitmap
! \item[ldt\_gds]
!   array description of the LDT grid
! \item[nc]
!   number of columns (in the east-west dimension) in the LDT grid
! \item[nr]
!   number of rows (in the north-south dimension) in the LDT grid
! \item[varfield]
!   output interpolated field
! \end{description} 
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

  real, dimension(nc*nr) :: ldt1d

  logical*1 :: lo(nc*nr)

!=== End variable declarations

!--------------------------------------------------------------------  
! Initialize output bitmap. 
!--------------------------------------------------------------------  
  lo = .true.
  mo = nc*nr
  cmap_struc(n)%mi = ncmap

  call conserv_interp(ldt_gds,lb,f,lo,ldt1d,cmap_struc(n)%mi,mo,&
       LDT_domain(n)%lat, LDT_domain(n)%lon, cmap_struc(n)%w112,&
       cmap_struc(n)%w122,cmap_struc(n)%w212,cmap_struc(n)%w222,&
       cmap_struc(n)%n112,cmap_struc(n)%n122,cmap_struc(n)%n212,&
       cmap_struc(n)%n222,LDT_rc%udef,iret)

!--------------------------------------------------------------------    
! assign the interpolated precip field to the output variable
!--------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = ldt1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_cmap

