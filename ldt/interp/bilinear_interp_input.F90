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
! !ROUTINE: bilinear_interp_input
! \label{bilinear_interp_input}
!
! !INTERFACE:    
subroutine bilinear_interp_input (nest, gridDesci,&
     n11,n12,n21,n22,w11,w12,w21,w22)

  use LDT_coreMod
  implicit none
! !ARGUMENTS:
  integer, intent(in) :: nest
  real, intent(in)    :: gridDesci(20)
  integer             :: n11(LDT_rc%lnc(nest)*LDT_rc%lnr(nest)),&
       n12(LDT_rc%lnc(nest)*LDT_rc%lnr(nest)),&
       n21(LDT_rc%lnc(nest)*LDT_rc%lnr(nest)),&
       n22(LDT_rc%lnc(nest)*LDT_rc%lnr(nest))
  real                :: w11(LDT_rc%lnc(nest)*LDT_rc%lnr(nest)),&
       w12(LDT_rc%lnc(nest)*LDT_rc%lnr(nest)),&
       w21(LDT_rc%lnc(nest)*LDT_rc%lnr(nest)),&
       w22(LDT_rc%lnc(nest)*LDT_rc%lnr(nest))
!
! !DESCRIPTION: 
!  This subprogram performs issues calls to compute the 
!  interpolation weights and neighbor information for bilinear 
!  interpolation,from any grid to any grid for scalar fields. 
!  The grids are defined by their grid description arrays. 
!  
!  The grid description arrays are based on the decoding 
!  schemes used by NCEP. However, in order to remove the integer
!  arithmetic employed in the original ipolates, the routines
!  are rewritten using real number manipulations. The general 
!  structure remains the same. 
!    
!  The current code recognizes the following projections: \newline
!             (gridDesc(1)=0) equidistant cylindrical \newline
!             (gridDesc(1)=1) mercator cylindrical \newline
!             (gridDesc(1)=3) lambert conformal conical \newline
!             (gridDesc(1)=4) gaussian cylindrical (spectral native) \newline
!             (gridDesc(1)=5) polar stereographic azimuthal \newline
!  where gridDesc could be defined for either the input grid or the 
!  output grid. 
!
!  The arguments are: 
!  \begin{description}
!    \item[gridDesci]
!     input grid description parameters 
!    \item[gridDesco]
!     output grid description parameters 
!    \item[npts] 
!     number of points to in the output field
!    \item[rlat]    
!     output latitudes in degrees
!    \item[rlon]    
!     output longitudes in degrees
!    \item[w11,w12,w21,w22]    
!     weights to be used for interpolation
!    \item[n11,n12,n21,n22]    
!     index of neighbor points 
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[compute\_earth\_coord](\ref{compute_earth_coord}) \newline
!     Computes the earth coordinates for the output grid
!   \item[compute\_grid\_coord](\ref{compute_grid_coord}) \newline
!     Computes the grid coordinates of the input grid, based
!     on the earth coordinates of the output grid. 
!   \item[get\_fieldpos](\ref{get_fieldpos}) \newline
!     computes the field position for a given point
!  \end{description}
!EOP
  integer             :: nc,nr
  real                :: rlat(LDT_rc%lnc(nest)*LDT_rc%lnr(nest))
  real                :: rlon(LDT_rc%lnc(nest)*LDT_rc%lnr(nest))

  integer             :: n
  integer             :: mo, nv 
  real, parameter     :: fill = -9999.0
  real                :: xpts(LDT_rc%lnc(nest)*LDT_rc%lnr(nest)), &
       ypts(LDT_rc%lnc(nest)*LDT_rc%lnr(nest))

  integer             :: i1, i2, j1, j2
  real                :: xi, xf, yi, yf
  integer             :: get_fieldpos

  nc = LDT_rc%lnc(nest)
  nr = LDT_rc%lnr(nest)
  rlat = LDT_domain(nest)%lat
  rlon = LDT_domain(nest)%lon
  mo = nc*nr
  !------------------------------------------------------------------------
  !  Calls the routines to decode the grid description and 
  !  calculates the weights and neighbor information to perform
  !  spatial interpolation. This routine eliminates the need to 
  !  compute these weights repeatedly during interpolation. 
  !------------------------------------------------------------------------
!  if(gridDesco(1).ge.0) then
!     call compute_earth_coord(gridDesco, mo,fill,xpts,ypts,rlon,rlat,nv,.true.)
!  endif

  call compute_grid_coord(gridDesci,mo,fill,xpts,ypts,rlon,rlat,nv)
  do n=1,mo
     xi=xpts(n)
     yi=ypts(n)
     if(xi.ne.fill.and.yi.ne.fill) then
        i1=xi
        i2=i1+1
        j1=yi
        j2=j1+1 
        xf=xi-i1
        yf=yi-j1
        n11(n)=get_fieldpos(i1,j1,gridDesci)
        n21(n)=get_fieldpos(i2,j1,gridDesci)
        n12(n)=get_fieldpos(i1,j2,gridDesci)
        n22(n)=get_fieldpos(i2,j2,gridDesci)
        if(min(n11(n),n21(n),n12(n),n22(n)).gt.0) then
           w11(n)=(1-xf)*(1-yf)
           w21(n)=xf*(1-yf)
           w12(n)=(1-xf)*yf
           w22(n)=xf*yf
        else
           n11(n)=0
           n21(n)=0
           n12(n)=0
           n22(n)=0
        endif
     else
        n11(n)=0
        n21(n)=0
        n12(n)=0
        n22(n)=0
     endif
  enddo
end subroutine bilinear_interp_input

