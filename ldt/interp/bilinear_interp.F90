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
! !ROUTINE: bilinear_interp
!  \label{bilinear_interp}
!        
! !REVISION HISTORY:
!   04-10-96  Mark Iredell; Initial Specification
!   05-27-04  Sujay Kumar : Modified verision with floating point arithmetic 
!
! !INTERFACE:
subroutine bilinear_interp(gridDesco,li,gi,lo,go,mi,mo, & 
     rlat,rlon,w11,w12,w21,w22,n11,n12,n21,n22,udef,iret)
!USES:
!$ use omp_lib
  implicit none
! !ARGUMENTS: 
  real      :: gridDesco(20)
  integer   :: mi
  integer   :: mo
  logical*1 :: li(mi)
  logical*1 :: lo(mo)
  real      :: gi(mi)
  real      :: go(mo)
  real      :: rlat(mo)
  real      :: rlon(mo)
  real      :: w11(mo),w12(mo), w21(mo),w22(mo)
  integer   :: n11(mo),n12(mo),n21(mo),n22(mo)
  real      :: udef
  integer   :: iret
!
! !DESCRIPTION: 
!  This subprogram performs bilinear interpolation
!  from any grid to any grid for scalar fields. The routine is based
!  on the spatial interpolation package ipolates from NCEP. 
!             
!  The algorithm simply computes (weighted) averages
!  of bilinearly interpolated points arranged in a square box
!  centered around each output grid point and stretching
!  nearly halfway to each of the neighboring grid points.
!  the grids are defined by their grid description arrays. 
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
!  output grid. The routine also returns the  
!  the number of output grid points
!  and their latitudes and longitudes are also returned.
!  The input bitmaps will be interpolated to output bitmaps.
!  output bitmaps will also be created when the output grid
!  extends outside of the domain of the input grid.
!  the output field is set to 0 where the output bitmap is off.
! 
!  The arguments are: 
!  \begin{description}
!    \item[gridDesco]
!     output grid description parameters 
!    \item[li]
!     logical input bitmaps
!    \item[gi]
!     real input fields to interpolate
!    \item[lo]
!     logical output bitmaps
!    \item[go]
!     real output fields interpolated
!    \item[mi]
!     integer dimension of input grid fields 
!    \item[mo]
!     integer dimension of output grid fields
!    \item[rlat]    
!     output latitudes in degrees
!    \item[rlon]    
!     output longitudes in degrees
!    \item[w11,w12,w21,w22]    
!     weights to be used for interpolation
!    \item[n11,n12,n21,n22]    
!     index of neighbor points 
!    \item[udef]
!     undefined value to be used
!    \item[iret]
!     return code (0-success)
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[polfixs](\ref{polfixs}) \newline
!    Apply corrections for poles
!  \end{description}
!EOP

  integer   :: n
  integer   :: ibo
  integer   :: nn
  real      :: wo(mo)
  real, parameter :: fill=-9999.

  nn = mo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  INTERPOLATE WITH OR WITHOUT BITMAPS
!$OMP PARALLEL
!$OMP DO
  do n=1, nn
     go(n)=0.
     wo(n)=0.

     if(n11(n).gt.0) then 
        if(li(n11(n))) then
           go(n)=go(n)+w11(n)*gi(n11(n))
           wo(n)=wo(n)+w11(n)
        endif
     endif
     if(n21(n).gt.0) then 
        if(li(n21(n))) then
           go(n)=go(n)+w21(n)*gi(n21(n))
           wo(n)=wo(n)+w21(n)
        endif
     endif
     if(n12(n).gt.0) then 
        if(li(n12(n))) then
           go(n)=go(n)+w12(n)*gi(n12(n))
           wo(n)=wo(n)+w12(n)
        endif
     endif
     if(n22(n).gt.0) then 
        if(li(n22(n))) then
           go(n)=go(n)+w22(n)*gi(n22(n))
           wo(n)=wo(n)+w22(n)
        endif
     endif
  enddo
!$OMP END DO
  ibo=1
!$OMP DO 
  do n=1,nn
     lo(n)=wo(n).ge.0.5
     if(lo(n)) then
        go(n)=go(n)/wo(n) 
     else
        ibo=1
        go(n)=udef
     endif
  enddo
!$OMP END DO 
!$OMP END PARALLEL
  if(gridDesco(1).eq.0) call polfixs(nn,mo,1,rlat,rlon,ibo,lo,go)
  iret = 0 

end subroutine bilinear_interp
