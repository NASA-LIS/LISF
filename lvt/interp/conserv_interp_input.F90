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
! !ROUTINE: conserv_interp_input
! \label{conserv_interp_input}
!
! !INTERFACE:
  subroutine conserv_interp_input(gridDesci,gridDesco,npts,&
       rlat2,rlon2,n112,n122,n212,n222,w112,w122,w212,w222)
    
    implicit none
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! 
! !ARGUMENTS: 
    real, intent(in)    :: gridDesci(50) 
    real                :: gridDesco(50)
    integer             :: npts
    real                :: rlat2(npts)
    real                :: rlon2(npts)
    integer             :: n112(npts,25),n122(npts,25),&
         n212(npts,25),n222(npts,25)
    real                :: w112(npts,25),w122(npts,25),&
         w212(npts,25),w222(npts,25)

!
! !DESCRIPTION: 
!  This subprogram performs issues calls to compute the 
!  interpolation weights and neighbor information for budget bilinear 
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
!    \item[rlat2]    
!     output latitudes in degrees
!    \item[rlon2]    
!     output longitudes in degrees
!    \item[w112,w122,w212,w222]    
!     weights to be used for interpolation
!    \item[n112,n122,n212,n222]    
!     index of neighbor points 
!    \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[compute\_earth\_coord](\ref{compute_earth_coord} \newline
!     Computes the earth coordinates for the output grid
!   \item[compute\_grid\_coord](\ref{compute_grid_coord} \newline
!     Computes the grid coordinates of the input grid, based
!     on the earth coordinates of the output grid. 
!   \item[get\_field\_pos](\ref{get_field_pos}) \newline
!     computes the field position for a given point
!  \end{description}
!EOP
    real                :: xpts(npts), ypts(npts)
    real                :: xptb(npts), yptb(npts)
    real                :: rlob(npts), rlab(npts)
    real, parameter     :: fill = -9999.0
    integer             :: ipopt(20)
    integer             :: nb1, nb2, mo
    integer             :: i1, i2, j1, j2
    real                :: xi, xf, yi, yf
    integer             :: get_fieldpos
    integer             :: iret,ib,nb,jb,n,nv,lb,wb
    integer             :: nb3,nb4

    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    mo = npts
    iret=0
    ipopt = 0
    ipopt(1) = -1
    ipopt(2) = -1
    if(gridDesco(1).ge.0) then
       call compute_earth_coord(gridDesco,mo,fill,xpts,ypts,rlon2,rlat2,nv)
       if(mo.eq.0) iret=3
    else
       iret=31
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    nb1=ipopt(1)
    if(nb1.eq.-1) nb1=2
    if(iret.eq.0.and.nb1.lt.0.) iret=32
    if(iret.eq.0.and.nb1.ge.20.and.ipopt(2).ne.-1) iret=32
    if(iret.eq.0) then
       nb2=2*nb1+1
       nb3=nb2*nb2
       nb4=nb3
       if(ipopt(2).ne.-1) then
          nb4=ipopt(2)
          do ib=1,nb1
             nb4=nb4+8*ib*ipopt(2+ib)
          enddo
       endif
    else
       nb2=0
       nb3=0
       nb4=0
    endif

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX
    do nb=1,nb3
       !  locate input points and compute their weights
       jb=(nb-1)/nb2-nb1
       ib=nb-(jb+nb1)*nb2-nb1-1
       lb=max(abs(ib),abs(jb))
       wb=1
       if(ipopt(2).ne.-1) wb=ipopt(2+lb)
       if(wb.ne.0) then
          do n=1,mo
             xptb(n)=xpts(n)+ib/real(nb2)
             yptb(n)=ypts(n)+jb/real(nb2)
          enddo
          call compute_earth_coord(gridDesco,mo,fill,xptb,yptb,rlob,rlab,nv)
          call compute_grid_coord(gridDesci,mo,fill,xptb,yptb,rlob,rlab,nv)
          if(iret.eq.0.and.nv.eq.0.and.lb.eq.0) iret=2
          do n=1,mo
             xi=xptb(n)
             yi=yptb(n)
             if(xi.ne.fill.and.yi.ne.fill) then
                i1=xi
                i2=i1+1
                j1=yi
                j2=j1+1
                xf=xi-i1
                yf=yi-j1
                n112(n,nb)=get_fieldpos(i1,j1,gridDesci)
                n212(n,nb)=get_fieldpos(i2,j1,gridDesci)
                n122(n,nb)=get_fieldpos(i1,j2,gridDesci)
                n222(n,nb)=get_fieldpos(i2,j2,gridDesci)
                if(min(n112(n,nb),n212(n,nb),n122(n,nb),n222(n,nb)).gt.0) then
                   w112(n,nb)=(1-xf)*(1-yf)
                   w212(n,nb)=xf*(1-yf)
                   w122(n,nb)=(1-xf)*yf
                   w222(n,nb)=xf*yf
                else
                   n112(n,nb)=0
                   n212(n,nb)=0
                   n122(n,nb)=0
                   n222(n,nb)=0
                endif
             else
                n112(n,nb)=0
                n212(n,nb)=0
                n122(n,nb)=0
                n222(n,nb)=0
             endif
!             print*, 'def ',n,nb,n113(n,nb),n213(n,nb)
          enddo
       endif
    enddo
  end subroutine conserv_interp_input
