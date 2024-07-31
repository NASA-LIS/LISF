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
! !ROUTINE: get_field_pos
!  \label{get_field_pos}
!
! !REVISION HISTORY:
!   04-10-96  Mark Iredell; Initial Specification
!   03-11-96  Mark Iredell; Allowed hemispheric grids to wrap over one pole
!   05-27-04  Sujay Kumar; Modified code with floating point arithmetic
!   01-17-09  Sujay Kumar; Specified the dlon value separately based on the
!             input grid projection
!
! !INTERFACE:
function get_fieldpos(i,j,gridDesc) result(field_pos)

  implicit none
! !ARGUMENTS: 

  integer     :: field_pos
  real        ::  gridDesc(20)
  integer     ::  i,j

! !DESCRIPTION: 
!  This subprogram returns the field position for a given grid point
!  based on the input grid definition.
!
!  The arguments are: 
!  \begin{description}
!    \item[i]
!     integer x grid point
!    \item[j]
!     integer y grid point
!    \item[gridDesc] 
!     grid description parameters 
!    \item[field\_pos]    
!     integer position in grid field to locate grid point
!    \end{description}
!EOP
 
  integer :: im, jm
  integer :: kscan
  integer :: is1
  integer :: nscan
  integer :: ii, jj
  real    :: rlat1,rlat2
  real    :: rlon1,rlon2
  integer :: iscan
  real    :: dlon
  real    :: dlat
  integer :: ig,jg
  
!  GET GRID DIMENSIONS
  im=gridDesc(2)
  jm=gridDesc(3)
  kscan=0
  is1=0
  nscan=mod(nint(gridDesc(20))/32,2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  ACCOUNT FOR WRAPAROUNDS IN EITHER DIRECTION
  ii=i
  jj=j
  if(gridDesc(1).eq.0) then 
     dlon = gridDesc(9)
  elseif(gridDesc(1).eq.1) then 
     dlon = gridDesc(8)
  elseif(gridDesc(1).eq.4) then 
     dlon = gridDesc(9)
  elseif(gridDesc(1).eq.9) then 
     dlon = gridDesc(10)
  endif

  if(gridDesc(1).eq.0.or.&
       gridDesc(1).eq.1.or.&
       gridDesc(1).eq.4.or.&
       gridDesc(1).eq.9) then
     rlon1=gridDesc(5)
     rlon2=gridDesc(8)
     iscan=mod(nint(gridDesc(20))/128,2)
     ig=nint(360/abs(dlon))
     if(im.ge.ig) then
        ii=mod(i-1+ig,ig)+1
        if((j.le.0.or.j.ge.jm+1).and.mod(ig,2).eq.0) then
           if(gridDesc(1).eq.0) then
              rlat1=gridDesc(4)
              rlat2=gridDesc(7)
              dlat=abs(rlat2-rlat1)/(jm-1)
              if(j.le.0.and.abs(rlat1).gt.90-0.25*dlat) then
                 jj=2-j
                 ii=mod(ii-1+ig/2,ig)+1
              elseif(j.le.0.and.abs(rlat1).gt.90-0.75*dlat) then
                 jj=1-j
                 ii=mod(ii-1+ig/2,ig)+1
              elseif(j.ge.jm+1.and.abs(rlat2).gt.90-0.25*dlat) then
                 jj=2*jm-j
                 ii=mod(ii-1+ig/2,ig)+1
              elseif(j.ge.jm+1.and.abs(rlat2).gt.90-0.75*dlat) then
                 jj=2*jm+1-j
                 ii=mod(ii-1+ig/2,ig)+1
              endif
           elseif(gridDesc(1).eq.4) then
              jg=gridDesc(10)*2
              if(j.le.0.and.jm.eq.jg) then
                 jj=1-j
                 ii=mod(ii-1+ig/2,ig)+1
              elseif(j.ge.jm+1.and.jm.eq.jg) then
                 jj=2*jm+1-j
                 ii=mod(ii-1+ig/2,ig)+1
              endif
           endif
        endif
     endif
  endif
  if(ii.ge.1.and.ii.le.im.and.jj.ge.1.and.jj.le.jm) then
     if(nscan.eq.0) then
        field_pos=ii+(jj-1)*im
     else
        field_pos=jj+(ii-1)*jm
     endif
  else
     field_pos=0
  endif

end function get_fieldpos
