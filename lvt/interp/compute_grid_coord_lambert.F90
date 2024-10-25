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
! !ROUTINE: compute_grid_coord_lambert
!  \label{compute_grid_coord_lambert}
!
! !INTERFACE:
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
!   04-10-96 Mark Iredell;  Initial Specification
!   07-15-05 Sujay Kumar; Modified verision with floating point arithmetic. 
! 
!EOP
!BOP
! 
! !INTERFACE:
subroutine compute_grid_coord_lambert(gridDesc,npts,fill,xpts,ypts,& 
     rlon,rlat,nret)
! !USES: 
  use map_utils

  implicit none
! !ARGUMENTS: 
  real            :: gridDesc(50)
  integer         :: npts
  real            :: fill
  real            :: xpts(npts),ypts(npts)
  real            :: rlat(npts)
  real            :: rlon(npts)
  integer         :: nret

! !DESCRIPTION:
!  This subroutine computes the grid coordinates of 
!  the specified domain for a lambert conformal projection.
!  This routine is based on the grid
!  decoding routines in the NCEP interoplation package. 
!  
!  \begin{description}
!    \item[gridDesc]
!     grid description parameters 
!    \item[npts]
!     integer maximum number of coordinates
!    \item[fill]
!     fill value to set invalid output data
!    \item[xpts]
!     output grid x point coordinates
!    \item[ypts]
!     output grid y point coordinates
!    \item[rlat]    
!     input latitudes in degrees
!    \item[rlon]    
!     input longitudes in degrees
!    \end{description}
!
! !NOTE:
!  This routine is currently unsupported. 
!EOP
  type(proj_info) :: proj
  integer  :: i
  
  if(griddesc(1).eq.3) then
     call map_set(PROJ_LC,gridDesc(4),gridDesc(5),&
          gridDesc(8)*1000,gridDesc(11),gridDesc(10),gridDesc(7),&
          nint(gridDesc(2)),nint(gridDesc(3)),proj)
     do i=1,npts
        call latlon_to_ij(proj,rlat(i), rlon(i), xpts(i),ypts(i))
     enddo
  endif  

end subroutine compute_grid_coord_lambert
