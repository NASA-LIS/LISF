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
! !MODULE: LVT_PRIV_gridMod
! \label(LVT_PRIV_gridMod)
!
! !INTERFACE:
module LVT_PRIV_gridMod 
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  The code in this file provides a description of the grid data structure in LIS
!
!  \subsubsection{Overview}
!  This module contains the grid data structure used in LIS. The data structure 
!  contains the variables specified for a single grid cell. It includes:
!  \begin{description}
!   \item[lat] 
!    latitude of the grid cell
!   \item[lon] 
!    longitude of the grid cell
!   \item[col] 
!     index of the grid point along the East-West grid dimension
!   \item[row] 
!    index of the grid point along the North-South grid dimension
!   \item[elev] 
!    Topological elevation of the grid cell
!   \item[slope] 
!    Topological slope of the grid cell
!   \item[aspect] 
!    Topological aspect of the grid cell
!   \item[curv] 
!    Topological curvature of the grid cell
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  14 Nov 2008: Sujay Kumar; Optimized version of grid representation
! 
!EOP

  implicit none
  public griddec
  type griddec
     real            :: lat    
     real            :: lon    
     integer         :: col   
     integer         :: row  
  end type griddec
end module LVT_PRIV_gridMod
