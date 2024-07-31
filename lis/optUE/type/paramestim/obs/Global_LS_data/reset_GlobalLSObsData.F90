!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: reset_GlobalLSObsdata
! \label{reset_GlobalLSObsdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine reset_GlobalLSObsdata()
! !USES: 
  use GlobalLSDataMod, only : globallsobs_struc

  implicit none
! !ARGUMENTS: 

!
! !DESCRIPTION:
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP

  integer             :: n

  n = 1


  globallsobs_struc(n)%flag = .false. 

end subroutine reset_GlobalLSObsdata

  
  

  

