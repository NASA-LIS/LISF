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
!
! !ROUTINE: setup_MODISlai
! \label{setup_MODISlai}
!
! !REVISION HISTORY:
!  16 Jul 2008: Sujay Kumar; Initial Specification
!  17 Oct 2011: Yudong Tian; Modified from gfrac to lai 
!
! !INTERFACE:
subroutine setup_MODISlai(n)
! !USES:
  use LIS_vegDataMod,  only : LIS_lai

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n

! !DESCRIPTION:
!  This subroutine sets up the variables required for reading the LAI 
!  climatology data from MODIS
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!EOP      

! This is a monthly climatology. Setting the interval type

  LIS_lai(n)%laiIntervalType = "monthly"

  call ESMF_ConfigFindLabel(LIS_config,&
       "MODIS LAI data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_lai(n)%laifile,rc=rc)
  call LIS_verify(rc,'MODIS LAI data directory: not defined')

end subroutine setup_MODISlai
