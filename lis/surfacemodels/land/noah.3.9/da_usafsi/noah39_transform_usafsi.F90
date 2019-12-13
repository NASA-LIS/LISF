!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah39_transform_ldtsi
! \label{noah39_transform_ldtsi}
!
! !REVISION HISTORY:
! 25Jun2006: Sujay Kumar: Initial Specification
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!  09 Apr 2019: Eric Kemp; Modified for Noah 3.9 and LDT-SI
!
! !INTERFACE:
subroutine noah39_transform_ldtsi(n,OBS_State)

! !USES:
  use ESMF
  use LIS_logMod,  only : LIS_verify
  use noah39_lsmMod
!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the LDTSI state
!  (meters) to the lsm state
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP

  ! Since LDT-SI is already in meters, no work is needed here.

end subroutine noah39_transform_ldtsi
