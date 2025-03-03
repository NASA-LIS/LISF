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
! !ROUTINE: reset_FLUXNETdata
! \label{reset_FLUXNETdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine reset_FLUXNETdata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use FLUXNETdata_module,     only : FLUXNETdata_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  write the Walnut Gulch PBMR soil moisture data
!  to disk
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP

  integer            :: n 
  
  do n=1,LIS_rc%nnest
     FLUXNETdata_struc(n)%yr = -1
     FLUXNETdata_struc(n)%mo = LIS_rc%mo
  enddo

end subroutine reset_FLUXNETdata
