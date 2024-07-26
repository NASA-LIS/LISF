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
! !ROUTINE: reset_ARMdata
! \label{reset_ARMdata}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine reset_ARMdata(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use ARMdata_module,     only : ARMdata_struc

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
     ARMdata_struc(n)%da = -1
     ARMdata_struc(n)%startFlag = .true.      
  enddo

end subroutine reset_ARMdata
