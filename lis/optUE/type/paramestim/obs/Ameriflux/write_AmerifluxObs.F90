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
! !ROUTINE: write_AmerifluxObs
! \label{write_AmerifluxObs}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine write_AmerifluxObs(Obj_Space)
! !USES: 
  use ESMF

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

end subroutine write_AmerifluxObs

!BOP
! 
! !ROUTINE: Amerifluxobs_filename
! \label{Amerifluxobs_filename}
! 
! !INTERFACE: 
subroutine Amerifluxobs_filename(obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)     :: obsname
! 
! !DESCRIPTION: 
!  This method generates a timestamped filename for the processed
!  PBMR observations. 
! 
!EOP

  character(len=12)    :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/PEOBS/'//cdate1//'.1gs4r'
  
end subroutine Amerifluxobs_filename



