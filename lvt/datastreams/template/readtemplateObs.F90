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
! !ROUTINE: readtemplateObs
! \label{readtemplateObs}
!
! !INTERFACE: 
subroutine readtemplateObs(i)
! 
! !USES:   
  use ESMF

  implicit none
  integer, intent(in) :: i 
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine should provide the data reader for the specific station data. 
! This routine is called by LVT at each timestep and is expected to return the 
! observations corresponding to the current time. The variables are returned 
! by simply mapping to the observational data structures of LVT. Two strategies
! can be followed: (1) Read the entire datan and keep it in memory at the first 
! timestep and then at future times, simply index into the data (2) Read the
! data corresponding to the current time - this approach works better if the 
! observations are organized as timestamped files. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  16 Feb 2008: Sujay Kumar, Initial Specification
! 
!EOP

! read the data and then use the log-method (see below) to map the relevant variable(s)
! to the LVT data structures. 

!  call LVT_logSingleVar(LVT_obsData(k)%snowdepth_obs,snowdepth)
  
end subroutine readtemplateObs
