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
! !ROUTINE: reset_LPRM_AMSREsm_obs_data
! \label{reset_LPRM_AMSREsm_obs_data}
!
! !REVISION HISTORY:
!  19 Feb 2013: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine reset_LPRM_AMSREsm_obs_data(Obj_Space)
! !USES: 
  use ESMF
  use LPRM_AMSREsm_obsMod, only : LPRM_AMSREsm_obs_struc
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod, only : LIS_logunit
  use LIS_timeMgrMod

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
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
  integer :: n

  n=1
       
  call LIS_registerAlarm("LPRM AMSRE soil moisture read alarm",&
       86400.0, 86400.0)
  LPRM_AMSREsm_obs_struc(n)%startMode = .true. 
  
  write(LIS_logunit,*) 'Reset AMSRE_SR lprm_sm observations'
 
end subroutine reset_LPRM_AMSREsm_obs_data

