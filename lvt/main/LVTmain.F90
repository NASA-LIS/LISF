!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: LVTmain
!  \label{LVTmain}
!   Main program for LIS Verification Toolkit
!
! !INTERFACE:
program LVTmain
! 
! !USES:       
  use LVT_coreMod, only : LVT_rc, LVT_configinit
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  Main driver program for LIS Verification Toolkit. The code first reads
!  the run time configuration file, then issues the calls to performs the 
!  initialization, followed by the run steps. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 08    Sujay Kumar  Initial Specification
! 
!EOP
  implicit none

  character*500 :: configfile = 'lvt.config'
  integer       :: i 
  integer       :: iargc
!BOC  

  i = iargc()
  if(i.ne.1) then 
     print*, 'Usage: '
     print*, 'LVT <lvtconfigfile>'
     error stop 1
  endif
  
  call getarg(1, configfile)
  LVT_rc%configfile = configfile
  call LVT_configinit(configfile)
  call LVTinit(trim(LVT_rc%runmode)//char(0))
  call LVTrun(trim(LVT_rc%runmode)//char(0))
!EOC
end program LVTmain


