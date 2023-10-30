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
! !ROUTINE: LDTmain
!  \label{LDTmain}
!   Main program for Land Data Toolkit 
!  
! !DESCRIPTION: 
!  Main driver program for Land Data Toolkit (LDT). The code first reads
!  the run time configuration file, then issues the calls to performs the 
!  initialization, followed by the run steps. 
!
! !REVISION HISTORY: 
!  31 Mar 2012    Sujay Kumar  Initial Specification
!
program LDTmain
! !USES:       
  use LDT_coreMod, only : LDT_rc, LDT_configinit
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

!EOP
  implicit none

  character(len=LDT_CONST_PATH_LEN) :: configfile
  integer       :: i 
  integer       :: iargc
!BOC  

  i = iargc()
  if(i.ne.1) then 
     print*, 'Usage: '
     print*, 'LDT <ldtconfigfile>'
     error stop 1
  endif
  
  call getarg(1, configfile)
  call LDT_configinit(configfile)
  call LDTinit(trim(LDT_rc%runmode)//char(0))
  call LDTrun(trim(LDT_rc%runmode)//char(0))
!EOC
end program LDTmain


