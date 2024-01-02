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
! !ROUTINE: LVTpreproc
!  \label{LVTpreproc}
!  
! !DESCRIPTION: 
!   Program that generates the required input "LIS" file for LVT. This
!   is a simplied version of the Land Data Toolkit (LDT). 
!   For full functionality and for supporting LIS version 7 outputs, 
!   please download and use LDT. 
!
! !REVISION HISTORY: 
!  31 Mar 2012    Sujay Kumar  Initial Specification
!
program LVTpreproc
! !USES:       
  use preprocMod
!EOP
  implicit none

  character*100 :: configfile
  integer       :: i 
  integer       :: iargc

  
  i = iargc()
  if(i.ne.1) then 
     print*, 'Usage: '
     print*, 'LVTpreproc <configfile>'
     stop
  endif
  
  call getarg(1, configfile)
  call readConfigFile(configfile)
  call paramProcInit()
  call paramProcWrite()

!EOC
end program LVTpreproc


