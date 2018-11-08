!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LVT) V1.0
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


