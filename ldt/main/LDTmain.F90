!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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

!EOP
  implicit none

  character*100 :: configfile
  integer       :: i 
  integer       :: iargc
!BOC  

  i = iargc()
  if(i.ne.1) then 
     print*, 'Usage: '
     print*, 'LDT <ldtconfigfile>'
     stop
  endif
  
  call getarg(1, configfile)
  call LDT_configinit(configfile)
  call LDTinit(trim(LDT_rc%runmode)//char(0))
  call LDTrun(trim(LDT_rc%runmode)//char(0))
!EOC
end program LDTmain


