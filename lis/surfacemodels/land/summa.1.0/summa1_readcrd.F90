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
! !ROUTINE: summa1_readcrd
! \label{summa1_readcrd}
!
! !REVISION HISTORY:
! 21 Jul 2004; Sujay Kumar, Initial Code
! 23 Oct 2007; Kristi Arsenault, Updated for use in LISv5.0
! 23 May 2014: David Mocko, Updated for LIS-7.0
!
! !INTERFACE:    
subroutine summa1_readcrd()

! !USES:
  use ESMF
  use LIS_coreMod,     only : LIS_rc, LIS_config
  use LIS_timeMgrMod, only : LIS_parseTimeString
  use LIS_logMod,      only : LIS_logunit, LIS_verify
  use summa1_lsmMod, only : summa1_struc

!
! !DESCRIPTION:
!
!  This routine reads the options specific to summa1 LSM 
!  option from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: rc
  integer :: n
  character*10  :: time

  call ESMF_ConfigFindLabel(LIS_config,"SUMMA.1.0 model timestep:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'SUMMA.1.0 model timestep: not defined')

     call LIS_parseTimeString(time,summa1_struc(n)%ts)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SUMMA.1.0 master file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          summa1_struc(n)%summaFileManagerFile,rc=rc)
     call LIS_verify(rc,'SUMMA.1.0 master file: not defined')
  enddo

  write(LIS_logunit,*)'[INFO] Running SUMMA LSM Option:'
  do n=1,LIS_rc%nnest
     summa1_struc(n)%summa1open=0
  enddo

end subroutine summa1_readcrd
