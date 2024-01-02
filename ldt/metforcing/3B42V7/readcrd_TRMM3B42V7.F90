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
! !ROUTINE: readcrd_TRMM3B42V7
! \label{readcrd_TRMM3B42V7}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 25 Aug 2006; Yudong Tian, Modification for 3B42, LDT 4.2 release 
! 21 Jun 2013: Soni Yatheendradas; change from earlier 3B42V6 code,
!              ported to 3B42V7
!
! !INTERFACE:
subroutine readcrd_TRMM3B42V7()
! !USES:
  use ESMF
  use TRMM3B42V7_forcingMod, only : TRMM3B42V7_struc
  use LDT_coreMod, only : LDT_config,LDT_rc
  use LDT_logMod, only : LDT_logunit

!
! !DESCRIPTION:
!
!  This routine reads the options specific to TRMM 3B42V7 forcing from
!  the LDT configuration file.
!
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LDT_config, "TRMM 3B42V7 forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,TRMM3B42V7_struc(n)%TRMM3B42V7dir,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'TRMM 3B42V7 forcing directory :',&
           trim(TRMM3B42V7_struc(n)%TRMM3B42V7DIR)
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure
! data is read in during first time step
!------------------------------------------------------------------------
     !TRMM3B42V7_struc(n)%TRMM3B42V7time = 0.0 ! SY
  enddo

end subroutine readcrd_TRMM3B42V7


