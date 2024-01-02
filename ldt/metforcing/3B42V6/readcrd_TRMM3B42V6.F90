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
! !ROUTINE: readcrd_TRMM3B42V6
! \label{readcrd_TRMM3B42V6}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 25 Aug 2006; Yudong Tian, Modification for 3B42, LDT 4.2 release 
! 21 Jun 2013: Soni Yatheendradas; change from earlier code
!
! !INTERFACE:
subroutine readcrd_TRMM3B42V6()

! !USES:
  use ESMF
  use TRMM3B42V6_forcingMod, only : TRMM3B42V6_struc
  use LDT_coreMod, only : LDT_config,LDT_rc
  use LDT_logMod,  only : LDT_logunit

!
! !DESCRIPTION:
!
!  This routine reads the options specific to TRMM 3B42V6 forcing from
!  the LDT configuration file.
!
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LDT_config, "TRMM 3B42V6 forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,TRMM3B42V6_struc(n)%TRMM3B42V6dir,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*)'Using TRMM 3B42V6 forcing (nest):',n
     write(LDT_logunit,*) 'TRMM 3B42V6 forcing directory :', &
           trim(TRMM3B42V6_struc(n)%TRMM3B42V6DIR)
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure
! data is read in during first time step
!------------------------------------------------------------------------
     !TRMM3B42V6_struc(n)%TRMM3B42V6time = 0.0 ! SY
  enddo

end subroutine readcrd_TRMM3B42V6


