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
! !ROUTINE: readcrd_stg4
! \label{readcrd_stg4}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 05 Jun 2006; Kristi Arsenault, Code and data implementation
!
! !INTERFACE:    
subroutine readcrd_stg4()

! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_config, LDT_rc
  use LDT_logMod,  only : LDT_logunit
  use stg4_forcingMod, only : stg4_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to STAGE4 forcing from 
!  the LDT configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

! - Retrieve Stage IV Forcing Dataset Directory Name Location
    call ESMF_ConfigFindLabel(LDT_config, "STAGE4 forcing directory:",rc=rc)

    do n=1,LDT_rc%nnest   ! Loop over different nested LDT_domains

       call ESMF_ConfigGetAttribute(LDT_config, stg4_struc(n)%stg4dir,rc=rc)
       write(LDT_logunit,*) 'STAGEIV forcing directory :', &
             trim(stg4_struc(n)%STG4DIR)

    !- Setting observed precip times to zero to ensure data is read in
    !   at first time step
       stg4_struc(n)%stg4time = 0.0

    enddo

end subroutine readcrd_stg4
