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
! !ROUTINE: readcrd_stg2
! \label{readcrd_stg2}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 05 Jun 2006; Kristi Arsenault, Code and data implementation
!
! !INTERFACE:    
subroutine readcrd_stg2()

! !USES:
  use ESMF
  use stg2_forcingMod, only : stg2_struc
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod, only : LIS_logunit
!
! !DESCRIPTION:
!
!  This routine reads the options specific to STAGE2 forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

! - Retrieve Stage II Forcing Dataset Directory Name Location
    call ESMF_ConfigFindLabel(LIS_config, "STAGE2 forcing directory:",rc=rc)

    do n=1,LIS_rc%nnest   ! Loop over different nested LIS_domains

       call ESMF_ConfigGetAttribute(LIS_config, stg2_struc(n)%stg2dir,rc=rc)

       write(LIS_logunit,*) 'Using STAGEII forcing'
       write(LIS_logunit,*) 'STAGEII forcing directory :', trim(stg2_struc(n)%STG2DIR)

    !- Setting observed precip times to zero to ensure data is read in
    !   at first time step
       stg2_struc(n)%stg2time = 0.0

    enddo

end subroutine readcrd_stg2
