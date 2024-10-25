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
! !ROUTINE: readcrd_AWAP
! \label{readcrd_AWAP}
!
! !REVISION HISTORY:
!
! 30 Jan 2017: Sujay Kumar, Initial version
!
! !INTERFACE:    
subroutine readcrd_AWAP()

! !USES:
  use ESMF
  use AWAP_forcingMod, only : AWAP_struc
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod, only : LIS_logunit
!
! !DESCRIPTION:
!
!  This routine reads the options specific to STAGE4 forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

    call ESMF_ConfigFindLabel(LIS_config, "AWAP forcing directory:",rc=rc)

    do n=1,LIS_rc%nnest   ! Loop over different nested LIS_domains

       call ESMF_ConfigGetAttribute(LIS_config, AWAP_struc(n)%AWAPdir,rc=rc)

       write(LIS_logunit,*) '[INFO] Using AWAP forcing'
       write(LIS_logunit,*) '[INFO] AWAP forcing directory :', trim(AWAP_struc(n)%AWAPDIR)

    !- Setting observed precip times to zero to ensure data is read in
    !   at first time step
       AWAP_struc(n)%AWAPtime = 0.0

    enddo

end subroutine readcrd_AWAP
