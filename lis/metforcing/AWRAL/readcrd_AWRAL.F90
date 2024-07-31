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
! !ROUTINE: readcrd_AWRAL
! \label{readcrd_AWRAL}
!
! !REVISION HISTORY:
!
! 30 Jan 2017: Sujay Kumar, Initial version
!
! !INTERFACE:    
subroutine readcrd_AWRAL()

! !USES:
  use ESMF
  use AWRAL_forcingMod, only : AWRAL_struc
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

    call ESMF_ConfigFindLabel(LIS_config, "AWRAL forcing directory:",rc=rc)

    do n=1,LIS_rc%nnest   ! Loop over different nested LIS_domains

       call ESMF_ConfigGetAttribute(LIS_config, AWRAL_struc(n)%AWRALdir,rc=rc)

       write(LIS_logunit,*) '[INFO] Using AWRAL forcing'
       write(LIS_logunit,*) '[INFO] AWRAL forcing directory :', trim(AWRAL_struc(n)%AWRALDIR)

    !- Setting observed forcing times to zero to ensure data is read in
    !   at first time step
       AWRAL_struc(n)%AWRALtime = 0.0

    enddo

end subroutine readcrd_AWRAL
