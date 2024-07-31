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
! !ROUTINE: readcrd_HiMATGMU
! \label{readcrd_HiMATGMU}
!
! !REVISION HISTORY:
! 28 Jul 2017; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_HiMATGMU()

! !USES:
  use ESMF
  use HiMATGMU_forcingMod, only : HiMATGMU_struc
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod, only : LIS_logunit
!
! !DESCRIPTION:
!
!  This routine reads the options specific to the GMU HiMAT forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

! - Retrieve Stage IV Forcing Dataset Directory Name Location
    call ESMF_ConfigFindLabel(LIS_config, "HiMAT GMU forcing directory:",rc=rc)

    do n=1,LIS_rc%nnest   ! Loop over different nested LIS_domains

       call ESMF_ConfigGetAttribute(LIS_config, HiMATGMU_struc(n)%HiMATGMUdir,rc=rc)

       write(LIS_logunit,*) '[INFO] Using HiMAT GMU forcing'
       write(LIS_logunit,*) '[INFO] HiMAT GMU forcing directory :', trim(HiMATGMU_struc(n)%HIMATGMUDIR)

    !- Setting observed precip times to zero to ensure data is read in
    !   at first time step
       HiMATGMU_struc(n)%HiMATGMUtime = 0.0

    enddo

end subroutine readcrd_HiMATGMU
