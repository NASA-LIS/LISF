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
! !ROUTINE: readcrd_cmap
! \label{readcrd_cmap}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_cmap()
! !USES:
  use ESMF 
  use cmap_forcingMod, only : cmap_struc
  use LIS_coreMod, only : LIS_config,LIS_rc
  use LIS_logMod, only : LIS_logunit

!
! !DESCRIPTION:
!
!  This routine reads the options specific to CMAP forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  
  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config, "CMAP forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,cmap_struc(n)%cmapdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)"Using CMAP forcing"
     write(LIS_logunit,*)" CMAP forcing directory: ",trim(cmap_struc(n)%cmapdir)
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
     cmap_struc(n)%cmaptime = 0.0
  enddo

  cmap_struc(:)%ncold = 512
  cmap_struc(:)%nrold = 256

end subroutine readcrd_cmap
