!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_ALMIPII
!  \label{readcrd_ALMIPII}
!
! !REVISION HISTORY:
! 02Feb2004; Sujay Kumar, Initial Code
! 04Mar2007; Kristi Arsenault, Implemented EDAS Elevation Correction
!
! !INTERFACE:    
subroutine readcrd_ALMIPII()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only :LIS_logunit, LIS_verify
  use ALMIPII_forcingMod, only : ALMIPII_struc


! !DESCRIPTION:
!
!  This routine reads the options specific to ALMIPII forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"ALMIPII forcing directory:",rc=rc)
  call LIS_verify(rc, 'ALMIPII forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,ALMIPII_struc(n)%dir,rc=rc)
  enddo


  call ESMF_ConfigFindLabel(LIS_config,"ALMIPII filename prefix:",rc=rc)
  call LIS_verify(rc, 'ALMIPII filename prefix: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,ALMIPII_struc(n)%filename_prefix,rc=rc)
  enddo
  write(LIS_logunit,*) 'Using ALMIPII forcing ..'

end subroutine readcrd_ALMIPII
