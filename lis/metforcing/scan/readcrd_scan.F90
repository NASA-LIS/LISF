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
! !ROUTINE: readcrd_scan
! \label{readcrd_scan}
!
! !REVISION HISTORY:
! 13 Apr 2007; Bailing Li, Initial Code
!
! !INTERFACE:    
subroutine readcrd_scan()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use scan_forcingMod, only : scan_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to SCAN station data 
!  forcing from the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"SCAN forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,scan_struc(n)%scandir,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"SCAN metadata file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,scan_struc(n)%metadata,&
          rc=rc)
     write(LIS_logunit,*)'Using SCAN forcing'
     write(LIS_logunit,*) 'SCAN forcing directory :',trim(scan_struc(n)%scandir)
     write(LIS_logunit,*) 'SCAN metadata file : ',trim(scan_struc(n)%metadata)
     scan_struc(n)%scantime2  = 3000.0
     scan_struc(n)%scantime1  = 0.0
     scan_struc(n)%startRead = .false.
  enddo

end subroutine readcrd_scan
