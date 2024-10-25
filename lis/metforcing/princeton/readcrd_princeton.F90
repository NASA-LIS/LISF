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
! !ROUTINE: readcrd_princeton
! \label{readcrd_princeton}
!
! !REVISION HISTORY:
! 26 Jan 2007; Hiroko Kato, Initial Code adopted from readcrd_princeton.F90
! 25 Jun 2007; Hiroko Kato, upgraded to LISv5.0
! 15 May 2017: Bailing Li; Added changes for reading in version 2.2 data
! 22 Oct 2018: Daniel Sarmiento; Added changes to support version 3 data
!
! !INTERFACE:    
subroutine readcrd_princeton()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_endrun, LIS_verify
  use princeton_forcingMod, only : princeton_struc
! !DESCRIPTION:
!
!  This routine reads the options specific to PRINCETON forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n, rc

  write(LIS_logunit,*)'[INFO] Using PRINCETON forcing'

  call ESMF_ConfigFindLabel(LIS_config,"PRINCETON forcing directory:",rc=rc)
  call LIS_verify(rc, 'PRINCETON forcing directory: not defined')

  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,princeton_struc(n)%princetondir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"PRINCETON forcing version:",rc=rc)
  call LIS_verify(rc, 'PRINCETON forcing version: not defined')

  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,princeton_struc(n)%version,rc=rc)
     if (princeton_struc(n)%version .ne. "2" .AND. princeton_struc(n)%version .ne. "2.2" .AND. princeton_struc(n)%version .ne. "3") then
         write(LIS_logunit,*) "[ERR] The current version of the Princeton LIS reader only"
         write(LIS_logunit,*) "[ERR] supports versions 2 (1 deg), 2.2 (1 deg), or 3 (0.25 deg)."
         write(LIS_logunit,*) "[ERR] Please input a valid version in the configuration file."
         call LIS_endrun()
     endif
     if (princeton_struc(n)%version .eq. "3") then
         write(LIS_logunit,*) "[WARN] The version 3 of the Princeton driver data is not global."
         write(LIS_logunit,*) "[WARN] The southern latitude boundary is at -59.875 deg. LIS has"
         write(LIS_logunit,*) "[WARN] not been tested when the domain defined in LDT extends beyond"
         write(LIS_logunit,*) "[WARN] this boundary. Redefine the domain in LDT if you encounter"
         write(LIS_logunit,*) "[WARN] any issues."
     endif
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'[INFO] PRINCETON forcing directory :',trim(princeton_struc(n)%princetonDIR)
     write(LIS_logunit,*)'[INFO] PRINCETON forcing version   :',princeton_struc(n)%version
     princeton_struc(n)%princetontime1 = 3000.0
     princeton_struc(n)%princetontime2 = 0.0
  enddo

end subroutine readcrd_princeton
