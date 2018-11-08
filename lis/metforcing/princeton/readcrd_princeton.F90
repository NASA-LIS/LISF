!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
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
!
! !INTERFACE:    
subroutine readcrd_princeton()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit
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
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,princeton_struc(n)%princetondir,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"PRINCETON forcing version:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,princeton_struc(n)%version,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'[INFO] PRINCETON forcing directory :',princeton_struc(n)%princetonDIR
     princeton_struc(n)%princetontime1 = 3000.0
     princeton_struc(n)%princetontime2 = 0.0
  enddo

end subroutine readcrd_princeton
