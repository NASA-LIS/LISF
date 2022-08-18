!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_gdasbc
!  \label{readcrd_gdasbc}
!
! !REVISION HISTORY:
! 02Feb2004; Sujay Kumar, Initial Code
! 24 Aug 2007: Chuck Alonge; Modified for use with GDASBC data
! 22 Jan 2012: K. Arsenault; Accommodate GES DISC, NCEP filename conventions
! 14 Mar 2014: David Mocko: Removed elevation file line in config file,
!                           as this data is now in the LDT input file.
!
! !INTERFACE:    
subroutine readcrd_gdasbc()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use gdasbc_forcingMod, only : gdasbc_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to NLDAS-2 forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"GDAS bias corrected forcing directory:",rc=rc)
  call LIS_verify(rc, 'GDASBC forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,gdasbc_struc(n)%gdasbcdir,rc=rc)
  enddo

  write(unit=LIS_logunit,fmt=*)'[INFO] Using GDAS bias corrected forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] GDAS bias corrected forcing directory : ',trim(gdasbc_struc(n)%gdasbcdir)

     gdasbc_struc(n)%ncold = 3600
     gdasbc_struc(n)%nrold = 1500
     gdasbc_struc(n)%gdasbctime1 = 3000.0
     gdasbc_struc(n)%gdasbctime2 = 0.0
  enddo

end subroutine readcrd_gdasbc
