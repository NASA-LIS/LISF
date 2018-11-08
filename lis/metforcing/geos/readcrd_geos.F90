!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_geos
! \label{readcrd_geos}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_geos()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod, only : LIS_logunit
  use geos_forcingMod, only : geos_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GEOS forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"GEOS forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,geos_struc(n)%geosdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using GEOS forcing'
     write(LIS_logunit,*) 'GEOS forcing directory :',geos_struc(n)%GEOSDIR
     geos_struc(n)%geostime1 = 3000.0
     geos_struc(n)%geostime2 = 0.0

  enddo
end subroutine readcrd_geos
