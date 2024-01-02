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
! !ROUTINE: readcrd_gdasT1534
! \label{readcrd_gdasT1534}
!
!
! !REVISION HISTORY:
!  20 June 2014: Sujay Kumar; initial implementation
!
! !INTERFACE:    
subroutine readcrd_gdasT1534()
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_config
  use LIS_logMod
  use gdasT1534_forcingMod, only: gdasT1534_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GDAST1534 forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"GDAS T1534 forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gdasT1534_struc(n)%gdasT1534dir,rc=rc)
     call LIS_verify(rc,'GDAS T1534 forcing directory: not specified')
  enddo

  do n=1,LIS_rc%nnest
     gdasT1534_struc(n)%nmif = 16 
     write(LIS_logunit,*) 'Using GDAS T1534 forcing'
     write(LIS_logunit,*) 'GDAS T1534 forcing directory :',trim(gdasT1534_struc(n)%GDAST1534DIR)
     gdasT1534_struc(n)%GDAST1534TIME1  = 3000.0
     gdasT1534_struc(n)%GDAST1534TIME2  = 0.0
  enddo

end subroutine readcrd_gdasT1534
