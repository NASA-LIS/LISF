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
! !ROUTINE: readcrd_geis
!  \label{readcrd_geis}
!
! !REVISION HISTORY:
! 15 Nov 2023; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_geis()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use geis_forcingMod, only : geis_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to the global EIS forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"Global EIS forcing directory:",rc=rc)
  call LIS_verify(rc, 'Global EIS forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,geis_struc(n)%geisdir,rc=rc)
  enddo

  write(unit=LIS_logunit,fmt=*)'[INFO] Using Global EIS forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] Global EIS forcing directory : ',trim(geis_struc(n)%geisdir)

     geis_struc(n)%geistime1 = 3000.0
     geis_struc(n)%geistime2 = 0.0
  enddo

end subroutine readcrd_geis
