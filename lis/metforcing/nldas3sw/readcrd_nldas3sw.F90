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
! !ROUTINE: readcrd_nldas3sw
!  \label{readcrd_nldas3sw}
!
! !REVISION HISTORY:
! 27 Dec 2024: David Mocko, Initial Specification
!                           (derived from readcrd_nldas20.F90)
!
! !INTERFACE:
subroutine readcrd_nldas3sw()
! !USES:
  use ESMF
  use LIS_coreMod, only         : LIS_rc,LIS_config
  use LIS_logMod,  only         : LIS_logunit,LIS_verify,LIS_endrun
  use nldas3sw_forcingMod, only : nldas3sw_struc

  implicit none
!
! !DESCRIPTION:
!  This routine reads the options specific to NLDAS-3 SWdown
!  forcing from the LIS configuration file.
!
!EOP
  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,                                &
                           "CERES SWdown forcing directory:",rc=rc)
  call LIS_verify(rc,"CERES SWdown forcing directory: not defined")

  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                          &
                                 nldas3sw_struc(n)%nldas3swfordir,rc=rc)
  enddo

  write(unit=LIS_logunit,fmt=*)                                        &
       "[INFO] Using CERES 4-km SWdown data for NLDAS-3"

  do n = 1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*)                                     &
          "[INFO] CERES SWdown forcing directory: ",                   &
          trim(nldas3sw_struc(n)%nldas3swfordir)

     nldas3sw_struc(n)%nldas3swtime1 = 3000.0
     nldas3sw_struc(n)%nldas3swtime2 = 0.0
  enddo

end subroutine readcrd_nldas3sw

