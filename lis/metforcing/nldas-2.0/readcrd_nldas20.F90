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
! !ROUTINE: readcrd_nldas20
!  \label{readcrd_nldas20}
!
! !REVISION HISTORY:
! 11 Jul 2024: David Mocko, Initial Specification
!                           (derived from readcrd_nldas2.F90)
!
! !INTERFACE:
subroutine readcrd_nldas20()
! !USES:
  use ESMF
  use LIS_coreMod, only        : LIS_rc,LIS_config
  use LIS_logMod,  only        : LIS_logunit,LIS_verify,LIS_endrun
  use nldas20_forcingMod, only : nldas20_struc

  implicit none
!
! !DESCRIPTION:
!  This routine reads the options specific to NLDAS-2 forcing from
!  the LIS configuration file.
!
!EOP
  integer :: n,rc

  do n = 1,LIS_rc%nnest
     nldas20_struc(n)%model_level_data  = 0
     nldas20_struc(n)%model_level_press = 0
     nldas20_struc(n)%model_pcp_data    = 0
     nldas20_struc(n)%model_dswrf_data  = 0
  enddo

  call ESMF_ConfigFindLabel(LIS_config,                            &
       "NLDAS-2.0 FORA forcing directory:",rc=rc)
  call LIS_verify(rc,                                              &
       "NLDAS-2.0 FORA forcing directory: not defined")
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          nldas20_struc(n)%nldas20foradir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,                            &
       "NLDAS-2.0 use FORB model level data:",rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          nldas20_struc(n)%model_level_data,default=0,rc=rc)
     call LIS_verify(rc,                                           &
          "NLDAS-2.0 use FORB model level data: is not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,                            &
       "NLDAS-2.0 use FORB model-based SWdown:",rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          nldas20_struc(n)%model_dswrf_data,default=0,rc=rc)
     call LIS_verify(rc,                                           &
          "NLDAS-2.0 use FORB model-based SWdown: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,                            &
       "NLDAS-2.0 use FORB model-based precip:",rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          nldas20_struc(n)%model_pcp_data,default=0,rc=rc)
     call LIS_verify(rc,                                           &
          "NLDAS-2.0 use FORB model-based precip: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,                            &
       "NLDAS-2.0 use FORB model-based pressure:",rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          nldas20_struc(n)%model_level_press,default=0,rc=rc)
     call LIS_verify(rc,                                           &
          "NLDAS-2.0 use FORB model-based pressure: not defined")
  enddo

  write(unit=LIS_logunit,fmt=*)                                    &
       "[INFO] Using NLDAS-2.0 v020 netCDF-4 forcing"

  do n = 1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*)                                 &
          "[INFO] NLDAS-2.0 FORA forcing directory: ",             &
          trim(nldas20_struc(n)%nldas20foradir)
     nldas20_struc(n)%ncold = 464
     nldas20_struc(n)%nrold = 224
     nldas20_struc(n)%nldas20time1 = 3000.0
     nldas20_struc(n)%nldas20time2 = 0.0

     if ((nldas20_struc(n)%model_level_data.eq.1).or.              &
          (nldas20_struc(n)%model_dswrf_data.eq.1).or.             &
          (nldas20_struc(n)%model_pcp_data.eq.1).or.               &
          (nldas20_struc(n)%model_level_press.eq.1)) then
        call ESMF_ConfigFindLabel(LIS_config,                      &
             "NLDAS-2.0 FORB forcing directory:",rc=rc)
        call LIS_verify(rc,                                        &
             "NLDAS-2.0 FORB forcing directory: not defined")
        call ESMF_ConfigGetAttribute(LIS_config,                   &
             nldas20_struc(n)%nldas20forbdir,rc=rc)
        write(unit=LIS_logunit,fmt=*)                              &
             "[INFO] NLDAS-2.0 FORB forcing directory: ",          &
             trim(nldas20_struc(n)%nldas20forbdir)
     endif
  enddo

end subroutine readcrd_nldas20

