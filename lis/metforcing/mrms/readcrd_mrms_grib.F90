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
! !ROUTINE: readcrd_mrms_grib
! \label{readcrd_mrms_grib}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 05 Jun 2006; Kristi Arsenault, Code and data implementation
! 13 Feb 2015; Jonathan Case, Modified for MRMS QPE
! 05 Sep 2017; Jessica Erlingis, Modified for operational MRMS
! 22 Feb 2019; Jessica Erlingis, Modified to add masking option
!
! !INTERFACE:    
subroutine readcrd_mrms_grib()

! !USES:
  use ESMF
  use mrms_grib_forcingMod, only : mrms_grib_struc
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod, only : LIS_logunit, LIS_endrun, LIS_verify
!
! !DESCRIPTION:
!
!  This routine reads the options specific to MRMS forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

! - Retrieve MRMS Forcing Dataset Directory Name Location
  call ESMF_ConfigFindLabel(LIS_config, "MRMS forcing directory:",rc=rc)
  call LIS_verify(rc,'MRMS forcing directory: not defined')

  do n=1,LIS_rc%nnest   ! Loop over different nested LIS_domains

    call ESMF_ConfigGetAttribute(LIS_config, mrms_grib_struc(n)%mrms_grib_dir,rc=rc)

    write(LIS_logunit,*) 'Using MRMS forcing'
    write(LIS_logunit,*) 'MRMS forcing directory: ', trim(mrms_grib_struc(n)%MRMS_GRIB_DIR)

    !- Setting observed precip times to zero to ensure data is read in
    !   at first time step
    mrms_grib_struc(n)%mrms_grib_time = 0.0

  enddo

! call ESMF_ConfigFindLabel(LIS_config,"MRMS masking:",rc=rc)
! call LIS_verify(rc,'MRMS masking: not defined')

!  call ESMF_ConfigFindLabel(LIS_config,"MRMS mask directory:",rc=rc)
!  call LIS_verify(rc,'MRMS mask directory: not defined. If not using a mask, set to None.')

  do n=1,LIS_rc%nnest
    call ESMF_ConfigFindLabel(LIS_config,"MRMS masking:",rc=rc)
    call LIS_verify(rc,'MRMS masking: not defined')
    call ESMF_ConfigGetAttribute(LIS_config,mrms_grib_struc(n)%mrms_mask_opt,rc=rc)
    if ((mrms_grib_struc(n)%mrms_mask_opt.ne.0).AND.(mrms_grib_struc(n)%mrms_mask_opt.ne.1)) then
       write(LIS_logunit,*) "[ERR] Valid options for MRMS masking are 0=No or 1=Yes"
       call LIS_endrun()
    endif
    if (mrms_grib_struc(n)%mrms_mask_opt.eq.1) then
      call ESMF_ConfigFindLabel(LIS_config,"MRMS mask threshold:",rc=rc)
      call LIS_verify(rc,'MRMS mask threshold: not defined. Please specify cutoff for masking.')
      call ESMF_ConfigGetAttribute(LIS_config,mrms_grib_struc(n)%mrms_mask_thresh,rc=rc)
    endif
    call ESMF_ConfigFindLabel(LIS_config,"MRMS mask directory:",rc=rc)
    call LIS_verify(rc,'MRMS mask directory: not defined. If not using a mask, set to None.')
    call ESMF_ConfigGetAttribute(LIS_config,mrms_grib_struc(n)%mrms_mask_dir,rc=rc)
    write(LIS_logunit,*) 'MRMS mask directory: ', mrms_grib_struc(n)%MRMS_MASK_DIR
  end do

end subroutine readcrd_mrms_grib
