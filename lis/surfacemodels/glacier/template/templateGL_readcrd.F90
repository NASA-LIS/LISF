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
! !ROUTINE: templateGL_readcrd
! \label{templateGL_readcrd}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
! !INTERFACE:
subroutine templateGL_readcrd()
! !USES:
    use ESMF
    use LIS_coreMod, only    : LIS_rc , LIS_config
    use LIS_timeMgrMod, only : LIS_parseTimeString
    use LIS_logMod, only     : LIS_logunit , LIS_verify, LIS_endrun
    use templateGL_Mod, only       : templateGL_struc
    use netcdf 
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NoahMP36 model from
!  the LIS configuration file.
!
!EOP
    implicit none

    integer      :: rc 
    integer      :: n, i
    character*10 :: time 
    character*6  :: str_i
    integer :: ios

    write(LIS_logunit, *) "[INFO] Start reading LIS configuration file for template Glacier model "
    call ESMF_ConfigFindLabel(LIS_config, "template glacier model timestep:", rc = rc)
    do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, Time, rc = rc)
        call LIS_verify(rc, "template glacier model timestep: not defined")
        call LIS_parseTimeString(time, templateGL_struc(n)%ts)
    enddo

  end subroutine TemplateGL_readcrd
