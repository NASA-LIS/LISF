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
! !ROUTINE: readcrd_gddp
! \label{readcrd_gddp}
!
! !REVISION HISTORY:
! 03 Feb 2022; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_gddp()
! !USES:
  use LIS_logMod
  use LIS_coreMod
  use gddp_forcingMod
  use ESMF
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GDDP forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"GDDP forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gddp_struc(n)%odir,rc=rc)
     call LIS_verify(rc,"GDDP forcing directory: is not defined")
  enddo

    call ESMF_ConfigFindLabel(LIS_config,"GDDP forcing scenario:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gddp_struc(n)%scenario,rc=rc)
     call LIS_verify(rc,"GDDP forcing scenario: is not defined")
     if (trim(gddp_struc(n)%scenario).ne."historical") then
        write(LIS_logunit,*) "[ERR] This GDDP forcing scenario is not supported: ",&
                            trim(gddp_struc(n)%scenario)
        write(LIS_logunit,*) "[ERR] ... Please select: 'historical'"
        write(LIS_logunit,*) "[ERR] Program stopping ..."
        call LIS_endrun()
     endif
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"GDDP reference daily climatology directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gddp_struc(n)%ref_dclimodir,rc=rc)
     call LIS_verify(rc,"GDDP reference daily climatology directory: is not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GDDP reference hourly climatology directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gddp_struc(n)%ref_hclimodir,rc=rc)
     call LIS_verify(rc,"GDDP reference hourly climatology directory: is not defined")
  enddo

  
  write(LIS_logunit,*)'[INFO] Using GDDP forcing'

end subroutine readcrd_gddp
