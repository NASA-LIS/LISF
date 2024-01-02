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
! !ROUTINE: readcrd_ecmwf
! \label{readcrd_ecmwf}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_ecmwf(findex)
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod,  only : LDT_logunit
  use ecmwf_forcingMod, only : ecmwf_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to ECMWF forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  integer, intent(in) :: findex

  integer :: n, rc

  if( LDT_rc%met_ecor_parms(findex) .ne. "none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS23R4 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs23r4(n),rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS25R1 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs25r1(n),rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS30R1 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs30r1(n),rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS33R1 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs33r1(n),rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS35R2 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs35r2(n),rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS35R3 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs35r3(n),rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS36R1 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs36r1(n),rc=rc)
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"ECMWF IFS37R2 elevation map:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%elevfileifs37r2(n),rc=rc)
    enddo

    return

  endif

  call ESMF_ConfigFindLabel(LDT_config,"ECMWF forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,ecmwf_struc(n)%ecmwfdir,rc=rc)
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) 'ECMWF forcing directory : ',trim(ecmwf_struc(n)%ecmwfdir)
     ecmwf_struc(n)%ecmwftime1 = 3000.0
     ecmwf_struc(n)%ecmwftime2 = 0.0
  enddo

end subroutine readcrd_ecmwf
