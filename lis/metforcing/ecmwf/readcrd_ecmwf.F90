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
subroutine readcrd_ecmwf()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use ecmwf_forcingMod, only : ecmwf_struc
  use LIS_logMod, only : LIS_logunit
! !DESCRIPTION:
!
!  This routine reads the options specific to ECMWF forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"ECMWF forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%ecmwfdir,rc=rc)
  enddo
#if 0
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS23R4 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs23r4,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS25R1 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs25r1,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS30R1 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs30r1,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS33R1 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs33r1,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS35R2 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs35r2,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS35R3 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs35r3,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS36R1 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs36r1,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ECMWF IFS37R2 elevation map:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ecmwf_struc(n)%elevfileifs37r2,rc=rc)
  enddo
#endif

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'Using ECMWF forcing'
     write(LIS_logunit,*) 'ECMWF forcing directory :',trim(ecmwf_struc(n)%ecmwfDIR)
#if 0
     write(LIS_logunit,*) 'ECMWF IFS23R4 elevation map :',ecmwf_struc(n)%elevfileifs23r4
     write(LIS_logunit,*) 'ECMWF IFS25R1 elevation map :',ecmwf_struc(n)%elevfileifs25r1
     write(LIS_logunit,*) 'ECMWF IFS30R1 elevation map :',ecmwf_struc(n)%elevfileifs30r1
     write(LIS_logunit,*) 'ECMWF IFS33R1 elevation map :',ecmwf_struc(n)%elevfileifs33r1
     write(LIS_logunit,*) 'ECMWF IFS35R2 elevation map :',ecmwf_struc(n)%elevfileifs35r2
     write(LIS_logunit,*) 'ECMWF IFS35R3 elevation map :',ecmwf_struc(n)%elevfileifs35r3
     write(LIS_logunit,*) 'ECMWF IFS36R1 elevation map :',ecmwf_struc(n)%elevfileifs36r1
     write(LIS_logunit,*) 'ECMWF IFS37R2 elevation map :',ecmwf_struc(n)%elevfileifs37r2
#endif
     ecmwf_struc(n)%ecmwftime1 = 3000.0
     ecmwf_struc(n)%ecmwftime2 = 0.0
  enddo

end subroutine readcrd_ecmwf
