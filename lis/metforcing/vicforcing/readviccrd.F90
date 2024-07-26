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
! !ROUTINE: readviccrd
! \label{readviccrd}
!
!
! !REVISION HISTORY:
! 16 Nov 2011; James Geiger, Initial Code
!
! !INTERFACE:    
subroutine readviccrd()
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_config
  use LIS_logMod,     only : LIS_logunit
  use vic_forcingMod, only : vicforcing_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to VIC formatted forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%vicdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%forcingInterval,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC NC:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%NC,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC NR:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%NR,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing domain lower left lat:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%slat,rc=rc)

  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing domain lower left lon:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%slon,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing domain upper right lat:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%elat,rc=rc)

  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing domain upper right lon:",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%elon,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing domain resolution (dx):",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%dlon,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"VIC forcing domain resolution (dy):",rc=rc)
  do n=1,LIS_rc%nnest
    call ESMF_ConfigGetAttribute(LIS_config,vicforcing_struc(n)%dlat,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using VIC formatted forcing'
     write(LIS_logunit,*) 'VIC forcing directory :',trim(vicforcing_struc(n)%vicdir)
  enddo

end subroutine readviccrd
