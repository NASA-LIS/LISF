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
! !ROUTINE: readcrd_nam242
! \label{readcrd_nam242}
!
!
! !REVISION HISTORY:
!     Sep 2012: NOHRSC/NOAA: Initial specification
!
! !INTERFACE:    
subroutine readcrd_nam242()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use LIS_logMod,        only : LIS_logunit
  use nam242_forcingMod, only: nam242_struc

  implicit none
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NAM242 forcing from 
!  the LIS configuration file. 
!  
!EOP
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"NAM242 forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,nam242_struc(n)%namdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using NAM242 forcing'
     write(LIS_logunit,*) 'NAM242 forcing directory :',trim(nam242_struc(n)%NAMDIR)
     nam242_struc(n)%namtime1  = 3000.0
     nam242_struc(n)%namtime2  = 0.0
  enddo

end subroutine readcrd_nam242
