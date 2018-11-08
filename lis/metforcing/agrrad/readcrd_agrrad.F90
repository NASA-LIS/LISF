!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_agrrad
! \label{readcrd_agrrad}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_agrrad()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use agrrad_forcingMod, only : agrrad_struc
  use LIS_logMod,        only : LIS_logunit

  implicit none
!
! !DESCRIPTION:
!
!  This routine reads the options specific to AGRMET algorithms from 
!  the LIS configuration file. 
!  
!EOP
  integer   :: n 
  integer   :: rc

  write(LIS_logunit,*)'Using AGRRAD forcing'
  call ESMF_ConfigFindLabel(LIS_config, "AGRRAD forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrrad_struc(n)%agrdir,rc=rc)
     write(LIS_logunit,*) 'AGRRAD forcing directory :',&
                          agrrad_struc(n)%agrdir
  enddo

end subroutine readcrd_agrrad
