!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_capa
! \label{readcrd_capa}
!
! !REVISION HISTORY:
! 23 Nov 2009; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_capa()
! !USES:
  use ESMF
  use capa_forcingMod, only : capa_struc
  use LIS_coreMod,     only : LIS_config,LIS_rc
  use LIS_logMod,      only : LIS_logunit

  implicit none
!
! !DESCRIPTION:
!
!  This routine reads the options specific to CAPA forcing from 
!  the LIS configuration file. 
!  
!EOP
  
  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config, "CAPA forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,capa_struc(n)%capadir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)'Using CAPA forcing'
     write(LIS_logunit,*) 'CAPA forcing directory :',capa_struc(n)%CAPADIR
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
     capa_struc(n)%capatime = 0.0
  enddo

end subroutine readcrd_capa
