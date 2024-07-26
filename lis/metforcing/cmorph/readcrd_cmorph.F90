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
! !ROUTINE: readcrd_cmorph
! \label{readcrd_cmorph}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 29 Dec 2003; Luis Gustavo, adapted for CMORPH
! 06 Jan 2005; Yudong Tian, modified for LISv4.2
!
! !INTERFACE:    
subroutine readcrd_cmorph()
! !USES:
  use ESMF 
  use LIS_coreMod, only : LIS_config, LIS_rc
  use LIS_logMod, only : LIS_logunit
  use cmorph_forcingMod, only : cmorph_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to CMORPH forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"CMORPH forcing directory:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,cmorph_struc(n)%cmorphdir,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CMORPH domain x-dimension size:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,cmorph_struc(n)%ncold,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CMORPH domain y-dimension size:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,cmorph_struc(n)%nrold,rc=rc)
     
     write(LIS_logunit,*)'Using CMORPH forcing'
     write(LIS_logunit,*) 'CMORPH forcing directory :',trim(cmorph_struc(n)%CMORPHDIR)
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
     cmorph_struc(n)%cmorphtime = 0.0
  enddo
end subroutine readcrd_cmorph
