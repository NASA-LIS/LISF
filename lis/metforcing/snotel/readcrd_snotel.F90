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
! !ROUTINE: readcrd_snotel
! \label{readcrd_snotel}
!
! !REVISION HISTORY:
! 08Jun2011; Yuqiong Liu, Initial Specification
!
! !INTERFACE:    
subroutine readcrd_snotel()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use snotel_forcingMod, only : snotel_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to SNOTEL station data 
!  forcing from the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"SNOTEL forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snotel_struc(n)%snoteldir,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"SNOTEL metadata file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,snotel_struc(n)%metadata,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"SNOTEL coord file:",rc=rc)
  do n=1,LIS_rc%nnest
!     call ESMF_ConfigGetAttribute(LIS_config,snotel_struc(n)%metadata,&
!          rc=rc)
     call ESMF_ConfigGetAttribute(LIS_Config, snotel_struc(n)%coorddata,&
          rc=rc)

     write(LIS_logunit,*)'Using SNOTEL forcing'
     write(LIS_logunit,*) 'SNOTEL forcing directory :',trim(snotel_struc(n)%snoteldir)
     write(LIS_logunit,*) 'SNOTEL metadata file : ',trim(snotel_struc(n)%metadata)
     write(LIS_logunit,*) 'SNOTEL coord file : ',trim(snotel_struc(n)%coorddata)

     snotel_struc(n)%snoteltime2  = 3000.0
     snotel_struc(n)%snoteltime1  = 0.0
     snotel_struc(n)%startRead = .false.
  enddo

end subroutine readcrd_snotel
