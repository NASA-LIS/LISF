!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_ceop
! \label{readcrd_ceop}
!
! !REVISION HISTORY:
! 08 Dec 2004; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_ceop()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc,LIS_config
  use ceop_forcingMod, only : ceop_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to CEOP station data 
!  forcing from the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config, "CEOP location index:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ceop_struc(n)%location,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CEOP forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ceop_struc(n)%ceopdir,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"CEOP metadata file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,ceop_struc(n)%metadata,&
          rc=rc)
  enddo
  do n=1, LIS_rc%nnest
     write(LIS_logunit,*)'Using CEOP forcing'
     write(LIS_logunit,*) 'CEOP forcing directory :',ceop_struc(n)%ceopdir
     write(LIS_logunit,*) 'CEOP metadata file : ',ceop_struc(n)%metadata
     ceop_struc(n)%ceoptime2  = 3000.0
     ceop_struc(n)%ceoptime1  = 0.0
     ceop_struc(n)%startRead = .false.
  enddo

end subroutine readcrd_ceop
