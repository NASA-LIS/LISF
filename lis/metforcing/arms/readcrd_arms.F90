!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_arms
! \label{readcrd_arms}
!
! !REVISION HISTORY:
! 28 Jun 2009; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_arms()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc,LIS_config
  use arms_forcingMod, only : arms_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to ARMS station data 
!  forcing from the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"ARMS forcing file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,arms_struc(n)%armsfile,&
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"ARMS precip forcing file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,arms_struc(n)%armspcpfile,&
          rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"ARMS station metadata file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,arms_struc(n)%mdatafile,rc=rc)
  enddo

  do n=1, LIS_rc%nnest
     write(LIS_logunit,*)'Using ARMS forcing'
     write(LIS_logunit,*) 'ARMS forcing directory :',arms_struc(n)%armsfile
     arms_struc(n)%startRead = .true.

     arms_struc(n)%armstime1 = 3000.0
     arms_struc(n)%armstime2 = 0.0
     
     arms_struc(n)%nstns = 84
  enddo

end subroutine readcrd_arms
