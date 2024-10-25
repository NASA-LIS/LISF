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
! !ROUTINE: readcrd_Loobos
! \label{readcrd_Loobos}
!
! !REVISION HISTORY:
! 13 Apr 2007: Bailing Li, Initial Code
! 05 Oct 2010: David Mocko, Updated for Bondville test case
! 20 Feb 2018: Shugong Wang, Updated for Loobos test case 
! 
! !INTERFACE:    
subroutine readcrd_Loobos()
! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use Loobos_forcingMod, only : Loobos_struc
!
! !DESCRIPTION:
!  This routine reads the options specific to the Loobos
!  test case station data forcing from the LIS configuration file.
!
!EOP

  implicit none
  integer :: n, rc

!      write(LIS_logunit,*) 'starting readcrd_Loobos'
  call ESMF_ConfigFindLabel(LIS_config,"Loobos forcing file:", &
       rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          Loobos_struc(n)%Loobosfile,rc=rc)
  enddo
  do n = 1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using Loobos forcing'
     write(LIS_logunit,*) 'Loobos forcing file: ',             &
          trim(Loobos_struc(n)%Loobosfile)
     Loobos_struc(n)%Loobostime2 = 17520
     Loobos_struc(n)%Loobostime1 = 0.0
     Loobos_struc(n)%startRead = .false.
  enddo
  
end subroutine readcrd_Loobos

