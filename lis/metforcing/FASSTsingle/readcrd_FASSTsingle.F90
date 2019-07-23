!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_FASSTsingle
! \label{readcrd_FASSTsingle}
!
! !REVISION HISTORY:
! 13 Apr 2007: Bailing Li, Initial Code
! 05 Oct 2010: David Mocko, Updated for FASST single-point test case
!
! !INTERFACE:    
subroutine readcrd_FASSTsingle()
  ! !USES:
  use ESMF
  use LIS_logMod, only : LIS_logunit
  use LIS_coreMod, only : LIS_rc, LIS_config
  use FASSTsingle_forcingMod, only : FASSTsingle_struc
  !
  ! !DESCRIPTION:
  !  This routine reads the options specific to the FASST single-point
  !  test case station data forcing from the LIS configuration file.
  !
  !EOP

  implicit none
  integer :: n, rc

  !      write(LIS_logunit,*) 'starting readcrd_FASSTsingle'
  call ESMF_ConfigFindLabel(LIS_config,"FASST forcing file:",      &
       rc=rc)
  do n = 1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,                      &
          FASSTsingle_struc(n)%FASSTsinglefile,rc=rc)
  enddo
  do n = 1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using FASSTsingle forcing'
     write(LIS_logunit,*) 'FASSTsingle forcing file: ',             &
          trim(FASSTsingle_struc(n)%FASSTsinglefile)
     FASSTsingle_struc(n)%FASSTsingletime2 = 3000.0
     FASSTsingle_struc(n)%FASSTsingletime1 = 0.0
     FASSTsingle_struc(n)%startRead = .false.
  enddo

end subroutine readcrd_FASSTsingle

