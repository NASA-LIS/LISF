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
! !ROUTINE: readcrd_narr
! \label{readcrd_narr}
!
!
! !REVISION HISTORY:
! 30 APR 2009; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine readcrd_narr()
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_config
  use LIS_logMod,    only : LIS_logunit, LIS_verify
  use narr_forcingMod, only: narr_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to NARR forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  integer :: n, rc

  call ESMF_ConfigFindLabel(LIS_config,"NARR forcing directory:",rc=rc)
  call LIS_verify(rc,'NARR forcing directory: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,narr_struc(n)%narrdir,rc=rc)
  enddo

  narr_struc(:)%nc = 349
  narr_struc(:)%nr = 277
!  narr_struc(:)%nlevels = 29

#if 0
  call ESMF_ConfigFindLabel(LIS_config,"NARR domain x-dimension size:",rc=rc)
  call LIS_verify(rc,'NARR domain x-dimension size: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,narr_struc(n)%nc,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"NARR domain y-dimension size:",rc=rc)
  call LIS_verify(rc,'NARR domain y-dimension size: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,narr_struc(n)%nr,rc=rc)
  enddo
#endif

  call ESMF_ConfigFindLabel(LIS_config,"NARR domain z-dimension size:",rc=rc)
  call LIS_verify(rc,'NARR domain z-dimension size: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,narr_struc(n)%nlevels,&
          default=29,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) 'Using NARR forcing'
     write(LIS_logunit,*) 'NARR forcing directory :',trim(narr_struc(n)%NARRDIR)
     narr_struc(n)%NARRTIME1  = 3000.0
     narr_struc(n)%NARRTIME2  = 0.0
  enddo

end subroutine readcrd_narr
