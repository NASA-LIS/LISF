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
! !ROUTINE: readcrd_PALSmetdata
!  \label{readcrd_PALSmetdata}
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
!
!
! !INTERFACE:    
subroutine readcrd_PALSmetdata()
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_config
  use LIS_logMod,  only : LIS_logunit, LIS_verify
  use PALSmetdata_forcingMod, only : PALSmetdata_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to PALS met forcing from 
!  the LIS configuration file. 
!  
!EOP

  implicit none
  integer :: n,rc

  
  call ESMF_ConfigFindLabel(LIS_config,"PALS met forcing directory:",rc=rc)
  call LIS_verify(rc, 'PALS met forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,&
          PALSmetdata_struc(n)%PALSmetdatadir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"PALS met forcing station name:",rc=rc)
  call LIS_verify(rc, 'PALS met forcing station name: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,&
          PALSmetdata_struc(n)%stn_name,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "PALS met forcing data start year:",&
       rc=rc)
  call LIS_verify(rc, 'PALS met forcing data start year: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          PALSmetdata_struc(n)%syr, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "PALS met forcing data start month:",&
       rc=rc)
  call LIS_verify(rc, 'PALS met forcing data start month: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          PALSmetdata_struc(n)%smo, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "PALS met forcing data start day:",&
       rc=rc)
  call LIS_verify(rc, 'PALS met forcing data start day: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          PALSmetdata_struc(n)%sda, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "PALS met forcing data start hour:",&
       rc=rc)
  call LIS_verify(rc, 'PALS met forcing data start hour: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          PALSmetdata_struc(n)%shr, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "PALS met forcing data start minute:",&
       rc=rc)
  call LIS_verify(rc, 'PALS met forcing data start minute: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          PALSmetdata_struc(n)%smn, rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "PALS met forcing data start second:",&
       rc=rc)
  call LIS_verify(rc, 'PALS met forcing data start second: not defined')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          PALSmetdata_struc(n)%sss, rc=rc)
  enddo

  write(unit=LIS_logunit,fmt=*)'Using PALS station met forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) 'PALS met forcing directory :',&
          trim(PALSmetdata_struc(n)%PALSmetdatadir)

     PALSmetdata_struc(n)%fcsttime1 = 3000.0
     PALSmetdata_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine readcrd_PALSmetdata
