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
! !ROUTINE: readconfig_WRF_AKdom
!  \label{readconfig_WRF_AKdom}
!
! !REVISION HISTORY:
!  21 Jun 2021: K.R. Arsenault; Updated for different WRF output files
!
! !INTERFACE:    
subroutine readconfig_WRF_AKdom(findex)
! !USES:
  use ESMF
  use LDT_coreMod,       only : LDT_rc, LDT_config
  use LDT_logMod,        only : LDT_logunit, LDT_verify
  use WRF_AKdom_forcingMod, only : WRFAK_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to WRF AK forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none

  integer, intent(in) :: findex
  integer :: n,rc

  if( LDT_rc%met_ecor_parms(findex) .ne."none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

     call ESMF_ConfigFindLabel(LDT_config,"WRF AK terrain height map:",rc=rc)
     call LDT_verify(rc,"WRF AK terrain height map: not defined")
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,WRFAK_struc(n)%file_wrfelev,rc=rc)
     enddo
     return
  endif

  call ESMF_ConfigFindLabel(LDT_config,"WRF AK forcing directory:",rc=rc)
  call LDT_verify(rc, 'WRF AK forcing directory: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,WRFAK_struc(n)%WRFAKdir,rc=rc)
  enddo

  write(unit=LDT_logunit,fmt=*)'[INFO] Using WRF AK forcing'

  do n=1,LDT_rc%nnest
     write(unit=LDT_logunit,fmt=*) '[INFO] WRF AK forcing directory :',WRFAK_struc(n)%WRFAKdir
     WRFAK_struc(n)%nest_id = 1

     WRFAK_struc(n)%WRFouttime1 = 3000.0
     WRFAK_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine readconfig_WRF_AKdom
