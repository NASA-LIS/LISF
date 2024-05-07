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
! !ROUTINE: readconfig_WRFoutv2
!  \label{readconfig_WRFoutv2}
!
! !REVISION HISTORY:
! 14 Mar 2013; Sujay Kumar, Initial Code
! 20 Nov 2020; K.R. Arsenault, Updated for different WRF output files
!
! !INTERFACE:    
subroutine readconfig_WRFoutv2(findex)
! !USES:
  use ESMF
  use LDT_coreMod,       only : LDT_rc, LDT_config
  use LDT_logMod,        only : LDT_logunit, LDT_verify
  use WRFoutv2_forcingMod, only : WRFoutv2_struc

! !DESCRIPTION:
!
!  This routine reads the options specific to WRF output forcing from 
!  the LDT configuration file. 
!  
!EOP
  implicit none
  integer, intent(in) :: findex
  integer :: n,rc

  if( LDT_rc%met_ecor_parms(findex) .ne."none" .and. &
      LDT_rc%runmode == "LSM parameter processing" ) then

     call ESMF_ConfigFindLabel(LDT_config,"WRFoutv2 terrain height map:",rc=rc)
     call LDT_verify(rc,"WRFoutv2 terrain height map: not defined")
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,wrfoutv2_struc(n)%file_wrfelev,rc=rc)
     enddo
     return
  endif
  
  call ESMF_ConfigFindLabel(LDT_config,"WRF output v2 forcing directory:",rc=rc)
  call LDT_verify(rc, 'WRF output v2 forcing directory: not defined')
  do n=1,LDT_rc%nnest    
     call ESMF_ConfigGetAttribute(LDT_config,WRFoutv2_struc(n)%WRFoutv2dir,rc=rc)
  enddo

  write(unit=LDT_logunit,fmt=*)'[INFO] Using WRF output v2 forcing'

  do n=1,LDT_rc%nnest
     write(unit=LDT_logunit,fmt=*) '[INFO] WRF output v2 forcing directory :',&
          WRFoutv2_struc(n)%WRFoutv2dir
     WRFoutv2_struc(n)%WRFouttime1 = 3000.0
     WRFoutv2_struc(n)%WRFouttime2 = 0.0
  enddo

end subroutine readconfig_WRFoutv2
